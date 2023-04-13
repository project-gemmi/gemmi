.. highlight:: cpp

Math and unit cell
##################

.. _coordinates:

Vectors and coordinates
=======================

Coordinates are represented by two classes:

* ``Position`` for coordinates in Angstroms (orthogonal coordinates),
* ``Fractional`` for coordinates relative to the unit cell
  (fractional coordinates).

Both ``Position`` and ``Fractional`` are derived from ``Vec3``,
which has three numeric properties: ``x``, ``y`` and ``z``.

.. doctest::

    >>> import gemmi
    >>> v = gemmi.Vec3(1.2, 3.4, 5.6)
    >>> v.x, v.z = v.z, v.x

The elements can also be indexed:

.. doctest::

    >>> v[0]  # C++ equivalent: v.at(0)
    5.6

The only reason to have separate types is to prevent functions that
expect fractional coordinates from accepting orthogonal ones, and vice versa.

In Python, vectors can be created from list and exported to a list:

.. doctest::

    >>> v.fromlist([3.0, -4.5, 4.6])
    >>> v.tolist()
    [3.0, -4.5, 4.6]

Vec3 has a number of methods:

.. doctest::

    >>> -v
    <gemmi.Vec3(-3, 4.5, -4.6)>
    >>> v + v
    <gemmi.Vec3(6, -9, 9.2)>
    >>> v += gemmi.Vec3(0, 0, 0)
    >>> v - v
    <gemmi.Vec3(0, 0, 0)>
    >>> v -= gemmi.Vec3(0, 0, 0)
    >>> 2 * v
    <gemmi.Vec3(6, -9, 9.2)>
    >>> v / 2
    <gemmi.Vec3(1.5, -2.25, 2.3)>
    >>> v.dot(v)
    50.41
    >>> v.cross(v)
    <gemmi.Vec3(0, 0, 0)>
    >>> v.length()
    7.1
    >>> v.approx(v, epsilon=1e-9)
    True

These methods are inherited by Position and Fractional.
Some of them are overridden to return the derived type,
others are not overridden and return the base class:

.. doctest::

    >>> frac = gemmi.Fractional(0.5, 0.5, 0.5)
    >>> frac + frac
    <gemmi.Fractional(1, 1, 1)>
    >>> 2 * frac
    <gemmi.Vec3(1, 1, 1)>

Additionally, derived classes have own methods:

.. doctest::

    >>> gemmi.Fractional(0.3, -0.3, 1.5).wrap_to_unit()
    <gemmi.Fractional(0.3, 0.7, 0.5)>
    >>> gemmi.Position(60, 70, 70).dist(gemmi.Position(50, 50, 50))
    30.0

and we have non-member functions that calculate angles in Cartesian
coordinates:

.. doctest::

    >>> from math import degrees
    >>> p1 = gemmi.Position(0, 0, 0)
    >>> p2 = gemmi.Position(0, 0, 1)
    >>> p3 = gemmi.Position(0, 1, 0)
    >>> p4 = gemmi.Position(-1, 1, 0)
    >>> degrees(gemmi.calculate_angle(p1, p2, p3))
    45.00000000000001
    >>> degrees(gemmi.calculate_dihedral(p1, p2, p3, p4))
    90.0

There are more functions in C++.
See headers ``gemmi/math.hpp`` and ``gemmi/calculate.hpp``.


3x3 matrices
============

Gemmi has only 3x3 matrices.
``Mat33`` is a matrix that can represent rotation and scaling,
and ``SMat33`` is a symmetric matrix (6 elements) that can represent
anisotropic ADP tensor. Let us start with the former.

Similarly to vectors, Mat33 can be converted to and from a Python's list:

.. doctest::

  >>> mat33 = gemmi.Mat33()  # identity matrix
  >>> m = mat33.tolist()
  >>> m
  [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
  >>> m[1][2] = -5
  >>> mat33.fromlist(m)
  >>> mat33
  <gemmi.Mat33 [1, 0, 0]
               [0, 1, -5]
               [0, 0, 1]>

We have usual methods expected in a matrix class:

.. doctest::

  >>> mat33.trace()
  3.0
  >>> mat33.transpose()
  <gemmi.Mat33 [1, 0, 0]
               [0, 1, 0]
               [0, -5, 1]>
  >>> mat33 + mat33
  <gemmi.Mat33 [2, 0, 0]
               [0, 2, -10]
               [0, 0, 2]>
  >>> mat33 @ mat33  # in C++: .multiply()
  <gemmi.Mat33 [1, 0, 0]
               [0, 1, -10]
               [0, 0, 1]>
  >>> mat33 @ gemmi.Vec3(1, 2, 3)
  <gemmi.Vec3(1, -13, 3)>
  >>> mat33.determinant()
  1.0
  >>> mat33.inverse()
  <gemmi.Mat33 [1, 0, 0]
               [0, 1, 5]
               [0, 0, 1]>

(and a few others that are not documented yet).

----

Symmetric matrix SMat33 is implemented as a C++ template
that can work with either 32- or 64-bit floating point numbers.
In Python we have two corresponding classes: SMat33f (32-bit)
and SMat33d (64-bit).
These classes are used primarily for anisotropic ADP tensors;
their member variables are named ``u11``, ``u22``, ``u33``,
``u12``, ``u13`` and ``u23``.

.. doctest::

  >>> aniso = gemmi.read_small_structure('../tests/4003024.cif').sites[2].aniso
  >>> aniso.u11
  0.103
  >>> aniso.elements_pdb()    # (u11, u22, u33, u12, u13, u23)
  [0.103, 0.156, 0.156, 0.0, 0.0, 0.0]
  >>> aniso.elements_voigt()  # (u11, u22, u33, u23, u13, u12)
  [0.103, 0.156, 0.156, 0.0, 0.0, 0.0]

SMat33 provides about a dozen of methods,
including calculations of eigenvalues and eigenvectors.
(This documentation is not complete yet).

.. doctest::

  >>> aniso.trace()
  0.41500000000000004
  >>> aniso.determinant()
  0.002506608
  >>> aniso.calculate_eigenvalues()
  [0.103, 0.156, 0.156]


.. _transform:

Transformations
===============

Working with macromolecular coordinates involves 3D transformations,
such as crystallographic and non-crystallographic symmetry operations,
and fractionalization and orthogonalization of coordinates.

3D transformations tend to be represented either by a 4x4 matrix,
or by a dual quaternion, or by a 3x3 matrix and a translation vector.
Gemmi uses the latter. Transformations are represented by
the ``Transform`` class that has two member variables:
``mat`` (of type ``Mat33``) and ``vec`` (of type ``Vec3``).

.. doctest::

  >>> tr = gemmi.Transform()  # identity
  >>> tr.mat
  <gemmi.Mat33 [1, 0, 0]
               [0, 1, 0]
               [0, 0, 1]>
  >>> tr.vec
  <gemmi.Vec3(0, 0, 0)>


Here is an example that shows a transformation read from a PDB file:

.. doctest::
  :hide:

  >>> import math

.. doctest::

  >>> # get NCS transformation from an example pdb file
  >>> ncs_op = gemmi.read_structure('../tests/1lzh.pdb.gz').ncs[0].tr
  >>> type(ncs_op)
  <class 'gemmi.Transform'>
  >>> ncs_op.mat
  <gemmi.Mat33 [0.97571, -0.2076, 0.06998]
               [0.2156, 0.96659, -0.13867]
               [-0.03885, 0.15039, 0.98786]>
  >>> _.determinant()
  1.0000038877996669
  >>> ncs_op.mat.trace()
  2.93016
  >>> math.degrees(math.acos((_ - 1) / 2))  # calculate rotation angle
  15.186116047571074
  >>> ncs_op.vec
  <gemmi.Vec3(-14.1959, 0.72997, -30.5229)>

  >>> # is the 3x3 matrix above orthogonal?
  >>> mat = ncs_op.mat
  >>> identity = gemmi.Mat33()
  >>> mat.multiply(mat.transpose()).approx(identity, epsilon=1e-5)
  True

  >>> ncs_op.apply(gemmi.Vec3(20, 30, 40))
  <gemmi.Vec3(1.8895, 28.4929, 12.7262)>
  >>> ncs_op.inverse().apply(_)
  <gemmi.Vec3(20, 30, 40)>

To avoid mixing of orthogonal and fractional coordinates
Gemmi also has ``FTransform``, which is like ``Transform``,
but can be applied only to ``Fractional`` coordinates.

.. _box:

Box
===

``Box`` is a small utility for calculation of bounding boxes.
It comes in two variants: for ``Position`` and ``Fractional``:

.. doctest::

  >>> box = gemmi.PositionBox()
  >>> box.extend(gemmi.Position(-5, 5, 0))
  >>> box.extend(gemmi.Position(4, 4, -1))
  >>> box.minimum
  <gemmi.Position(-5, 4, -1)>
  >>> box.maximum
  <gemmi.Position(4, 5, 0)>
  >>> box.get_size()
  <gemmi.Position(9, 1, 1)>
  >>> box.add_margin(0.5)  # changes both minimum and maximum
  >>> box.get_size()
  <gemmi.Position(10, 2, 2)>

  >>> # Fractional variant works in the same way
  >>> box = gemmi.FractionalBox()

In C++ it is a template ``Box<T>`` defined in ``gemmi/math.hpp``.


.. _unitcell:

Unit cell
=========

When working with a structural model in a crystal we need to know
the unit cell. In particular, we use the unit cell to switch between
orthogonal (Cartesian) and fractional coordinates.

The UnitCell class stores the cell parameters
(``a``, ``b``, ``c``, ``alpha``, ``beta``, ``gamma``)
and other properties of the cell precalculated for efficiency
(orthogonalization and fractionalization transformations,
the volume, parameters of the reciprocal unit cell).

Here are the most important properties of UnitCell and a couple of methods
for switching between fractional and Cartesian coordinates:

**C++**

.. literalinclude:: code/cell.cpp

**Python**

.. doctest::

    >>> cell = gemmi.UnitCell(25.12, 39.50, 45.07, 90, 90, 90)
    >>> cell
    <gemmi.UnitCell(25.12, 39.5, 45.07, 90, 90, 90)>
    >>> cell.a, cell.b, cell.c
    (25.12, 39.5, 45.07)
    >>> cell.alpha, cell.beta, cell.gamma
    (90.0, 90.0, 90.0)
    >>> cell.volume
    44720.2568
    >>> cell.frac.mat  # fractionalization matrix
    <gemmi.Mat33 [0.0398089, 0, 0]
                 [0, 0.0253165, 0]
                 [0, 0, 0.0221877]>
    >>> cell.fractionalize(gemmi.Position(10, 10, 10))
    <gemmi.Fractional(0.398089, 0.253165, 0.221877)>
    >>> cell.orth.mat  # orthogonalization matrix
    <gemmi.Mat33 [25.12, 0, 0]
                 [0, 39.5, 0]
                 [0, 0, 45.07]>
    >>> cell.orthogonalize(gemmi.Fractional(0.5, 0.5, 0.5))
    <gemmi.Position(12.56, 19.75, 22.535)>
    >>> cell.orthogonalize_box(box)
    <gemmi.PositionBox object at 0x...>

A symmetry operation that works on fractional coordinates can also be
"orthogonalized" -- converted to :ref:`transformation <transform>`
that operates on Cartesian coordinates:

.. doctest::

    >>> cell.op_as_transform(gemmi.Op('-z,y+1/2,-x'))  #doctest: +ELLIPSIS
    <gemmi.Transform object at 0x...>
    >>> _.apply(gemmi.Position(0, 6, 2.1))
    <gemmi.Vec3(-1.17045, 25.75, 0)>

Cells can be compared with:

* ``approx()`` --- meant for almost identical cells that differ only
  due to numeric errors. It checks if the cell parameters differ
  by less than a given absolute tolerance ε:

  .. doctest::

      >>> cell2 = gemmi.UnitCell(25, 39, 45, 89, 90, 88)
      >>> cell.approx(cell2, epsilon=1e-6)
      False

* ``is_similar()`` --- uses relative tolerance to compare the edge lengths
  and absolute tolerance in degrees to compare the angles:

  .. doctest::

    >>> cell.is_similar(cell2, rel=0.03, deg=2.5)
    True

Next, we can obtain the reciprocal cell:

.. doctest::

    >>> cell.reciprocal()
    <gemmi.UnitCell(0.0398089, 0.0253165, 0.0221877, 90, 90, 90)>

and `metric tensors <https://dictionary.iucr.org/Metric_tensor>`_
in the direct and reciprocal space:

.. doctest::

    >>> cell.metric_tensor()
    <gemmi.SMat33d(631.014, 1560.25, 2031.3, 0, 0, 0)>
    >>> cell.reciprocal_metric_tensor()
    <gemmi.SMat33d(0.00158475, 0.000640923, 0.000492294, 0, 0, 0)>

If the lattice is centered, we can obtain a primitive cell.
We have a function that takes centring type (return value of
``SpaceGroup.centring_type()``), uses matrix from ``centred_to_primitive()``
and returns orthogonalization matrix of a primitive cell:
of the primitive cell:

.. doctest::

    >>> cell.primitive_orth_matrix('I')
    <gemmi.Mat33 [-12.56, 12.56, 12.56]
                 [19.75, -19.75, 19.75]
                 [22.535, 22.535, -22.535]>


This matrix can be used to obtain the G\ :sup:`6` and S\ :sup:`6` vectors,
which are used in Niggli and Selling-Delaunay :ref:`cell reduction <niggli>`.

Function ``is_compatible_with_spacegroup`` checks if the space group
operations don't change the metric tensor elements by more than *ε*
(*ε*\ =0.001 by default):

.. doctest::

    >>> cell.is_compatible_with_spacegroup(gemmi.SpaceGroup('I 2 2 2'))
    True
    >>> cell.is_compatible_with_spacegroup(gemmi.SpaceGroup('P 3'), eps=0.01)
    False


The UnitCell object stores internally (in `UnitCell.images``) a list of
symmetry transformations -- crystallographic symmetry and, in case of
macromolecules, also NCS -- that transform asymmetric unit (ASU) into
the complete unit cell. This list is populated by the class that contains
the UnitCell. It is done automatically when reading a coordinate file.
If you set the unit cell, space group or NCS manually,
call Structure.setup_cell_images() or SmallStructure.setup_cell_images()
to update ``images``.
(The NCS operarations in this list are only those marked as not "given"
in the MTRIX record in the PDB format or in _struct_ncs_oper in mmCIF).

UnitCell.images are used for searching neighbors,
calculating structure factors, and a few other things.
The following functions also rely on it:

* ``UnitCell::volume_per_image() -> double`` -- returns ``UnitCell::volume``
  divided by the number of the molecule images in the unit cell,

  .. doctest::

    >>> st = gemmi.read_structure('../tests/1pfe.cif.gz')
    >>> st.spacegroup_hm
    'P 63 2 2'
    >>> st.cell.volume / st.cell.volume_per_image()
    12.0

* ``UnitCell::is_special_position(const Position& pos, double max_dist=0.8) -> int`` --
  returns the number of nearby symmetry mates of an atom.
  Non-zero only for atoms on special positions.
  For example, returns 3 for an atom on 4-fold symmetry axis.

  .. doctest::

    >>> # chloride ion in 1PFE is significantly off the special position
    >>> cl = st[0].sole_residue('A', gemmi.SeqId('20'))[0]
    >>> cl
    <gemmi.Atom CL at (-0.3, 23.0, -19.6)>
    >>> round(1.0 / cl.occ)
    6
    >>> st.cell.is_special_position(cl.pos, max_dist=0.5)
    0
    >>> st.cell.is_special_position(cl.pos, max_dist=0.8)
    3
    >>> st.cell.is_special_position(cl.pos, max_dist=1.2)
    5

* ``UnitCell::find_nearest_image(const Position& ref, const Position& pos, Asu asu) -> NearestImage`` --
  with the last argument set to ``Asu::Any``,
  it returns the symmetric image of ``pos`` that is nearest to ``ref``.
  The last argument can also be set to ``Asu::Same`` or ``Asu::Different``.

* ``UnitCell::find_nearest_pbc_image(const Position& ref, const Position& pos, int image_idx)`` --
  similar to the function above, but takes the index of symmetry transformation
  as an argument and finds only the unit cell shift. The section about
  :ref:`neighbor search <neighbor_search>` has an example of usage.

The unit cell can be used to determine interplanar spacing *d*:sub:`hkl`
in the reciprocal space (the resolution corresponding to a reflection):

.. doctest::

    >>> cell.calculate_d([0, 1, 0])
    39.5

Computationally, *d* is calculated from 1/*d*:sup:`2`, so if you
need the latter you can calculate it directly:

.. doctest::

    >>> cell.calculate_1_d2([8, -9, 10])
    0.20256818878283983

When changing a symmetry setting of coordinates or reindexing reflections
we need a new unit cell, which can be obtained with one of functions
``changed_basis_forward()`` and ``changed_basis_backward()``:

.. doctest::

    >>> cell.changed_basis_backward(gemmi.Op('y,z,x'), set_images=True)
    <gemmi.UnitCell(45.07, 25.12, 39.5, 90, 90, 90)>

With ``set_images=False`` the ``images`` list in the new unit cell is empty.
With ``True`` -- it contains transformed original list
(but it doesn't work correctly when the cell volume changes).


.. _niggli:

Unit cell reduction
===================

Note: it is about a specialized functionality that few people will ever need.

The reduction finds special bases of lattices. In practice, these bases are
found by removing lattice centering (i.e. obtaining a primitive cell)
and using a prescribed iterative procedure.
As it is worded in the International Tables for Crystallography A 3.1.1.4 (2016),
"the reduction procedures employ metrical properties to develop a sequence
of basis transformations which lead to a *reduced basis* and *reduced cell*".

There are three popular unit cell reductions:

- the Minkowski-Buerger reduction, which minimizes *a*\ +\ *b*\ +\ *c*
  (in special cases multiple, up to 6 different bases have the same
  minimal sum *a*\ +\ *b*\ +\ *c*),
- the Eisenstein-Niggli reduction, which adds extra conditions
  to the previous one and makes the result unique,
- the Selling-Delaunay reduction (the second name is alternatively
  transliterated as Delone), which minimizes
  *a*:sup:`2`\ +\ *b*:sup:`2`\ +\ *c*:sup:`2`\ +\ (*a*\ +\ *b*\ +\ *c*)\ :sup:`2`.

First names here (Minkowski, Eisenstein, Selling) belong to mathematicians
working on the reduction of quadratic forms.
The second names -- to people applying this math to crystallography.
Usually, we use only the second name. The Niggli reduction is the most
popular of the three.

Niggli and Buerger reductions
-----------------------------

Gemmi implements separately the Niggli and Buerger reductions.
The procedures are iterative. Most of the unit cells from the PDB
need only 1-2 iterations to get reduced (1.3 on average, not counting
the *normalization* steps as separate iterations).
On the other hand, one can always construct a primitive cell with extremely
long basis vectors that would require hundreds of iterations.
The Buerger reduction is simpler and faster than Niggli,
but Niggli is also fast -- one iteration takes less than 1μs.

Gemmi implementation is based on the algorithms published by B. Gruber
in the 1970's: Gruber,
`Acta Cryst. A29, 433 <https://doi.org/10.1107/S0567739473001063>`_ (1973)
for the Buerger reduction, and Křivý & Gruber,
`Acta Cryst. A32, 297 <https://doi.org/10.1107/S0567739476000636>`_ (1976)
for the Niggli reduction.
Additionally, the Niggli reduction is using ε to compare numbers, as proposed
by Grosse-Kunstleve *et al*,
`Acta Cryst. A60, 1 <https://doi.org/10.1107/S010876730302186X>`_ (2004).

Gruber's algorithms use vector named G\ :sup:`6`, which is
`somewhat similar <https://dictionary.iucr.org/Metric_tensor>`_
to the metric tensor. G\ :sup:`6` has six elements named:
A, B, C, ξ (xi), η (eta) and ζ (zeta), which correspond to:

    (**a**:sup:`2`, **b**:sup:`2`, **c**:sup:`2`, 2\ **b**\ ⋅\ **c**, 2\ **a**\ ⋅\ **c**, 2\ **a**\ ⋅\ **b**)

Gemmi has a class named GruberVector that contains these six numbers
and reduction algorithms implemented as methods.
This class can be initialized with UnitCell and SpaceGroup:

.. doctest::

  >>> cell = gemmi.UnitCell(63.78, 63.86, 124.40, 90.0, 90.0, 90.0)
  >>> sg = gemmi.SpaceGroup('I 2 2 2')
  >>> gv = gemmi.GruberVector(cell, sg)

or with 6-tuple corresponding to G\ :sup:`6` of a primitive cell:

.. doctest::

  >>> g6_param = gv.parameters  # obtain such a tuple
  >>> gemmi.GruberVector(g6_param)
  <gemmi.GruberVector((5905.34, 5905.34, 5905.34, -7742.79, -7732.57, 3664.69))>

We can check if G\ :sup:`6` already corresponds to a Buerger and Niggli cell:

.. doctest::

  >>> gv.is_niggli()
  False
  >>> gv.is_buerger()
  False

We can access the G\ :sup:`6` parameters as a tuple:

.. doctest::

  >>> gv.parameters
  (5905.337, 5905.337, 5905.337, -7742.7856, -7732.5744, 3664.686)

and obtain the corresponding cell parameters (with angles in degrees):

.. doctest::
  :skipif: sys.platform == 'win32'  # the last digit differs with MSVC

  >>> gv.cell_parameters()  # primitive cell
  (76.84619053668177, 76.84619053668177, 76.84619053668177, 130.96328311485175, 130.89771578326727, 71.92353702711762)

And most importantly, we can reduce the cell.
``niggli_reduce()`` performs the Niggli reduction on G\ :sup:`6`,
returning the number of iterations it took:

.. doctest::

  >>> gv.niggli_reduce()
  3

Now G\ :sup:`6` contains smaller numbers:

.. doctest::

  >>> gv
  <gemmi.GruberVector((4067.89, 4078.10, 5905.34, -4078.10, -4067.89, -0.00))>

To create a new UnitCell with reduced parameters do:

.. doctest::

  >>> gemmi.UnitCell(* gv.cell_parameters())
  <gemmi.UnitCell(63.78, 63.86, 76.8462, 114.551, 114.518, 90)>

or use a helper method:

.. doctest::

  >>> gv.get_cell()
  <gemmi.UnitCell(63.78, 63.86, 76.8462, 114.551, 114.518, 90)>

Similarly, we can perform the Buerger reduction:

.. doctest::

  >>> gv = gemmi.GruberVector(g6_param)
  >>> gv.buerger_reduce()
  3

In this case both functions gave the same result.

.. doctest::

  >>> gv.is_niggli()
  True
  >>> gv.get_cell()
  <gemmi.UnitCell(63.78, 63.86, 76.8462, 114.551, 114.518, 90)>

Functions ``niggli_reduce``, ``is_niggli`` and ``is_buerger`` can take optional
parameter ``epsilon`` (default: 1e-9) that is used for comparing numbers.
Additionally, ``niggli_reduce`` can take ``iteration_limit`` (default: 100).
To check how the computations would work without ε we can set it to 0:

.. doctest::

  >>> gv.is_buerger(epsilon=0)
  True
  >>> gv.is_niggli(epsilon=0)
  False
  >>> gv.niggli_reduce(epsilon=0, iteration_limit=100)
  6
  >>> gv.get_cell()
  <gemmi.UnitCell(63.78, 63.86, 76.8462, 114.551, 114.518, 90)>

Here, the Niggli conditions were initially found not fulfilled, because
one expression that should be non-negative was about -5e-13.
A few extra iterations sorted it out (without any real changes),
but it's not always the case -- that's why we have ``iteration_limit``
to prevent infinite loop.

The original Křivý-Gruber algorithm doesn't calculate the change-of-basis
transformation that leads to the reduced cell. In gemmi,
this transformation can be obtained as proposed in the 2004 paper
of Grosse-Kunstleve *et al*: the change-of-basis matrix is updated
in each step together with the Gruber vector.
Updating this matrix makes the reduction twice slower
(but it's still in tens of ns, so it's fast enough for any purpose).
To track the change of basis, pass the following option:

.. doctest::

  >>> gv = gemmi.GruberVector(cell, sg, track_change_of_basis=True)

After the Niggli reduction, the transformation will be available
in the ``change_of_basis`` property:

.. doctest::

  >>> gv.niggli_reduce()
  3
  >>> cob = gv.change_of_basis
  >>> cob
  <gemmi.Op("x-z/2,y-z/2,z/2")>

This operator transforms Niggli cell to the original cell
(so it's actually *the inverse* of the reduction change-of-basis):

.. doctest::

  >>> gv.get_cell().changed_basis_forward(cob, set_images=False)
  <gemmi.UnitCell(63.78, 63.86, 124.4, 90, 90, 90)>

and the other way around:

.. doctest::

  >>> cell.changed_basis_backward(cob, set_images=False)
  <gemmi.UnitCell(63.78, 63.86, 76.8462, 114.551, 114.518, 90)>

Currently, tracking is implemented only for the Niggli reduction,
not for the Buerger reduction.

Selling-Delaunay reduction
--------------------------

Gemmi implementation is based on

- section `3.1.2.3 <https://onlinelibrary.wiley.com/iucr/itc/Ac/ch3o1v0001/>`_
  "Delaunay reduction and standardization" in the Tables vol. A (2016),
- Patterson & Love (1957), "Remarks on the Delaunay reduction",
  `Acta Cryst. 10, 111 <https://doi.org/10.1107/S0365110X57000328>`_,
- Andrews *et al* (2019),
  "Selling reduction versus Niggli reduction for crystallographic lattices",
  `Acta Cryst. A75, 115 <https://doi.org/10.1107/S2053273318015413>`_.

Similarly to the GruberVector, here we have a class named SellingVector
that contains the six elements of S\ :sup:`6` -- the inner products
among the four vectors **a**, **b**, **c**
and **d**\ =–(\ **a**\ +\ **b**\ +\ **c**\ ):

    *s*\ :sub:`23`\ =\ **b**\ ⋅\ **c**,
    *s*\ :sub:`13`\ =\ **a**\ ⋅\ **c**,
    *s*\ :sub:`12`\ =\ **a**\ ⋅\ **b**,
    *s*\ :sub:`14`\ =\ **a**\ ⋅\ **d**,
    *s*\ :sub:`24`\ =\ **b**\ ⋅\ **d**,
    *s*\ :sub:`34`\ =\ **c**\ ⋅\ **d**.

SellingVector can be initialized with UnitCell and SpaceGroup:

.. doctest::

  >>> sv = gemmi.SellingVector(cell, sg)

or with a tuple of six numbers S\ :sup:`6`:

.. doctest::

  >>> sv.parameters
  (-3871.3928, -3866.2872, 1832.343, -3871.3928, -3866.2872, 1832.343)
  >>> gemmi.SellingVector(_)
  <gemmi.SellingVector((-3871.39, -3866.29, 1832.34, -3871.39, -3866.29, 1832.34))>

Similarly as in the previous section, we can check if S\ :sup:`6`
already corresponds to a Delaunay cell:

.. doctest::

  >>> sv.is_reduced()
  False

Each reduction step decreases Σ\ **b**\ :sub:`i`:sup:`2`
(**b**\ :sub:`1`, **b**\ :sub:`2`, **b**\ :sub:`3` and **b**\ :sub:`4`
are alternative symbols for **a**, **b**, **c** and **d**).
The sum Σ\ **b**\ :sub:`i`:sup:`2` can be calculated with:

.. doctest::

  >>> sv.sum_b_squared()
  23621.348

Similarly to ``niggli_reduce()``, the Selling reduction procedure takes
optional arguments ``epsilon`` and ``iteration_limit``
and returns the iteration count:

.. doctest::

  >>> sv.reduce()
  2

Now we can check the result:

.. doctest::

  >>> sv
  <gemmi.SellingVector((-2033.94, -2033.94, -1832.34, -2039.05, -2039.05, 0.00))>
  >>> sv.is_reduced()
  True
  >>> sv.sum_b_squared()
  19956.662

Now, the corresponding four vectors can be in any order.
We may sort them so that *a*\ ≤\ *b*\ ≤\ *c*\ ≤\ *d*:

.. doctest::

  >>> sv.sort()
  >>> sv
  <gemmi.SellingVector((-2039.05, -2033.94, 0.00, -2033.94, -2039.05, -1832.34))>

Finally, we can get the corresponding UnitCell:

.. doctest::

  >>> gemmi.UnitCell(* sv.cell_parameters())
  <gemmi.UnitCell(63.78, 63.86, 76.8462, 114.551, 114.518, 90)>
  >>> sv.get_cell()  # helper function that does the same
  <gemmi.UnitCell(63.78, 63.86, 76.8462, 114.551, 114.518, 90)>

S\ :sup:`6` can be used to calculate G\ :sup:`6`, and the other way around:

.. doctest::

  >>> sv.gruber()
  <gemmi.GruberVector((4067.89, 4078.10, 5905.34, -4078.10, -4067.89, 0.00))>
  >>> _.selling()
  <gemmi.SellingVector((-2039.05, -2033.94, 0.00, -2033.94, -2039.05, -1832.34))>

TBC

