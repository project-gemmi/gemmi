Symmetry
########

This section describes functionality related to the 3D space groups.

In C++: ``#include <gemmi/symmetry.hpp>`` (it's header-only).

Space group table
=================

Gemmi tabulates 550+ settings of the 230 crystallographic space groups.
Each entry includes:

* ``number`` -- space group number (1-230),
* ``ccp4`` -- ccp4 number (assigned to particular settings; modulo 1000
  they give space group number: 3 and 1003 correspond to
  ``P 1 2 1`` and ``P 1 1 2``; 0 means none),
* ``hm`` -- Hermann-Mauguin (H-M) symbol a.k.a. the international notation
  (``I a -3 d``, ``C 1 2 1``),
* ``ext`` -- extension to the H-M notations (none, ``1``, ``2``, ``H``
  or ``R``) that make extended H-M symbol (``R 3:H``, ``P 4/n:1``),
* ``hall`` -- the Hall symbol (``-I 4bd 2c 3``, ``C 2y``, ``P 32 2c (0 0 -1)``)
  used to generate the space group operations,
* and ``basisop`` -- change of basis operator from the reference setting.

This data is derived primarily from the CCP4 :file:`syminfo.lib` file,
which in turn is based on the data from sgtbx_ that was augmented
with the old data from a CCP4 file named :file:`symop.lib`.

The data from sgtbx is also available in the older SgInfo_ library,
as well as in the `International Tables <http://it.iucr.org/>`_
for Crystallography Vol. B ch. 1.4 (in 2010 ed.). It has 530 entries
including 2 duplicates (3 different names for the same settings)
in the space group 68.

Gemmi includes also settings from OpenBabel_ that are absent in
:file:`syminfo.lib`. It does not include all possible settings,
but if needed, more entries can be added.
For example, we do not include all
the C- and F-centred tetragonal space groups featured in
`Crystallographic Space Group Diagrams and Tables <http://img.chem.ucl.ac.uk/sgp/mainmenu.htm>`_
by Jeremy Karl Cockcroft
(they are also listed on
`this page <https://cci.lbl.gov/cctbx/multiple_cell.html>`_
written by R.W. Grosse-Kunstleve,
and mentioned in the 2015 edition of ITfC Vol.A, Table 1.5.4.4).

We also tabulated alternative names.
For now this only includes new standard names introduced in 1990's by the IUCr
committee. For example, the space group no. 39 is now officially, in the
Volume A of the `International Tables <http://it.iucr.org/>`_,
named Aem2 not Abm2.
Most of the crystallographic software (as well as ITfC Vol.B)
still uses the old names.
(spglib_ uses new ones, sgtbx reads new names with the option ``ad_hoc_1992``).

The usual way to access a space group from the table is to search
it by name. In C++::

  #include <gemmi/symmetry.hpp>
  // ...
  const SpaceGroup* sg = find_spacegroup_by_name(name);

and in Python:

.. doctest::

  >>> import gemmi
  >>> gemmi.find_spacegroup_by_name('I2')
  <gemmi.SpaceGroup("I 1 2 1")>
  >>> # or just
  >>> gemmi.SpaceGroup('I2')
  <gemmi.SpaceGroup("I 1 2 1")>

.. note::

    The rest of this section has only Python examples mixed with the text.
    One longer C++ example is at the end.

The name above is expected to be either a full international Hermann-Mauguin
symbol or a short symbol (*I2* instead of *I 1 2 1*).
It can be an "alternative" name:

.. doctest::

    >>> gemmi.find_spacegroup_by_name('C m m e') # new names have 'e' and 'g'
    <gemmi.SpaceGroup("C m m a")>

Sometimes in the PDB, the setting of the hexagonal crystal family
is not clear from the symmetry symbol alone. For instance, the H-M symbol
"R 3" can mean either hexagonal or rhombohedral setting.
Furthermore, ITfC vol. A, PDB and SCALEPACK assign different meanings to
"H 3" and similar symbols, as described in
`Fuzzy space group symbols: H3 and H32 <http://www.phenix-online.org/newsletter/CCN_2011_01.pdf>`_.
This ambiguity can be resolved by comparing angles of the unit cell.
The ratio of *γ* to *α* angles is 120:90 in the hexagonal system
and 1:1 in rhombohedral. (Note: gemmi assumes that H3 refers to
space group 146 (R3), not 143 (P3)). Therefore, ``find_spacegroup_by_name()``
accepts also *α* and *γ* angles. If the angles are not passed,
the hexagonal system is returned:

.. doctest::

    >>> gemmi.find_spacegroup_by_name('R 3 2')
    <gemmi.SpaceGroup("R 3 2:H")>
    >>> gemmi.find_spacegroup_by_name('R 3 2', alpha=92.02, gamma=92.02)
    <gemmi.SpaceGroup("R 3 2:R")>
    >>> gemmi.find_spacegroup_by_name('R 3 2', alpha=90, gamma=120)
    <gemmi.SpaceGroup("R 3 2:H")>
    >>> # of course, you do not need angles if you use extended H-M symbol
    >>> gemmi.find_spacegroup_by_name('R 3 2:R')
    <gemmi.SpaceGroup("R 3 2:R")>

You can also get space group by number:

.. doctest::

  >>> gemmi.find_spacegroup_by_number(5)
  <gemmi.SpaceGroup("C 1 2 1")>

The number is the ccp4 number mentioned above,
so some 4-digit numbers are also recognized:

.. doctest::

  >>> gemmi.find_spacegroup_by_number(4005)
  <gemmi.SpaceGroup("I 1 2 1")>

For values 1-230 the number corresponds to the first setting of
this space group in the table in ITfC vol. B.
Usually, it is the reference (standard) setting of the space group,
but unfortunately not always. So we have another function that returns
always the reference setting:

.. doctest::

  >>> gemmi.get_spacegroup_reference_setting(48)
  <gemmi.SpaceGroup("P n n n:2")>

The last function for searching the space group table checks
symmetry operations:

.. doctest::

  >>> gemmi.symops_from_hall('C 2y (x,y,-x+z)')  #doctest: +ELLIPSIS
  <gemmi.GroupOps object at 0x...>
  >>> gemmi.find_spacegroup_by_ops(_)
  <gemmi.SpaceGroup("I 1 2 1")>
  >>> _.hall
  'I 2y'

This example shows also how to find space group corresponding to a Hall symbol.
The Hall notation encodes all the group operations.
Unfortunately, in a non-unique way.
Different Hall symbols can be used to encode the same symmetry operations.
In the example above "C 2y (x,y,-x+z)" is equivalent to "I 2y".
That's why we compare operations not symbols.

Finally, we may iterate over the space group table:

.. doctest::

  >>> for sg in gemmi.spacegroup_table():
  ...     if sg.ext == 'H':
  ...         print(sg)
  <gemmi.SpaceGroup("R 3:H")>
  <gemmi.SpaceGroup("R -3:H")>
  <gemmi.SpaceGroup("R 3 2:H")>
  <gemmi.SpaceGroup("R 3 m:H")>
  <gemmi.SpaceGroup("R 3 c:H")>
  <gemmi.SpaceGroup("R -3 m:H")>
  <gemmi.SpaceGroup("R -3 c:H")>

If you would like to ignore entries that are absent in SgInfo, sgtbx, spglib
and in the International Tables vol. B, use only the first 530 entries
of the table. In Python, we have a helper function for this:

.. doctest::

  >>> for sg in gemmi.spacegroup_table_itb():
  ...     pass


SpaceGroup
==========

``gemmi.SpaceGroup`` represents an entry in the space group table.
It has the properties listed at the beginning of this section:

.. doctest::

  >>> sg = gemmi.SpaceGroup('R3')
  >>> sg.number
  146
  >>> sg.ccp4
  146
  >>> sg.hm
  'R 3'
  >>> sg.ext
  'H'
  >>> sg.hall
  'R 3'
  >>> sg.basisop
  <gemmi.Op("x,y,z")>

and a few methods:

.. doctest::

  >>> sg.xhm()                    # extended Hermann-Mauguin name
  'R 3:H'
  >>> sg.short_name()             # short name
  'H3'
  >>> sg.is_enantiomorphic()      # is it one of 22 chiral space groups?
  False
  >>> sg.is_sohncke()             # is it one of 65 Sohncke space groups?
  True
  >>> sg.is_symmorphic()          # is it one of 73 symmorphic space groups?
  True
  >>> sg.is_centrosymmetric()     # does it have inversion?
  False
  >>> sg.point_group_hm()         # H-M name of the point group
  '3'
  >>> sg.laue_str()               # name of the Laue class
  '-3'
  >>> sg.crystal_system_str()     # name of the crystal system
  'trigonal'
  >>> sg.is_reference_setting()
  True
  >>> sg.centring_type()          # lattice centering type
  'R'
  >>> sg.centred_to_primitive()   # change-of-basis operator to a primitive lattice
  <gemmi.Op("2/3*x-y/3-z/3,x/3+y/3-2/3*z,x/3+y/3+z/3")>
  >>> # and, most importantly, the symmetry operations
  >>> sg.operations()             #doctest: +ELLIPSIS
  <gemmi.GroupOps object at 0x...>

Chirality-related functions ``is_enantiomorphic()`` and ``is_sohncke()``
can be confusing. Here, we follow the IUCr dictionary:

* `Sohncke groups <https://dictionary.iucr.org/Sohncke_groups>`_
  (a.k.a. non-enantiogenic space groups)
  are "the three-dimensional space groups containing only operations
  of the first kind (rotations, rototranslations, translations)";
  65 groups in which chiral structures crystallize.
* Enantiomorphic space groups
  (a.k.a. `chiral groups <https://dictionary.iucr.org/Chiral_space_group>`_)
  are those whose *group* structure is chiral;
  22 groups forming 11 enantiomorphic pairs.
  (So chiral structures can crystallize not only in the chiral space groups
  but also in 43 of the achiral ones.)


Operations
==========

Crystallographic symmetry operations have a few notations.
Gemmi undertands only coordinate triplets (sometimes called
the Jones' faithful notation) such as ``x,x-y,z+1/2``.
The symmetry operation is represented by class Op.
It can be created in C++ as::

  #include <gemmi/symmetry.hpp>
  ...
  gemmi::Op op = gemmi::parse_triplet("-y,x-y,z+1/3");

and in Python as:

.. doctest::

    >>> import gemmi
    >>> op = gemmi.Op('-y,x-y,z+1/3')

We can also do the opposite:

.. doctest::

    >>> op.triplet()
    '-y,x-y,z+1/3'

Alternatively, the letters can be hkl or abc, or uppercase:

.. doctest::

    >>> op.triplet(style='h')
    '-k,h-k,l+1/3'
    >>> op.triplet(style='a')
    '-b,a-b,c+1/3'
    >>> op.triplet(style='X')
    '-Y,X-Y,Z+1/3'

The operation consists of a 3x3 rotation matrix and
a translation vector, both stored internally as integers that need to be
divided by ``DEN`` == 24 to get the actual values.

.. doctest::

    >>> op.rot
    [[0, -24, 0], [24, -24, 0], [0, 0, 24]]
    >>> op.tran
    [0, 0, 8]

Alternatively, the operation can be expressed as a single 4x4 transformation
matrix (which in crystallography is called *Seitz matrix*).

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> op.seitz()
    [[0, -1, 0, 0], [1, -1, 0, 0], [0, 0, 1, Fraction(1, 3)], [0, 0, 0, 1]]
    >>> op.float_seitz()
    [[0.0, -1.0, 0.0, 0.0],
     [1.0, -1.0, 0.0, 0.0],
     [0.0, 0.0, 1.0, 0.3333333333333333],
     [0.0, 0.0, 0.0, 1.0]]

The type of the rotation part of operation can be determined with rot_type(),
which it based on Table 1 from the
`RWGK's space-group algorithms paper <https://doi.org/10.1107/S0108767398010186>`_.
It returns integer N, meaning N-fold rotation for positive N
and rotoinversion for negative N:

.. doctest::

    >>> gemmi.Op('-y,x-y,z').rot_type()
    3
    >>> gemmi.Op('x,y,-z').rot_type()
    -2

Operations can be combined, inverted and wrapped:

.. doctest::

    >>> gemmi.Op('x-y,x,z+1/6') * '-x,-y,z+1/2'
    <gemmi.Op("-x+y,-x,z+2/3")>
    >>> _.inverse()
    <gemmi.Op("-y,x-y,z-2/3")>
    >>> _.wrap()
    <gemmi.Op("-y,x-y,z+1/3")>

Wrapping applies *modulo* 1 to the translational part.
Which is usually desirable for crystallographic symmetry.
But the triplets and matrices may represent also, for example,
generation of a biological assembly.
For this reason, ``x,y+1,z`` is not automatically reduced to the identity.
But when operations are combined,
they are assumed to be symmetry operations and the result is wrapped to [0,1):

.. doctest::

    >>> op
    <gemmi.Op("-y,x-y,z+1/3")>
    >>> op * op
    <gemmi.Op("-x+y,-x,z+2/3")>
    >>> op * op * op  # without wrapping we'd have z+1
    <gemmi.Op("x,y,z")>

The ``Op.rot`` matrix is called "rotation matrix", because that's the primary
purpose, but it can also represent different linear transformations:

.. doctest::

    >>> enlarging_op = gemmi.Op("-y+z,x+z,-x+y+z")
    >>> enlarging_op.inverse()
    <gemmi.Op("-x/3+2/3*y-z/3,-2/3*x+y/3+z/3,x/3+y/3+z/3")>
    >>> _ * enlarging_op
    <gemmi.Op("x,y,z")>

In the real space, a crystal symmetry operation can be applied
to the fractional atom position to get an equivalent position:

.. doctest::

    >>> op.apply_to_xyz([0.25, 0.21875, 0.3])
    [-0.21875, 0.03125, 0.6333333333333333]


In the reciprocal space, the same operation relates equivalent reflections.
The rotational part determines Miller indices,
the translational part -- phase shift.

.. doctest::

    >>> hkl = [3, 0, 1]
    >>> op.apply_to_hkl(hkl)
    [0, -3, 1]
    >>> op.phase_shift(hkl)  # -120 degrees in radians
    -2.0943951023931953


Groups of Operations
====================

Each space group setting corresponds to a unique set of operations.
This set is represented by class ``GroupOps``.

Symmetry operations (rotation + translation) and
centring vectors (translation only) are stored separately:

.. doctest::

  >>> sg = gemmi.SpaceGroup('I2')
  >>> ops = sg.operations()
  >>> ops  #doctest: +ELLIPSIS 
  <gemmi.GroupOps object at 0x...>
  >>> list(ops.sym_ops)
  [<gemmi.Op("x,y,z")>, <gemmi.Op("-x,y,-z")>]
  >>> list(ops.cen_ops)
  [[0, 0, 0], [12, 12, 12]]

but they can be combined on the fly:

.. doctest::

  >>> len(ops)
  4
  >>> for op in ops:
  ...   print(op.triplet())
  ...
  x,y,z
  -x,y,-z
  x+1/2,y+1/2,z+1/2
  -x+1/2,y+1/2,-z+1/2

We can apply a change-of-basis operator to GroupOps:

.. doctest::

  >>> sg.basisop  # I2 -> C2
  <gemmi.Op("x,y,-x+z")>
  >>> ops.change_basis_forward(sg.basisop)
  >>> gemmi.find_spacegroup_by_ops(ops)
  <gemmi.SpaceGroup("C 1 2 1")>
  >>> ops.change_basis_backward(sg.basisop)  # and the other way around
  >>> gemmi.find_spacegroup_by_ops(ops)
  <gemmi.SpaceGroup("I 1 2 1")>

In particular, we can switch between enantiomorphic pairs using inversion:

.. doctest::

  >>> ops = gemmi.SpaceGroup('P 41').operations()
  >>> ops.change_basis_forward(gemmi.Op('-x,-y,-z'))
  >>> gemmi.find_spacegroup_by_ops(ops)
  <gemmi.SpaceGroup("P 43")>

and change to primitive space group using ``centred_to_primitive()``:

.. doctest::

  >>> sg = gemmi.SpaceGroup('R 3:H')
  >>> sg.centring_type()
  'R'
  >>> ops = sg.operations()
  >>> ops.change_basis_backward(sg.centred_to_primitive())
  >>> gemmi.find_spacegroup_by_ops(ops)
  <gemmi.SpaceGroup("R 3:R")>
  >>> _.centring_type()
  'P'

We can create GroupOps from a list of operations
(in C++ use function ``split_centering_vectors``):

.. doctest::

  >>> op_list = ['x,y,z', 'x,-y,z+1/2', 'x+1/2,y+1/2,z', 'x+1/2,-y+1/2,z+1/2']
  >>> new_ops = gemmi.GroupOps([gemmi.Op(o) for o in op_list])

or from a Hall symbol:

.. doctest::

  >>> gemmi.symops_from_hall('P 4w 2c')  #doctest: +ELLIPSIS
  <gemmi.GroupOps object at 0x...>
  >>> len(_)
  8

The Hall symbols encode *generators* which are then used to obtain
all the operations. If you'd wonder what generators are encoded, use:

.. doctest::

    >>> gen = gemmi.generators_from_hall('P 4w 2c')  # no.93
    >>> gen.cen_ops  # centering P
    [[0, 0, 0]]
    >>> gen.sym_ops
    [<gemmi.Op("x,y,z")>, <gemmi.Op("-y,x,z+1/4")>, <gemmi.Op("x,-y,-z+1/2")>]

Combining these 3 generators reconstructs all the 8 symmetry operations:

.. doctest::

    >>> gen.add_missing_elements()
    >>> len(gen)
    8

A GroupOps object can be used to search the space group table for a matching
space group:

.. doctest::

  >>> gemmi.find_spacegroup_by_ops(new_ops)
  <gemmi.SpaceGroup("C 1 c 1")>

To check only the lattice centering we can use:

.. doctest::

  >>> new_ops.find_centering()
  'C'

We can check if the operations contain inversion:

.. doctest::

  >>> new_ops.is_centrosymmetric()
  False

and we can tell which reflections are centric (as opposed to acentric; a reflection
is centric in the given space group if its Friedel mate (-*h*,-*k*,-*l*) is equivalent
to it by symmetry):

.. doctest::

  >>> new_ops.is_reflection_centric([1, 2, 3])
  False
  >>> new_ops.is_reflection_centric([0, 2, 0])
  True

Similarly, we can check for systematic absences:

.. doctest::

  >>> new_ops.is_systematically_absent([1, 2, 3])
  True
  >>> new_ops.is_systematically_absent([1, 3, 2])
  False

We can calculate the epsilon factor ε, which tells how many times the symmetry
operations map the reflection onto itself:

.. doctest::

  >>> new_ops.epsilon_factor([1, 3, 2])
  2
  >>> new_ops.epsilon_factor([2, 0, 2])
  4

We also have a function that calculates ε ignoring centering vectors
(equivalent to the ``epsilon()`` function in cctbx):

.. doctest::

  >>> new_ops.epsilon_factor_without_centering([1, 3, 2])
  1
  >>> new_ops.epsilon_factor_without_centering([2, 0, 2])
  2

We can create a new GroupOps object with the translational components of
glide planes and screw axes removed. Such a new set of operations corresponds
to a `symmorphic <https://dictionary.iucr.org/Symmorphic_space_groups>`_
space group.

.. doctest::

  >>> new_ops.derive_symmorphic()  #doctest: +ELLIPSIS
  <gemmi.GroupOps object at 0x...>
  >>> gemmi.find_spacegroup_by_ops(_)
  <gemmi.SpaceGroup("C 1 m 1")>

And can add inversion (-x,-y,-z). This function either doubles the number
of operations and returns True, or (if the group already has inversion)
it returns False leaving the group unchanged:

.. doctest::

  >>> new_ops.add_inversion()
  True
  >>> gemmi.find_spacegroup_by_ops(new_ops)
  <gemmi.SpaceGroup("C 1 2/c 1")>


ASU
===

The asymmetric unit (ASU) of a space group is a non-redundant part of the unit cell,
a part that can be used to generate the complete unit cell by application
of the symmetry operations.

Direct space
------------

ASU as a geometric shape can be defined by a set of inequalities -- planes
that cut out a convex volume. Usually, crystallographic computations
do not require ASU to be convex or even contiguous, so gemmi implements
a more simplistic approach.

.. _asu_brick:

Similarly to other software packages, Gemmi uses ASU bricks --
parallelepipeds (cuboids in fractional coordinates) that contains an ASU.
A brick may contain more than one ASU if a brick-shaped ASU is not possible.
Both the ASU and the brick can be chosen in many ways.
For example, for P 21 21 21 Gemmi will use the following brick:

.. doctest::

  >>> brick = gemmi.find_asu_brick(gemmi.SpaceGroup('P 21 21 21'))
  >>> brick.str()
  '0<=x<1/2; 0<=y<=1/2; 0<=z<1'

while the CCP4 and CCTBX libraries use, respectively, these definitions:

|  0<=x<1; 0<=y<1; 0<=z<=1/4
|  0<=x<1; 0<=y<=1/4; 0<=z<1

The brick description is stored in an object of class AsuBrick.
The lower boundary in always at 0 (inclusive), the upper one is defined by:

.. doctest::

  >>> brick.size  # multiplied by 24 to get integers
  [12, 12, 24]
  >>> brick.incl  # which upper boundaries are included
  [False, True, False]

This corresponds to the human-readable string:

.. doctest::

  >>> print(brick.str())
  0<=x<1/2; 0<=y<=1/2; 0<=z<1

AsuBrick has also function get_extent() that can be used to
:ref:`set extent of a CCP4 map <set_extent>`:

.. doctest::

  >>> brick.get_extent()  #doctest: +ELLIPSIS
  <gemmi.FractionalBox object at 0x...>
  >>> _.minimum, _.maximum
  (<gemmi.Fractional(-1e-09, -1e-09, -1e-09)>, <gemmi.Fractional(0.5, 0.5, 1)>)

The exact maximum values of the FractionalBox are slightly above or below
the values printed here, depending on whether the boundary is included.

In crystallography, we often work with real-space data
(such as electron density) on a 3D grid.
In Gemmi, to make calculations simpler, such a grid (class Grid)
spans over the whole unit cell.
To process only points in an ASU we use an ASU mask -- this is described
in the section about :ref:`MaskedGrid <masked_grid>`.
Such an ASU, generated by ``masked_asu()``, is inside the ASU brick,
but for some space groups it is not contiguous (which shouldn't be a problem).

.. _reciprocal_asu:

Reciprocal space
----------------

The reciprocal asymmetric unit can also be chosen in various ways.
Gemmi is consistent with both CCP4 and cctbx here (but some programs,
such as TNT, use different ASU choice).

.. doctest::

  >>> p2 = gemmi.SpaceGroup('P 1 2 1')
  >>> asu = gemmi.ReciprocalAsu(p2)
  >>> asu.condition_str()
  'k>=0 and (l>0 or (l=0 and h>=0))'

The condition returned from ``condition_str()`` refers to the standard settings
of the space group, so it is the same for P 1 2 1 and P 1 1 2.
Other functions depend on the settings used.

One can check if a reflection is in the ASU:

.. doctest::

  >>> asu.is_in([1, -2, 3])
  False
  >>> gemmi.ReciprocalAsu(gemmi.SpaceGroup('P 1 1 2')).is_in([1, -2, 3])
  True

and what is the equivalent reflection in the ASU:

.. doctest::

  >>> asu.to_asu([1, -2, 3], p2.operations())
  ([1, 2, 3], 4)

``to_asu()`` returns also index of the symmetry operation between
the original and the returned reflection -- the same number as ISYM in
the MTZ format, odd for reflections in the positive asu (I+),
even for negative (I-).
The second argument (GroupOps) is passed explicitly to avoid determining
space group operations many times when ``to_asu()`` is in a loop.

Twinning
========

Merohedral twinning (more correctly:
`twinning by merohedry <https://dictionary.iucr.org/Twinning_by_merohedry>`_)
has the twin operator belonging to the point group of the lattice
but not to the point group of the crystal. So to find potential twinning
operators we first need to find the symmetry of the lattice.

Determination of the lattice symmetry in gemmi is based on the section 2.1
in `P.H. Zwart et al (2006) <http://legacy.ccp4.ac.uk/newsletters/newsletter44/articles/explore_metric_symmetry.html>`_,
which in turn is based on methods from
`A. Lebedev et al. (2006) <https://doi.org/10.1107/S0907444905036759>`_
and `Y. Le Page (1982) <https://doi.org/10.1107/S0021889882011959>`_.
They:

* determine a :ref:`reduced cell <niggli>`
  (usually a Niggli cell is used, but it can be any Buerger cell),
* find exact and approximate lattice symmetries, with a user-specified
  threshold on the deviation from perfect symmetry,
* use the found symmetries and the change-of-basis operator from the cell
  reduction to obtain the lattice symmetry group in the original unit cell.


In gemmi, lattice symmetry can be obtained with function
find_lattice_symmetry() that takes a unit cell, centering, and obliquity
threshold in degrees.
The threshold limits the obliquity (δ angle as defined by Le Page)
of 2-fold rotations that are used as generators for the symmetry group.
The function returns a GroupOps object with all the lattice symmetry operations
except inversion (you can call GroupOps.add_inversion() to add it).
Rotations are in ``sym_ops`` and the cell centering vectors
(unimportant) are in ``cen_ops``.

.. doctest::

  >>> cell = gemmi.UnitCell(174.22, 53.12, 75.17, 90.0, 115.29, 90.0)
  >>> gemmi.find_lattice_symmetry(cell, centring='C', max_obliq=3)  #doctest: +ELLIPSIS
  <gemmi.GroupOps object at 0x...>
  >>> for op in _.sym_ops: print(op.triplet())
  x,y,z
  -x,y,-z
  -x,-y,-2*x+z
  x,-y,2*x-z

Internally, find_lattice_symmetry() reduces the unit cell and calls
lower-level function find_lattice_symmetry_r().
Here is how we could do it manually:

.. doctest::

  >>> gv = gemmi.GruberVector(cell, 'C', track_change_of_basis=True)
  >>> gv.niggli_reduce()
  2
  >>> reduced_cell = gv.get_cell()
  >>> reduced_ops = gemmi.find_lattice_symmetry_r(reduced_cell, 3)
  >>> reduced_ops.change_basis_forward(gv.change_of_basis)
  >>> for op in reduced_ops.sym_ops: print(op.triplet())
  x,y,z
  -x,y,-z
  -x,-y,-2*x+z
  x,-y,2*x-z
  >>> reduced_ops.cen_ops  # just for the record, centering is here:
  [[0, 0, 0], [12, 12, 0]]

The work of find_lattice_symmetry_r() begins with calling
find_lattice_2fold_ops() that checks 2-fold symmetries enumerated by P. Zwart
and returns the matching ones together with the δ angles:

.. doctest::

  >>> for (op, delta) in gemmi.find_lattice_2fold_ops(reduced_cell, max_obliq=3):
  ...   print(op.triplet(), ' -> obliquity %.2f degrees' % delta)
  x-z,-y,-z  -> obliquity 0.00 degrees
  -x,y-z,-z  -> obliquity 0.27 degrees
  -x+z,-y+z,z  -> obliquity 0.27 degrees

The last two values (δ ≠ 0) are pseudo-symmetries that may result in
twinning by pseudo-merohedry.

After this introduction, here is the function that determines
potential (pseudo-)merohedral twinning operators:

.. doctest::

  >>> gemmi.find_twin_laws(cell, gemmi.SpaceGroup('C2'), max_obliq=3, all_ops=False)
  [<gemmi.Op("-x,-y,-2*x+z")>]

With all_ops=False, this function returns only non-redundant twinning
operators (coset representatives), which is usually what is needed.
Note that since any of the equivalent operators can be returned,
other implementations (class twin_laws in CCTBX, libT.f in CCP4)
may return different ones.

Above, we have one possible twinning. It is by pseudo-merohedry.
If we wanted to see only twinning by merohedry, the threshold should be
near zero, allowing only for numerical imprecision:

.. doctest::

  >>> gemmi.find_twin_laws(cell, gemmi.SpaceGroup('C2'), max_obliq=1e-6, all_ops=False)
  []
  >>> gemmi.find_twin_laws(cell, gemmi.SpaceGroup('P1'), max_obliq=1e-6, all_ops=False)
  [<gemmi.Op("-x,y,-z")>]


C++ Example
===========

Since the code snippets above were only in Python,
here we compensate it with one longer example in C++.

.. literalinclude:: code/sym.cpp
   :language: cpp


Implementation notes
====================

When creating a space group table, one can either store an explicit list
of all the operations for each space group, or generate the operations
from a smaller set of symbols. It is a
trade-off between simplicity of the code and the amount of the tabulated data.

The number of symmetry operations per space group is between 1 and 192,
but they can be split into symmetry operations (max. 48 for point group m-3m)
and so-called *centring vectors* (max. 4 for the face-centered lattice).

* The simplest way of storing the operations is to list them all (i.e. 192
  triplets for no. 228) in a text file.
  This approach is used by OpenBabel_ (in :file:`space-groups.txt`).

* It makes sense to keep the centring vectors separately
  (192 becomes 48 + 4).
  This is done in the CCP4 :file:`syminfo.lib` file (539 entries),
  which is used by csymlib_ (part of libccp4) and a few other projects.

* The operations can be as well tabulated in the code.
  This approach is used (with one or two layers of indirection to reduce
  the data size) by spglib_ (C library) and NGL_.

* The inversion center (if applicable) can be kept separately,
  to reduce the maximum number of stored operations (48 -> 24+1).
  This is how the symmetry data is encoded in Fityk_.

* Actually all the operations can be generated from only a few generators,
  at the expense of more complex code.
  Example: Mantid_ (:file:`SpaceGroupFactory.cpp`).

* Finally, one can use one of the two computer-adapted descriptions from ITfC.
  The so-called explicit notation (``ICC$I3Q000$P4C393$P2D933``) is the
  longer of the two, but easier to parse by the computer.
  It is used in the SPGGEN_ program.

* The Hall notation (``-I 4bd 2c 3``), first proposed by Sydney R. Hall
  in 1981, is shorter and more popular.
  It can be interpreted by a few libraries:

  * SgInfo_ and SgLite_ (old C libraries from Ralf W. Grosse-Kunstleve
    recently re-licensed to BSD),
  * sgtbx_ (successor of SgInfo written in C++/Python, part of cctbx),
  * CCP4 Clipper_,

  and by many programs.
  On the bad side, the conciseness is achieved by complex
  `rules <http://cci.lbl.gov/sginfo/hall_symbols.html>`_ of interpreting
  the symbols; the choice of a Hall symbol for given settings
  is not unambiguous and the symbols differ
  between editions of ITfC, and between sgtbx and :file:`syminfo.lib`.

After contemplating all the possibilities we ended up implementing
the most complex solution: Hall symbols. The relative complexity
does not mean it is slow: translating a Hall notation to generators
takes less than a microsecond on a typical desktop machine.
Closing a group is also below a microsecond for most of the groups,
and up to a few microseconds for the highest symmetry group Fm3̅m
(Gemmi uses Dimino's algorithm for this).

.. _SgInfo: https://github.com/rwgk/sginfo
.. _SgLite: https://github.com/rwgk/sglite
.. _sgtbx: https://github.com/rwgk/sglite
.. _csymlib: http://legacy.ccp4.ac.uk/html/C_library/csymlib_8h.html
.. _spglib: https://atztogo.github.io/spglib/
.. _Clipper: http://www.ysbl.york.ac.uk/~cowtan/clipper/doc/
.. _OpenBabel: https://github.com/openbabel/openbabel
.. _Mantid: https://github.com/mantidproject/mantid
.. _Shmueli: http://dx.doi.org/10.1107/S0108767384001161
.. _NGL: https://github.com/arose/ngl
.. _Fityk: https://github.com/wojdyr/fityk
.. _SPGGEN: http://dx.doi.org/10.1107/S1600576716007330
