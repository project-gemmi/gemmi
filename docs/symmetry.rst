Symmetry
########

The Gemmi symmetry module provides space group related functionality
needed in other parts of the library -- when working with coordinate
files, electron density maps and reflections.

Although the Gemmi project is developed for macromolecular crystallography,
for which only 65 space groups are relevant,
we cover all the 230 crystallographic space groups
for the sake of completeness.

For C++: this part of Gemmi has no dependencies,
all is in a single header :file:`symmetry.hpp`.

Space group table
=================

Gemmi tabulates 550+ settings of the 230 crystallographic space groups.
Each entry includes:

* ``number`` -- space group number (1-230),
* ``ccp4`` -- ccp4 number (assigned to particular settings; modulo 1000
  they give space group number: 3 and 1003 correspond to
  ``P 1 2 1`` and ``P 1 1 2``; 0 means none),
* ``hm`` -- Herman-Mauguin (H-M) symbol a.k.a. the international notation
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
including 3 duplicates (different names for the same settings)
in the space group 68.

Gemmi includes also settings from OpenBabel_ that are absent in
:file:`syminfo.lib`. If needed we will add more entries in the future.
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
Most of the crystallographic software (as well as the Volume B of the Tables)
still use the old names.
(spglib_ uses new ones,
sgtbx reads new names with the option ``ad_hoc_1992``).

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

.. note::

    The rest of this section has only Python examples mixed with the text.
    One longer C++ example is at the end.

The name above is expected to be either a full international Herman-Mauguin
symbol or a short symbol (*I2* instead of *I 1 2 1*).
This functions also searches the tabulated alternative names:

.. doctest::

    >>> gemmi.find_spacegroup_by_name('C m m e') # new names have 'e' and 'g'
    <gemmi.SpaceGroup("C m m a")>

Sometimes in the PDB, the setting of the hexagonal crystal system
is not clear from the symmetry symbol alone. For instance, the H-M symbol
"R 3" can mean either hexagonal or rhombohedral setting.
This ambiguity can be resolved by comparing angles of the unit cell.
The ratio of gamma to alpha angles is 120:90 in the hexagonal system
and 1:1 in rhombohedral. Therefore, ``find_spacegroup_by_name()``
accepts also alpha and gamma angles. If the angles are not passed,
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

The next function for searching the space group table checks
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

The last function for searching space group is also comparing operations.
It takes two arguments: a space group and a change-of-basis operator,
and searches for space group settings that match the transformed operations
of the original space group:

.. doctest::

  >>> # I2 -> C2
  >>> gemmi.find_spacegroup_by_change_of_basis(gemmi.SpaceGroup('I2'), gemmi.Op('x,y,x+z'))
  <gemmi.SpaceGroup("C 1 2 1")>
  >>> # enantiomorphic pair
  >>> gemmi.find_spacegroup_by_change_of_basis(gemmi.SpaceGroup('P 41'), gemmi.Op('-x,-y,-z'))
  <gemmi.SpaceGroup("P 43")>

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

``gemmi.SpaceGroup`` represents an entry in the space group table.
It has the properties listed at the beginning of this section
(``number``, ``ccp4``, ``hm``, ``ext``, ``hall``) and a few methods:

.. doctest::

  >>> sg = gemmi.SpaceGroup('R3') # equivalent to find_spacegroup_by_name()
  >>> sg.xhm()                    # extended Hermann-Mauguin name
  'R 3:H'
  >>> sg.short_name()             # short name
  'H3'
  >>> sg.is_enantiomorphic()      # is it one of 22 chiral space groups?
  False
  >>> sg.is_sohncke()             # is it one of 65 Sohncke space groups?
  True
  >>> sg.point_group_hm()         # H-M name of the point group
  '3'
  >>> sg.laue_str()               # name of the Laue class
  '-3'
  >>> sg.crystal_system_str()     # name of the crystal system
  'trigonal'
  >>> sg.is_reference_setting()
  True
  >>> # and the most important...
  >>> sg.operations()             #doctest: +ELLIPSIS
  <gemmi.GroupOps object at 0x...>

Categories related to chirality can be confusing.
Here, we follow the IUCr dictionary:

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

If you would like to ignore entries that are absent in SgInfo, sgtbx, spglib
and in the International Tables vol. B, use only the first 530 entries
of the gemmi table. In Python, we have a helper function for this:

.. doctest::

  >>> for sg in gemmi.spacegroup_table_itb():
  ...     pass
  >>> sg.ccp4  # sg here is the last space group that was iterated
  230


Implementation notes
--------------------

The choice between having an explicit list of all
the operations and generating them from a smaller set of symbols is a
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
.. _csymlib: http://www.ccp4.ac.uk/html/C_library/csymlib_8h.html
.. _spglib: https://atztogo.github.io/spglib/
.. _Clipper: http://www.ysbl.york.ac.uk/~cowtan/clipper/doc/
.. _OpenBabel: https://github.com/openbabel/openbabel
.. _Mantid: https://github.com/mantidproject/mantid
.. _Shmueli: http://dx.doi.org/10.1107/S0108767384001161
.. _NGL: https://github.com/arose/ngl
.. _Fityk: https://github.com/wojdyr/fityk
.. _SPGGEN: http://dx.doi.org/10.1107/S1600576716007330


Operations
==========

Crystallographic symmetry operations have a few notations.
Gemmi undertands only coordinate triplets (sometimes called
the Jones' faithful notation) such as ``x,x-y,z+1/2``.
The symmetry operation is represented by class Op.
It can be created in C++ as::

  #include <gemmi/symmetry.hpp>
  ...
  gemmi::Op op = gemmi::parse_triplet("-y,x-y,z");

and in Python as:

.. doctest::

    >>> import gemmi
    >>> op = gemmi.Op('-y,x-y,z+1/3')

We can also do the opposite:

.. doctest::

    >>> op.triplet()
    '-y,x-y,z+1/3'

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
    <gemmi.Op("-1/3*x+2/3*y-1/3*z,-2/3*x+1/3*y+1/3*z,1/3*x+1/3*y+1/3*z")>
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

  >>> ops = gemmi.find_spacegroup_by_name('I2').operations()
  >>> ops  #doctest: +ELLIPSIS 
  <gemmi.GroupOps object at 0x...>
  >>> list(ops.sym_ops)
  [<gemmi.Op("x,y,z")>, <gemmi.Op("-x,y,-z")>]
  >>> list(ops.cen_ops)
  [[0, 0, 0], [12, 12, 12]]

but they are be combined on the fly:

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

  >>> ops.change_basis(gemmi.Op('x,y,x+z'))  # I2 -> C2
  >>> gemmi.find_spacegroup_by_ops(ops)
  <gemmi.SpaceGroup("C 1 2 1")>

We can create GroupOps from a list of operations:

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

    >>> list(gemmi.generators_from_hall('P 4w 2c'))  # no.93
    [<gemmi.Op("x,y,z")>, <gemmi.Op("-y,x,z+1/4")>, <gemmi.Op("x,-y,-z+1/2")>]

Combining these 3 generators reconstructs all the 8 symmetry operations.

The GroupOps object has a couple of functions:

.. doctest::

  >>> new_ops.is_centric()
  False
  >>> new_ops.find_centering()
  'C'

and, again, it can be used to search in the space group table:

.. doctest::

  >>> gemmi.find_spacegroup_by_ops(new_ops)
  <gemmi.SpaceGroup("C 1 c 1")>

C++ Example
===========

Since the code snippets above were only in Python,
here we compensate it with one longer example in C++.

.. literalinclude:: code/sym.cpp
   :language: cpp

