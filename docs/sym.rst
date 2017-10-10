Symmetry
########

The Gemmi symmetry module provides space-group related functionality
needed in other parts of the library -- when working with coordinate
files, electron density maps and reflections.

Although the Gemmi project is developed for macromolecular crystallography,
for which only 65 space groups are relevant,
we cover all the 230 crystallographic space groups
for the sake of completeness.

This part of Gemmi has no dependencies:
all is in a single C++ header :file:`symmetry.hpp`.

Space group table
=================

Gemmi tabulates about 540 settings of the 230 crystallographic space groups,
including:

* space group numbers (1-230),
* ccp4 numbers (assigned to particular settings; modulo 1000
  they give space group number: 3 and 1003 correspond to
  ``P 1 2 1`` and ``P 1 1 2``; 0 means none),
* Herman-Mauguin (H-M) symbols a.k.a. the international notation
  (``I a -3 d``, ``C 1 2 1``),
* extensions to the H-M notations (none, ``1``, ``2``, ``H`` or ``R``)
  that make extended H-M symbols (``R 3:H``, ``P 4/n:1``),
* and the Hall symbols (``-I 4bd 2c 3``, ``C 2y``, ``P 32 2c (0 0 -1)``)
  that are used to generate operations.

This data is derived primarily from the CCP4 :file:`syminfo.lib` file,
which in turn is based on the data from sgtbx_ that was augmented
with the old data from a CCP4 file named :file:`symop.lib`.

The data from sgtbx is also available in the older SgInfo_ library,
as well as in the `International Tables <http://it.iucr.org/>`_
for Crystallography Vol. B ch. 1.4 (in 2010 ed.). It has 530 entries
including 3 duplicates (different names for the same settings)
in space-group 68.

We are aware that there are many other practical settings, and if needed we
will add more entries in the future. For example, we do not include
the C- and F-centred tetragonal space groups that are featured in
`Crystallographic Space Group Diagrams and Tables <http://img.chem.ucl.ac.uk/sgp/mainmenu.htm>`_
by Jeremy Karl Cockcroft
and which are also mentioned in the ITfC Vol.A (Table 1.5.4.4
in the 2015 edition).

We also tabulated alternative names.
For now this only includes new standard names introduced in 1990's by the IUCr
committee. For example, the space-group no. 39 is now officially, in the
Volume A of the `International Tables <http://it.iucr.org/>`_,
named Aem2 not Abm2.
Most of the crystallographic software (as well as the Volume B of the Tables)
still use the old names.
(spglib_ uses new ones,
sgtbx reads new names with the option ``ad_hoc_1992``).


Triplets and matrices
=====================

Crystallographic symmetry operations have a few notations.
Gemmi undertands only coordinate triplets (sometimes called
the Jones' faithful notation) such as ``x,x-y,z+1/2``.

The operations are equivalent to

* either a 3x3 rotation matrix and a translation
  vector,
* or to a single 4x4 transformation matrix (which in crystallography
  is called *Seitz matrix*).

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> import gemmi
    >>> gemmi.Op('x,-y,z+1/3').seitz()
    [[1, 0, 0, 0], [0, -1, 0, 0], [0, 0, 1, Fraction(1, 3)], [0, 0, 0, 1]]
    >>> gemmi.Op('x,-y,z+1/3').float_seitz()
    [[1.0, 0.0, 0.0, 0.0],
     [0.0, -1.0, 0.0, 0.0],
     [0.0, 0.0, 1.0, 0.3333333333333333],
     [0.0, 0.0, 0.0, 1.0]]

The triplets and matrices may not represent crystallographic symmetry,
but for example, generation of the biological assembly.
For this reason, ``x,y+1,z`` is not reduced to the identity.
But when the operations are transformed by the multiplication operator
in C++ and Python they are assumed to be symmetry operations
and the translation is wrapped to [0,1).

Operations and Hall symbols
===========================

Each space-group settings correspond to a unique set of operations.

For a software developer, the choice between having an explicit list of all
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
  This approach is used in Mantid_ (triplets are tabulated in
  :file:`SpaceGroupFactory.cpp`).

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

Operations generated from a Hall symbol are stored in an object
(``struct GroupOps``) with two lists: symmetry operations
(rotation + translation) and centring vectors (translation only).
``GroupOps`` has functions to iterate over all operations
(symops combined with centering ops), to apply change-of-basis,
to add missing group elements (if any), etc.

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


C++
===

.. literalinclude:: doc_sym.cpp
   :language: cpp

Python
======

.. doctest::

    >>> import gemmi
    >>>
    >>> # operators
    >>> gemmi.Op('x-y,x,z+1/6') * '-x,-y,z+1/2'
    <gemmi.Op("-x+y,-x,z+2/3")>
    >>> _.inverse()
    <gemmi.Op("-y,x-y,z-2/3")>
    >>> _.wrap()
    <gemmi.Op("-y,x-y,z+1/3")>
    >>>
    >>> # Hall symbols
    >>> list(gemmi.generators_from_hall("P 4w 2c"))  # no.93
    [<gemmi.Op("x,y,z")>, <gemmi.Op("-y,x,z+1/4")>, <gemmi.Op("x,-y,-z+1/2")>]
    >>> len(gemmi.symops_from_hall("P 4w 2c"))
    8
    >>>
    >>> # iterate over tabulated space-group settings
    >>> for sg in gemmi.spacegroup_table():
    ...   if sg.ccp4 != 0:
    ...     assert sg.ccp4 % 1000 == sg.number
    ...
    >>> max(sg.ccp4 for sg in gemmi.spacegroup_table())
    5005
    >>>
    >>> # select a space group
    >>> gemmi.find_spacegroup_by_number(5)
    <gemmi.SpaceGroup("C 1 2 1")>
    >>> gemmi.find_spacegroup_by_name('C m m e') # new names have 'e' and 'g'
    <gemmi.SpaceGroup("C m m a")>
    >>>
    >>> ops = gemmi.find_spacegroup_by_name('I2').operations()
    >>> for op in ops:
    ...   print(op.triplet())
    ...
    x,y,z
    -x,y,-z
    x+1/2,y+1/2,z+1/2
    -x+1/2,y+1/2,-z+1/2
    >>>
    >>> ops.change_basis(gemmi.Op('x,y,x+z'))  # I2 -> C2
    >>> gemmi.find_spacegroup_by_ops(ops)
    <gemmi.SpaceGroup("C 1 2 1")>

