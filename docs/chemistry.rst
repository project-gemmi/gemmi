
Chemistry
#########

This section covers:

* working with small molecule models in structural chemistry
* and working with chemical components in structural biology
  (chemical components describe parts of macromolecular models).

.. _elements:

Elements
========

When working with molecular structures, it is good to have basic data
from the periodic table at hand.

**C++**

.. literalinclude:: code/elem.cpp


**Python**

.. doctest::

    >>> import gemmi
    >>> gemmi.Element('Mg').weight
    24.305
    >>> gemmi.Element(118).name
    'Og'
    >>> gemmi.Element('Mo').atomic_number
    42

.. _covalent_radius:

We also included covalent radii of elements from a
`Wikipedia page <https://en.wikipedia.org/wiki/Covalent_radius>`_,
which has data from
Cordero *et al* (2008), *Covalent radii revisited*, Dalton Trans. 21, 2832.

.. doctest::

    >>> gemmi.Element('Zr').covalent_r
    1.75

van der Waals radii taken from Wikipedia and cctbx:

.. doctest::

    >>> gemmi.Element('K').vdw_r
    2.75

and a flag for metals:

.. doctest::

    >>> gemmi.Element('Mg').is_metal
    True
    >>> gemmi.Element('C').is_metal
    False

The classification into metals and non-metals is somewhat arbitrary.
It can be adjusted using the function `set_is_metal`:

.. doctest::

    >>> gemmi.Element('Sb').is_metal
    True
    >>> gemmi.set_is_metal('Sb', False)
    >>> gemmi.Element('Sb').is_metal
    False

The scattering properties of elements are covered
in the :ref:`Scattering <scattering>` section.

.. _small_molecules:

Small Molecules
===============

CIF files that describe small-molecule and inorganic structures
can be read into a `SmallStructure` object.
Unlike macromolecular :ref:`Structure <structure>`,
`SmallStructure` has no hierarchy.
It is a flat list of atomic sites (`SmallStructure::Site`)
together with the unit cell and symmetry.

.. literalinclude:: code/smcif.cpp

.. doctest::

    >>> import gemmi
    >>> SiC = gemmi.read_small_structure('../tests/1011031.cif')
    >>> SiC.cell
    <gemmi.UnitCell(4.358, 4.358, 4.358, 90, 90, 90)>
    >>> # content of _symmetry_space_group_name_H-M or _space_group_name_H-M_alt
    >>> SiC.spacegroup
    <gemmi.SpaceGroup("F -4 3 m")>
    >>> list(SiC.sites)
    [<gemmi.SmallStructure.Site Si1>, <gemmi.SmallStructure.Site C1>]
    >>> len(SiC.get_all_unit_cell_sites())
    8

Each atomic site has the following properties:

.. doctest::

    >>> site = SiC.sites[0]
    >>> site.label
    'Si1'
    >>> site.type_symbol
    'Si4+'
    >>> site.fract
    <gemmi.Fractional(0, 0, 0)>
    >>> site.occ
    1.0
    >>> site.u_iso  # not specified here
    0.0
    >>> site.element  # obtained from type_symbol 'Si4+'
    gemmi.Element('Si')
    >>> site.charge   # obtained from type_symbol 'Si4+'
    4

The occupancies in small molecules normally represent the actual chemical
occupancy.
This differs from macromolecular crystallography, where models normally store
"crystallographic" occupancy -- atoms on special positions have their occupancy
divided by the number of symmetry images in the same place.
This reduction of occupancy simplifies the calculation of structure factors.

.. doctest::

    >>> 1 / site.occ
    1.0
    >>> SiC.change_occupancies_to_crystallographic()
    >>> 1 / site.occ
    24.0

We will need another cif file to show anisotropic ADPs and disorder_group:

.. doctest::

    >>> perovskite = gemmi.read_small_structure('../tests/4003024.cif')
    >>> for site in perovskite.sites:
    ...   print(site.label, site.aniso.nonzero(), site.disorder_group or 'n/a')
    Cs1 True n/a
    Sn2 False 1
    Cl1 True n/a
    In False 2
    >>> perovskite.sites[2].aniso.u11
    0.103
    >>> perovskite.sites[2].aniso.u22
    0.156
    >>> perovskite.sites[2].aniso.u33
    0.156
    >>> perovskite.sites[2].aniso.u12
    0.0
    >>> perovskite.sites[2].aniso.u13
    0.0
    >>> perovskite.sites[2].aniso.u23
    0.0

----

The Python examples above read CIF files using `read_small_structure()`.
Alternatively, the same can be done in two steps:

.. doctest::

    >>> cif_doc = gemmi.cif.read('../tests/1011031.cif')
    >>> SiC = gemmi.make_small_structure_from_block(cif_doc.sole_block())

Now you also have access to the CIF document.

.. _small_spacegroup:

SmallStructure::spacegroup
--------------------------

When reading a small-molecule CIF file, a few CIF items that describe
the space group are read and stored in member variables:

.. doctest::

    >>> st = gemmi.read_small_structure('../tests/2013551.cif')
    >>> st.symops
    ['x, y, z', '-y, x-y, z', 'y, x, -z', '-x+y, -x, z', '-x, -x+y, -z', 'x-y, -y, -z', '-x, -y, -z', 'y, -x+y, -z', '-y, -x, z', 'x-y, x, -z', 'x, x-y, z', '-x+y, y, z']
    >>> st.spacegroup_hall
    '-P 3 2"'
    >>> st.spacegroup_hm
    'P -3 m 1'
    >>> st.spacegroup_number
    164

and the function `determine_and_set_spacegroup("S.H2")` is automatically
run to set `spacegroup`:

.. doctest::

    >>> st.spacegroup
    <gemmi.SpaceGroup("P -3 m 1")>

`determine_and_set_spacegroup()` takes one argument, a string in which characters
specify what to use, and in what order, for space group determination:

* `S` = symmetry operations stored in `symops`,
* `H` = Hall symbol from `spacegroup_hall` (we compare symmetry operations
  encoded in the Hall symbol, not the strings),
* `1` = H-M symbol; for space groups such as "P n n n" that have two origin
  choices listed in the International Tables, use *Origin Choice 1*,
* `2` = H-M symbol, with *Origin Choice 2* where applicable,
* `N` = the space group number,
* `.` (after S or H) = if the symmetry operations pass sanity checks,
  stop and use them regardless of whether they correspond to one of
  the settings tabulated in Gemmi.

If a symbol or operations match one of the 560+ space group settings tabulated
in Gemmi, `spacegroup` is set to this setting. Otherwise, if `.` is encountered
and the previous character (`S` or `H`) was evaluated to a valid set of symops,
it is assumed that these operations were correct: `spacegroup` is left null
and `cell.images` are set from the list of operations.
About 350 (out of 500,000+) entries in the COD use such settings.
Most of them have an unconventional choice of the origin
(e.g. "P 1 21 1 (a,b,c-1/4)").

To use a different order of items than "S.H2",
call determine_and_set_spacegroup() again:

.. doctest::

    >>> st.determine_and_set_spacegroup('H.1')

Errors such as an incorrect format of the symop triplets or of the Hall
symbol are silently ignored, and the consistency between different items
is not checked. That's because this function is run when reading a file;
throwing an exception at that stage would prevent reading a file.
We have a separate function to check for errors and inconsistencies.
It returns a string, one line -- one error:

.. doctest::

    >>> st.check_spacegroup()
    ''

If the spacegroup setting used in a file is not tabulated in Gemmi,
you can still create a GroupOps object with symmetry operations:

.. doctest::

    >>> gemmi.GroupOps([gemmi.Op(o) for o in st.symops])  #doctest: +ELLIPSIS
    <gemmi.GroupOps object at 0x...>
    >>> # or
    >>> gemmi.symops_from_hall(st.spacegroup_hall)  #doctest: +ELLIPSIS
    <gemmi.GroupOps object at 0x...>

In C++ it would be similar, except that the following function
would be used to make gemmi::GroupOps from symops::

    GroupOps split_centering_vectors(const std::vector<Op>& ops)


without CIF file
----------------

If your structure is stored in a macromolecular format (PDB, mmCIF)
you can read it first as macromolecular :ref:`Structure <structure>`
and convert it to `SmallStructure`:

.. doctest::

  >>> gemmi.mx_to_sx_structure(gemmi.read_structure('../tests/HEM.pdb'))
  <gemmi.SmallStructure: HEM>

You could also create `SmallStructure` from scratch:

.. doctest::

    >>> small = gemmi.SmallStructure()
    >>> small.spacegroup_hm = 'F -4 3 m'
    >>> small.cell = gemmi.UnitCell(4.358, 4.358, 4.358, 90, 90, 90)
    >>> small.determine_and_set_spacegroup("2")
    >>> # add a single atom
    >>> site = gemmi.SmallStructure.Site()
    >>> site.label = 'C1'
    >>> site.element = gemmi.Element('C')
    >>> site.fract = gemmi.Fractional(0.25, 0.25, 0.25)
    >>> site.occ = 1
    >>> small.add_site(site)


.. _chemcomp:

Chemical Components
===================

Residues (monomers) and small molecule components of macromolecular models
are called *chemical components*.
Gemmi can use three sources of knowledge about chemical components:

* built-in basic data about 350+ popular components,
* the Chemical Component Dictionary (CCD) maintained by the PDB
  (25,000+ components),
* so-called CIF files compatible with the format of the CCP4 Monomer Library
  (more about monomer libraries in the :ref:`next section <monlib>`).

.. _find_tabulated_residue:

Built-in data
-------------

The built-in data is accessed through the function `find_tabulated_residue`.
It contains only minimal information about each residue:
assigned category, the "standard" flag (non-standard residues are marked
as HETATM in the PDB, even in polymer), one-letter code,
the number of hydrogens and molecular weight:

.. literalinclude:: code/resinfo.cpp

.. doctest::

    >>> gln = gemmi.find_tabulated_residue('GLN')
    >>> gln.is_amino_acid()
    True
    >>> gln.one_letter_code
    'Q'
    >>> round(gln.weight, 3)
    146.144
    >>> gln.hydrogen_count
    10
    >>> gemmi.find_tabulated_residue('DOD').is_water()
    True
    >>> # PDB marks "non-standard" residues as HETATM.
    >>> # Pyrrolysine is standard - some microbes have it.
    >>> gemmi.find_tabulated_residue('PYL').is_standard()
    True
    >>> gemmi.find_tabulated_residue('MSE').is_standard()
    False

One-letter code is an upper case letter if it is a standard residue.
Otherwise, it can be the letter for the parent residue in lower case,
or a space. It is common to use `X` for non-standard residue --
for this we have helper function `fasta_code()`:

.. doctest::

    >>> gemmi.find_tabulated_residue('MET').one_letter_code
    'M'
    >>> gemmi.find_tabulated_residue('MSE').one_letter_code
    'm'
    >>> gemmi.find_tabulated_residue('HOH').one_letter_code
    ' '
    >>> gemmi.find_tabulated_residue('MET').fasta_code()
    'M'
    >>> gemmi.find_tabulated_residue('MSE').fasta_code()
    'X'

The table includes only 362 entries, selected from the most popular residues in the PDB.
Residue kind is sometimes debatable, the user may change it.

.. doctest::

    >>> pst = gemmi.find_tabulated_residue('PST')
    >>> pst.kind
    ResidueKind.UNKNOWN
    >>> pst.kind = gemmi.ResidueKind.DNA

.. _CCD_etc:

CCD and monomer CIF files
-------------------------

To get more complete information, including atoms and bonds in the monomer,
we need to first read either the `CCD <https://www.wwpdb.org/data/ccd>`_
or a :ref:`monomer library <monlib>`.

The CCD :file:`components.cif` file describes all the monomers
(residues, ligands, solvent molecules) from the PDB entries.
Importantly, it contains information about bonds.

.. note::

    The absence of bond information in mmCIF files from wwPDB is a
    `well-known problem <https://www.cgl.ucsf.edu/chimera/data/mmcif-oct2013/mmcif.html>`_.
    This information is included in so-called
    `updated mmCIF <https://doi.org/10.1093/nar/gkv1047>`_ files from PDBe,
    as well as in BinaryCIF and mmJSON files.

Macromolecular refinement programs need to know more about monomers
than the CCD can tell: they need to know how to restrain the structure.
Therefore, they have their own dictionaries of monomers (a.k.a monomer
libraries), such as the CCP4 Monomer Library (for Refmac),
where each monomer is described by one cif file.
These libraries are often complemented by user's own cif files.

Gemmi provides a `ChemComp` class that corresponds to a monomer
from either the CCD or a cif file.

.. literalinclude:: ../examples/with_bgl.cpp
   :lines: 13-14,42-47

.. doctest::

    >>> # SO3.cif -> gemmi.ChemComp
    >>> block = gemmi.cif.read('../tests/SO3.cif')[-1]
    >>> so3 = gemmi.make_chemcomp_from_block(block)

This class is not fully documented yet.

The examples in :ref:`graph_analysis`
show how to access `ChemComp`'s atoms and bonds.

.. _monlib:

Monomer library
===============

Structural biologists routinely use prior knowledge about biomolecules
to augment the data obtained in an experiment.
This prior knowledge is what we know about preferred geometries in molecules
(distances between atoms, etc.). This knowledge is extracted primarily from
experimental small molecule databases (COD and CSD) and QM calculations.
One way to store that prior knowledge is in a so-called *monomer library*.
In addition to describing monomers (chemical components
from the previous section), the monomer library also describes links
between monomers and may contain various other data useful in
macromolecular refinement.

In Gemmi, data from a monomer library is stored in the class `MonLib`.
Currently, `MonLib` is modeled after and works only with
the `CCP4 monomer library <https://github.com/MonomerLibrary/monomers>`_.
This library was introduced in the
`early 2000s <https://doi.org/10.1107/S0907444904023510>`_
to provide restraint templates for Refmac.
There are only two other popular MX refinement programs: PHENIX and BUSTER.
PHENIX provides `geostd <https://github.com/MonomerLibrary/monomers>`_,
which was forked from CCP4 ML and is still quite similar.
In `BUSTER <https://www.globalphasing.com/buster/>`_ the prior knowledge
is organized differently.

The restraints we use are similar to those used in molecular dynamics
(bond, angle, dihedral, improper dihedral and planarity restraints).
Originally, the monomer library was created because MD potentials
were deemed inadequate for refinement. Since then,
both the restraints in experimental structural biology and MD potentials
have improved, independently of each other. In recent years there have been
a few examples of using AMBER and OpenMM potentials for MX refinement.
Currently, there is no clear advantage to one approach over the other.

`MonLib` has the following member variables (`std::` omitted for readability):

* `string monomer_dir` -- the top-level directory, in CCP4 it's `$CLIBD_MON`.
* `map<string, ChemComp> monomers` -- chemical components read from
  the monomer library. Usually, we read only the components present in a model.
  Chemical components contains restraint templates
* `map<string, ChemLink> links` --  link descriptions. Each ChemLink
  contains rules for determining to what links it is applicable,
  restrain templates for restraining a link, and the names of modifications
  to be applied to linked monomers.
* `map<string, ChemMod> modifications` -- each modification is a set of rules
  for changing `ChemComp`. The rules can add, remove or modify
  atoms (they often remove a hydrogen atom) and restraints.
* `map<string, ChemComp::Group> cc_groups` -- groups defined in the library
  to classify monomers. Each monomer is assigned a group (such as `peptide`
  or `non-polymer`). The groups are then used to match link templates against
  the links between monomers.
* `EnerLib ener_lib` -- data from `$CLIBD_MON/ener_lib.cif`.

TBC

For now, here is an example of how to read the CCP4 monomer library
(Refmac dictionary):

.. doctest::
  :skipif: ccp4_path is None

  >>> monlib_path = os.environ['CCP4'] + '/lib/data/monomers'
  >>> # Usually, residue names from a model are obtained by calling:
  >>> #resnames = gemmi.Model.get_all_residue_names()
  >>> resnames = ['LEU', 'VAL', 'HIS', 'SER', 'ASN', 'HOH']
  >>> monlib = gemmi.MonLib()
  >>> monlib.read_monomer_lib(monlib_path, resnames, logging=sys.stderr)
  True

The `logging` argument above is described in the next section.

`MonLib` can be used to prepare :ref:`Topology <topology>`.

TBC

