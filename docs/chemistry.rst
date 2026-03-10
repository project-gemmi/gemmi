
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

.. _chemcomp-smarts:

SMARTS Matching
---------------

Gemmi includes a lightweight SMARTS-subset matcher that can be used to find
functional groups or specific motifs within a `ChemComp`.

It supports:

* Atomic symbols (e.g. `C`, `[N]`, `[Fe]`) and wildcards (`*`).
* Aromaticity (`[c]`, `[n]`, etc.).
* Bond types: single (`-`), double (`=`), and any (`~`).
* Connectivity and H-count constraints (e.g. `[CX4]`, `[OH1]`).
* Branching using parentheses.

.. doctest::

    >>> # Match aromatic carbons in benzene
    >>> benzene_cif = "data_comp_BEN\nloop_\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n" + \
    ...               "\n".join(f"C{i} C" for i in range(1,7)) + \
    ...               "\nloop_\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\n" + \
    ...               "\n".join(f"C{i} C{i%6+1} arom" for i in range(1,7))
    >>> cc = gemmi.make_chemcomp_from_block(gemmi.cif.read_string(benzene_cif).sole_block())
    >>> matches = cc.match_smarts("[c]")
    >>> len(matches)
    6

.. _chemcomp-chemical-normalization:

Chemical normalization
======================

Gemmi's `ChemComp` normalization is performed by
`apply_chemical_adjustments()` (declared in `gemmi/cc_adj.hpp`,
implemented in `src/cc_adj.cpp`).

This is a deterministic local-graph normalization stage that can change:

* formal charges,
* protonation state (add/remove hydrogens),
* selected bond orders (for resonance-style normalization),
* affected restraints around edited atoms (bonds/angles tied to added or
  removed hydrogens).

It is intentionally rule-based and motif-driven; it is not a general pKa
predictor or tautomer enumerator.

apply_chemical_adjustments()
---------------------

Rules are applied in a fixed order (inside of `apply_chemical_adjustments()`):
examples below show representative components.

Acid/oxoacid deprotonation
^^^^^^^^^^^^^^^^^^^^^^^^^^

`oxoacid_phosphate` (step 1, SMARTS: ``OP(=O)(O)(*)`` / ``O=P(O)(O)O``)

   .. list-table::
      :widths: 45 10 45
      :class: borderless

      * - .. figure:: img/adj_oxoacid_phosphate_before.svg
            :alt: oxoacid_phosphate before
            :width: 100%
        - ➡
        - .. figure:: img/adj_oxoacid_phosphate_after.svg
            :alt: oxoacid_phosphate after
            :width: 100%

   Example: `ATP <https://www.rcsb.org/ligand/ATP>`_
   For phosphate motifs of the form R-O-PO3, this rule deprotonates
   phosphoryl oxygens (AceDRG-style doubly deprotonated representation).
   The alkoxy R-O part is not deprotonated by this rule.

`oxoacid_sulfate` (step 2, SMARTS: ``O=S(=O)(O)O``)

   .. list-table::
      :widths: 45 10 45
      :class: borderless

      * - .. figure:: img/adj_oxoacid_sulfate_before.svg
            :alt: oxoacid_sulfate before
            :width: 100%
        - ➡
        - .. figure:: img/adj_oxoacid_sulfate_after.svg
            :alt: oxoacid_sulfate after
            :width: 100%

   Example: `0SG <https://www.rcsb.org/ligand/0SG>`_

`oxoacid_sulfite` (step 3, SMARTS: ``OS(=O)O``)

   .. list-table::
      :widths: 45 10 45
      :class: borderless

      * - .. figure:: img/adj_oxoacid_sulfite_before.svg
            :alt: oxoacid_sulfite before
            :width: 100%
        - ➡
        - .. figure:: img/adj_oxoacid_sulfite_after.svg
            :alt: oxoacid_sulfite after
            :width: 100%

   This rule currently follows AceDRG-style sulfite handling by deprotonating
   one hydroxyl oxygen in matching sulfate-like tri-oxygen sulfur motifs.

`single_bond_oxide` (step 5, SMARTS: ``[O]-[*]`` + internal constraints)

   .. list-table::
      :widths: 45 10 45
      :class: borderless

      * - .. figure:: img/adj_single_bond_oxide_before.svg
            :alt: single_bond_oxide before
            :width: 100%
        - ➡
        - .. figure:: img/adj_single_bond_oxide_after.svg
            :alt: single_bond_oxide after
            :width: 100%

   Example: `H1T <https://www.rcsb.org/ligand/H1T>`_
   Oxygen must be neutral, have no hydrogen neighbors, and have exactly one
   heavy single-bond neighbor.
   This is not general alcohol deprotonation. It normalizes pre-existing
   single-bond oxide-like oxygens that are already non-protonated in the graph.

`carboxy_asp` (step 7, SMARTS: ``O=C(O)C(*)`` / ``O=C(O)c(*)`` / ``O=C(O)CN(*)`` / ``O=C(O)C(N)(*)`` / ``O=C(O)CN``)

   .. list-table::
      :widths: 45 10 45
      :class: borderless

      * - .. figure:: img/adj_carboxy_asp_before.svg
            :alt: carboxy_asp before
            :width: 100%
        - ➡
        - .. figure:: img/adj_carboxy_asp_after.svg
            :alt: carboxy_asp after
            :width: 100%

   Example: `ASP <https://www.rcsb.org/ligand/ASP>`_

`terminal_carboxylate` (step 8, SMARTS: ``O=C(O)CN(*)`` with `OXT`/`HXT` naming context)

   .. list-table::
      :widths: 45 10 45
      :class: borderless

      * - .. figure:: img/adj_terminal_carboxylate_before.svg
            :alt: terminal_carboxylate before
            :width: 100%
        - ➡
        - .. figure:: img/adj_terminal_carboxylate_after.svg
            :alt: terminal_carboxylate after
            :width: 100%

   Example: `A0G <https://www.rcsb.org/ligand/A0G>`_
   Unlike `carboxy_asp`, this rule targets the terminal carboxylate motif
   identified via `OXT`/`HXT` context.

Resonance normalization
^^^^^^^^^^^^^^^^^^^^^^^

`nitro_group` (step 4, SMARTS: ``[N+](=O)(O)C``)

   .. list-table::
      :widths: 45 10 45
      :class: borderless

      * - .. figure:: img/adj_nitro_group_before.svg
            :alt: nitro_group before
            :width: 100%
        - ➡
        - .. figure:: img/adj_nitro_group_after.svg
            :alt: nitro_group after
            :width: 100%

   Example: `NE5 <https://www.rcsb.org/ligand/NE5>`_
   Shown as charge-separated nitro resonance (`R-N+(=O)-O-`).
   The exact AceDRG representation (formal charges vs valence-only
   representation) is under review.

Targeted handling
^^^^^^^^^^^^^^^^^

`hexafluorophosphate` (step 6, SMARTS: ``F[P](F)(F)(F)(F)F``)

   .. list-table::
      :widths: 45 10 45
      :class: borderless

      * - .. figure:: img/adj_hexafluorophosphate_before.svg
            :alt: hexafluorophosphate before
            :width: 100%
        - ➡
        - .. figure:: img/adj_hexafluorophosphate_after.svg
            :alt: hexafluorophosphate after
            :width: 100%

   Example: `A9J <https://www.rcsb.org/ligand/A9J>`_

Cationic nitrogen completion/protonation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`guanidinium` (step 9, SMARTS: ``CNC(=N)N``)

   .. list-table::
      :widths: 45 10 45
      :class: borderless

      * - .. figure:: img/adj_guanidinium_before.svg
            :alt: guanidinium before
            :width: 100%
        - ➡
        - .. figure:: img/adj_guanidinium_after.svg
            :alt: guanidinium after
            :width: 100%

   Example: `00L <https://www.rcsb.org/ligand/00L>`_

`amino_ter_amine` (step 10, SMARTS: ``NCC(=O)N(*)``)

   .. list-table::
      :widths: 45 10 45
      :class: borderless

      * - .. figure:: img/adj_amino_ter_amine_before.svg
            :alt: amino_ter_amine before
            :width: 100%
        - ➡
        - .. figure:: img/adj_amino_ter_amine_after.svg
            :alt: amino_ter_amine after
            :width: 100%

   Example: `00K <https://www.rcsb.org/ligand/00K>`_

`terminal_amine` (step 11, SMARTS: ``[*]C[NH2]``)

    .. list-table::
       :widths: 45 10 45
       :class: borderless

       * - .. figure:: img/adj_terminal_amine_before.svg
             :alt: terminal_amine before
             :width: 100%
         - ➡
         - .. figure:: img/adj_terminal_amine_after.svg
             :alt: terminal_amine after
             :width: 100%

   Example: `LYS <https://www.rcsb.org/ligand/LYS>`_
   Current implementation applies context-dependent protonation of terminal
   amine-like motifs.

`protonated_amide_n` (step 12, SMARTS: ``CCC(=O)N``)

    .. list-table::
       :widths: 45 10 45
       :class: borderless

       * - .. figure:: img/adj_protonated_amide_n_before.svg
             :alt: protonated_amide_n before
             :width: 100%
         - ➡
         - .. figure:: img/adj_protonated_amide_n_after.svg
             :alt: protonated_amide_n after
             :width: 100%

    Example: `BJS <https://www.rcsb.org/ligand/BJS>`_

The order is part of behavior: earlier edits can affect pattern matching in
later steps.

prepare_chemcomp()
------------------

`prepare_chemcomp()` is the full restraint-preparation pipeline. The order
below is important because later stages depend on graph/charge edits made
earlier.

The function now also has an explicit options object,
`PrepareChemcompOptions`, which makes the runtime controls part of the API
instead of relying only on environment variables. The legacy overload is still
available.

Execution order
^^^^^^^^^^^^^^^

`prepare_chemcomp()` follows this control flow:

0. detect carborane branches before the standard path:

   * full carborane mode (early return): enabled when the component has at
     least one non-hydrogen atom with 4+ boron neighbors, and all
     non-hydrogen atoms are only B/C/metal. In this case,
     `apply_carborane_mode()` is used, CCP4 types are assigned, and the
     function returns immediately.
     Example: `1KW <https://www.rcsb.org/ligand/1KW>`_.
   * otherwise continue with the standard path, but keep a carborane-seed flag
     for possible mixed-mode post-processing.
     Example: `9UK <https://www.rcsb.org/ligand/9UK>`_.

1. seed missing angles from existing bonds (unless `--no-angles`);
2. run `apply_chemical_adjustments()`;
3. run `add_n_terminal_h3()` (may add `H3` and corresponding N-centered angles);

   `add_n_terminal_h3()` adds the third proton only for matching
   N-terminus-like motifs. `sync_n_terminal_h3_angles()` is run later (after
   table fill) to keep the newly added H3 angles consistent with sibling
   N-H/N-H2 restraints.

   .. list-table::
      :widths: 45 10 45
      :class: borderless

      * - .. figure:: img/adj_add_n_terminal_h3_before.svg
            :width: 100%
        - ➡
        - .. figure:: img/adj_add_n_terminal_h3_after.svg
            :alt: add_n_terminal_h3 after
            :width: 100%

   Example: `ALA <https://www.rcsb.org/ligand/ALA>`_
4. run charge corrections;
5. fill missing bond/angle values from AceDRG tables;
6. if `H3` was added, run `sync_n_terminal_h3_angles()` to align H3-angle
   values with existing N-H/N-H2 geometry;
7. add missing chirality, torsion and plane restraints;
8. if the carborane-seed flag is set, apply mixed carborane post-processing:
   cluster-local part.
   Example: `9UK <https://www.rcsb.org/ligand/9UK>`_.
9. assign CCP4 atom types.

Python example:

.. doctest::
  :skipif: ccp4_path is None

    >>> import os
    >>> path = '../tests/ccd/ASP.cif'
    >>> block = gemmi.cif.read(path).sole_block()
    >>> cc = gemmi.make_chemcomp_from_block(block)
    >>> tables = gemmi.AcedrgTables()
    >>> tables.load_tables(ccp4_path + '/share/acedrg/tables')
    >>> before = {a.id for a in cc.atoms}
    >>> gemmi.prepare_chemcomp(cc, tables)
    >>> after = {a.id for a in cc.atoms}
    >>> removed = sorted(before - after)
    >>> added = sorted(after - before)
    >>> atoms = {a.id: a for a in cc.atoms}
    >>> removed
    ['HD2', 'HXT']
    >>> added
    ['H3']
    >>> atoms['N'].charge
    1.0
    >>> atoms['OD2'].charge
    -1.0
    >>> atoms['OXT'].charge
    -1.0

Explicit options
^^^^^^^^^^^^^^^^

`PrepareChemcompOptions` currently exposes:

* `atom_stereo` -- maps atom names to `pdbx_stereo_config` strings used by
  chirality generation,
* `no_angles` -- skips the initial angle-seeding pass,
* `strict_mode`, `compat_mode`, `trace_mode` -- tri-state overrides with
  values `PrepareOverride.Auto`, `PrepareOverride.Disable`, or
  `PrepareOverride.Enable`.

In `Auto` mode, the corresponding environment variable is still consulted
(`GEMMI_ACE_STRICT`, `GEMMI_ACE_COMPAT`, `GEMMI_ACE_TRACE`). `Enable` and
`Disable` override the environment explicitly for that call.

Mode semantics
""""""""""""""

`strict_mode`
    Enables extra validation checks during the pipeline. At present, this
    mainly checks for invalid restraint references, duplicate atom ids in
    planes, and NaN bond/angle targets at the final stage. In strict mode,
    such conditions raise an exception instead of being tolerated or silently
    cleaned up.

`compat_mode`
    Enables AceDRG-like behavior in places where Gemmi now has a more
    chemistry-first default.

    Current effects:

    * peptide mode is gated by the canonical-AA-backbone test
      (`N`, `CA`, `C`, `O`, `OXT` and associated bonding / H naming),
    * torsion candidate selection uses the older AceDRG-like selector instead
      of Gemmi's scored selector,
    * pyranose chair torsion rewriting is always applied, instead of only when
      coordinates are incomplete.

    With `compat_mode` disabled, Gemmi prefers chemically reasonable defaults
    over matching AceDRG quirks.

`trace_mode`
    Emits phase-by-phase summaries to stderr. The trace reports how many atoms,
    bonds, angles, torsions, chiralities and planes were added or removed in
    each pipeline phase, plus net charge changes. This is intended for
    debugging and regression analysis rather than normal use.

Python example:

.. doctest::
  :skipif: ccp4_path is None

    >>> import gemmi
    >>> path = '../tests/ccd/ASP.cif'
    >>> block = gemmi.cif.read(path).sole_block()
    >>> cc = gemmi.make_chemcomp_from_block(block)
    >>> tables = gemmi.AcedrgTables()
    >>> tables.load_tables(ccp4_path + '/share/acedrg/tables')
    >>> options = gemmi.PrepareChemcompOptions()
    >>> options.compat_mode = gemmi.PrepareOverride.Disable
    >>> options.strict_mode = gemmi.PrepareOverride.Disable
    >>> gemmi.prepare_chemcomp(cc, tables, options)

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

`read_monomer_lib()` takes a monomer directory path, a list of residue names,
and an optional :ref:`Logger <logger>` argument.
Example -- reading the CCP4 monomer library (Refmac dictionary):

.. doctest::
  :skipif: ccp4_path is None

  >>> monlib_path = os.environ['CCP4'] + '/lib/data/monomers'
  >>> # Usually, residue names from a model are obtained by calling:
  >>> #resnames = gemmi.Model.get_all_residue_names()
  >>> resnames = ['LEU', 'VAL', 'HIS', 'SER', 'ASN', 'HOH']
  >>> monlib = gemmi.MonLib()
  >>> monlib.read_monomer_lib(monlib_path, resnames, logging=sys.stderr)
  True

EnerLib
-------

`EnerLib` is the in-memory representation of selected data from
`ener_lib.cif` (the same schema used in CCP4 monomer-library and
AceDRG-tables layouts).
It is attached to `MonLib` as `MonLib.ener_lib`.

The atom types here are :ref:`CCP4 energy types <energy-types>`
described further below.

Loaded categories:

* `_lib_atom` -> atom-type properties (`element`, `hb_type`, `vdw_radius`,
  `vdwh_radius`, `ion_radius`, `valency`, `sp`),
* `_lib_bond` -> bond-type entries (`atom_type_1`, `atom_type_2`, bond order,
  length, esd).

How it is used:

* `MonLib.read_monomer_lib()` reads `ener_lib.cif` and populates `ener_lib`,
* `MonLib.find_ideal_distance()` consults it to estimate ideal distances
  between atoms, including metal-aware radius handling and fallback tiers
  when exact type-pair entries are missing.

Reading only `ener_lib.cif` into `EnerLib` (Python):

.. doctest::
  :skipif: ccp4_path is None

  >>> ener = gemmi.EnerLib()
  >>> ener.read(os.path.join(monlib_path, 'ener_lib.cif'))
  >>> len(ener.atoms) > 0
  True

`MonLib` can be used to prepare :ref:`Topology <topology>`.

.. _chemistry-gemmi-drg-overview:

gemmi drg: high-level overview
===============================

This section is a conceptual overview of `gemmi drg` capabilities.
API-level documentation will be added separately.

`gemmi drg` generates monomer restraint dictionaries by combining
chemical rules with statistical knowledge derived from AceDRG tables.

Workflow and data sources
-------------------------

At a high level, `gemmi drg`:

* reads a monomer definition (typically a CCD-like CIF),
* builds a molecular graph (atoms, bonds, valence context),
* consults AceDRG-derived tables for bond/angle targets and sigmas,
* emits restraint categories suitable for refinement workflows.

Input expectations and normalization
------------------------------------

Input is typically a CCD-style component definition with atom and bond
information. In practice, `gemmi drg` must also normalize chemistry before
table lookup, because small differences in protonation state or local bond
annotation can move an atom into a different type bucket.

For details on the ChemComp-level rule set, see
:ref:`chemcomp-chemical-normalization`.

Normalization includes:

* protonation/deprotonation rules aligned with AceDRG behavior where
  feasible,
* selected functional-group corrections used to stabilize atom typing and
  avoid type drift from equivalent input representations,
* graph cleanup steps that make downstream typing and fallback selection
  deterministic.

These functional-group corrections are intentionally pragmatic rather than
fully general pKa modeling. Current examples include:

* oxoacid normalization (for phosphate/sulfate-like motifs), including
  deterministic O-H deprotonation in qualifying local patterns,
* carboxylate normalization (side-chain and terminal forms), including
  removal of acidic hydrogens where AceDRG would represent a deprotonated
  carboxylate,
* amine/guanidinium normalization, including AceDRG-style protonation and
  hydrogen completion for selected terminal or resonance-stabilized motifs,
* selected special cases such as PF6-like phosphorus environments and
  metal-adjacent non-metal charge correction used by downstream typing.

Core pipeline
-------------

Given a monomer CIF, the current implementation performs these stages:

* chemistry normalization, including protonation/deprotonation handling
  and selected functional-group corrections,
* atom-environment derivation and atom typing (both AceDRG-style
  signatures and CCP4-compatible energy types),
* bond and angle assignment from AceDRG-style reference statistics,
  with progressively broader fallback levels,
* inference of missing stereochemical categories from topology, including
  torsions, chiral centers and planarity restraints.

Restraint assignment strategy
-----------------------------

Bond/angle assignment is driven by hierarchical matching. Conceptually, it
tries the most specific local environment first, then relaxes constraints
in controlled steps until a statistically supported target is found.

Bond-table matching levels
^^^^^^^^^^^^^^^^^^^^^^^^^^

For bonds, matching is keyed by progressively simpler context:

* atom hashes, hybridization pair and same-ring flag,
* neighbor descriptors (`nb2` and `nb1nb2`-style summaries),
* atom-type descriptors (`cod_main`) and full COD class for exact matching.

The level sequence mirrors AceDRG-style logic:

* level 0: exact full COD-class match (most specific),
* levels 1-2: type-relaxed matches/aggregates,
* levels 3-6: neighborhood-relaxed matches/aggregates,
* levels 7-8: broader `nb2`-level aggregation,
* levels 9-11: hash/hybrid/ring summaries from HRS-style hierarchy.

The search does not always start at level 0. It computes a dynamic start
level from key availability, skipping levels that cannot match for the
current atom pair.

Angle-table matching levels
^^^^^^^^^^^^^^^^^^^^^^^^^^^

For angles, keys are centered on the middle atom hash plus sorted flank
hashes, ring-size/hybridization value key, and then progressively reduced
context (roots, neighbor summaries, atom types).

The lookup ladder is:

* 1D: full key including types (approx level 0),
* 2D: no type component (approx 1),
* 3D: no `nb` component (approx 2),
* 4D: roots-only beyond hash/value key (approx 3),
* 5D: hash + value key (approx 4),
* 6D: hash-only summary (approx 6).

When center/flank hashes are equal, an additional swapped-center pass is
used to reproduce AceDRG table-orientation behavior. There is also a
wildcard partial-hash fallback used only when the exact hash triple has
no loaded table entries.

Thresholding and final fallbacks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Acceptance is thresholded by observation support (defaults are AceDRG-like;
for many paths the default minimum is 3 observations).

In particular:

* bond level 0 uses observation count from the matched record,
* bond levels 1-8 use the number of contributing entries in the aggregate,
* bond levels 9-11 (HRS hierarchy) accept presence of compatible entries,
* angle 1D-5D paths are accepted only when count threshold is met,
* angle wildcard/HRS paths are used only after specific keyed paths fail.

If type-specific statistical matching fails:

* bonds fall back through HRS/element-hybrid routes, then CCP4 `ener_lib`
  compatibility lookup,
* angles fall back through HRS and then geometry defaults based on center
  hybridization/coordination (with dedicated metal and high-coordination
  handling).

This ordering keeps as much chemistry context as possible before using broad
fallbacks. It improves robustness on unusual ligands while still favoring
AceDRG-like targets whenever data are available.

Ring aromaticity and fused-ring context
---------------------------------------

Ring handling is central to output quality and compatibility. It affects
both atom typing and the final restraint lookup.

This area received substantial tuning to match AceDRG conventions:

* ring membership and aromaticity are propagated into environment labels,
* fused systems are treated as connected ring networks, not isolated
  independent rings,
* shared atoms in fused systems can keep mixed labels such as `[5a,6a]`,
* these labels are used directly in AceDRG signatures and therefore
  influence which bond/angle statistics are selected.

For aromaticity assignment, the implementation follows AceDRG's
electron-counting logic based on a Huckel-like `4n+2` criterion on
ring pi-electron totals, with AceDRG-style strict/permissive phases:

* a strict phase (used for primary statistical-table lookup),
* a permissive phase (used for output typing in edge cases),
* plus AceDRG-specific ring exceptions for selected 5-member systems.

Fallback selection also tries to preserve ring/aromatic context as long as
possible before dropping to broader generic classes. Small differences in
ring/aromatic labeling can cascade into different types and restraints,
so matching AceDRG behavior here is important for practical parity.

Special chemistry handling
--------------------------

Some chemistries need dedicated logic beyond generic rules. For example,
carborane-like systems have specialized handling aimed at reproducing
AceDRG-like typing and restraint targets more closely.

Output restraint categories
---------------------------

The generated dictionary includes standard geometric restraint families used
in crystallographic refinement:

* bond restraints (target distances + sigmas),
* angle restraints (target angles + sigmas),
* torsion restraints (including automatically inferred torsions),
* chirality restraints,
* planarity restraints.

These categories are generated from the molecular graph and assigned types,
not from a single hard-coded template per residue.

Compatibility focus and scope
-----------------------------

A major goal of this implementation is practical compatibility with
AceDRG output and conventions (not only broad chemical plausibility).
In particular, substantial effort has been invested in matching
AceDRG-like behavior for atom typing, protonation logic, and restraint
selection in edge-case chemistries.

The produced restraint values are refinement targets (ideal values and
sigmas/esds). They are empirical/statistical restraints, not a QM
geometry optimization.

For background on AceDRG algorithms, see:
`Long et al. (2017), Acta Cryst. D73, 112-122 <https://doi.org/10.1107/S2059798317000067>`_.
For the project rationale and compatibility goals, see
`Gemmi discussion #401 <https://github.com/project-gemmi/gemmi/discussions/401>`_.

.. _energy-types:

Atom typing: CCP4 energy types
------------------------------

One key concept in restraint generation is the CCP4 "energy type"
(`_chem_comp_atom.type_energy`). This is a chemistry-aware atom class
used by monomer-library restraint tables.

These types are not elements. They encode local environment features
such as hydrogen count, local bonding pattern and ring/aromatic context.
The full list of non-element types from `ener_lib.cif`
(metals and other elements use their element symbol directly):

**Carbon sp** --
``CSP`` (triple bond or two double bonds),
``CSP1`` (triple bond + 1H).

**Carbon sp2** --
``C`` (no H, carbonyl C),
``C1`` (1H),
``C2`` (2H),
``CR15`` (1H in 5-ring),
``CR16`` (1H in 6-ring),
``CR5`` (no H in 5-ring),
``CR6`` (no H in 6-ring),
``CR1`` (between two pyrrole units),
``CR1H`` (CR1 + 1H),
``CR55`` (between 5-5 rings),
``CR56`` (between 5-6 rings),
``CR66`` (between 6-6 rings).

**Carbon sp3** --
``CH1`` (1H),
``CH2`` (2H),
``CH3`` (3H),
``CT`` (no H).

**Nitrogen sp** --
``NSP``/``NS`` (no H, triple bond),
``NSP1``/``NS1`` (1H, triple bond).

**Nitrogen sp2, 3 bonds** --
``NH0`` (no H),
``NH1`` (1H, main-chain N),
``NH2`` (2H),
``NC1`` (1H, charged),
``NC2`` (2H, charged),
``NPA``/``NPB`` (no H, heme).
**Nitrogen sp2, 2 bonds** --
``N20`` (no H),
``N21`` (1H).
**Nitrogen sp2, ring** --
``NRD5`` (lone pair, 5-ring),
``NRD6`` (lone pair, 6-ring),
``NR15`` (1H, 5-ring),
``NR16`` (1H, 6-ring),
``NR5`` (3 non-H bonds, 5-ring),
``NR6`` (3 non-H bonds, 6-ring),
``NR55`` (between 5-5 rings),
``NR56`` (between 5-6 rings),
``NR66`` (between 6-6 rings).

**Nitrogen sp3, 4 bonds** --
``NT`` (no H),
``NT1`` (1H),
``NT2`` (2H),
``NT3`` (3H),
``NT4`` (4H).
**Nitrogen sp3, 3 bonds** --
``N``/``N30`` (no H),
``N31`` (1H),
``N32`` (2H),
``N33`` (3H).

**Oxygen sp2** --
``O`` (no charge, main-chain O),
``OC`` (charged),
``OP`` (charged, bonded to P),
``OS`` (charged, bonded to S),
``OB`` (charged, bonded to B).
**Oxygen sp3** --
``O2`` (2 bonds),
``OC2`` (2 bonds, charged),
``OH1`` (alcohol),
``OH2`` (water).

**Phosphorus** --
``P`` (4 bonds),
``P1`` (3 bonds).

**Sulfur** --
``S`` (no H),
``SH1`` (1H),
``S1`` (1 double bond),
``S2`` (2 bonds),
``S3`` (3 bonds, sulphoxide).

**Silicon** --
``SI`` (tetragonal),
``SI1`` (other).

**Hydrogen** --
``H`` (generic),
plus specific subtypes: ``HCH``, ``HCR``, ``HNC1``, ``HNC2``,
``HNH1``, ``HNH2``, ``HNR5``, ``HNR6``, ``HOH1``, ``HOH2``, ``HSH1``.

Historically, energy types were also used to look up ideal bond lengths
and angles. In practice this turned out to be too complex for too little gain.
Currently, their main purpose is to provide hydrogen-bonding type,
van der Waals radii, and ionic radii -- i.e. parameters for
non-bonding interactions.

AceDRG environment types (described :ref:`below <acedrg-env-signatures>`)
are far richer, but there are around 800 000 of them.
In principle they could also be used for non-bonding interaction analysis,
but this has not been done yet.

.. _acedrg-env-signatures:

AceDRG environment signatures
-----------------------------

AceDRG environment types are explicit local-neighborhood signatures.
Examples from AceDRG tables include `C(C)(H)3`, `N(CC)`, `O(C)(H)`,
`S(CC)(O)3`, `P(CC)3(O)`, and ring/aromatic-aware forms such as
`C(C[6a]C[6a]2)` and `N(C[6a]C[6a]2)`.

More complex "full" signatures can be substantially richer, for example:

* `C[5a,6a](C[5a,6a]C[6a]N[5a])(N[5a]C[5a]C[5])(N[6a]C[6a]){1|C<4>,1|N<2>,1|N<3>,1|O<2>,3|H<1>}`,
* `N[5a](C[5a,6a]C[5a,6a]N[6a])(C[5]C[5]O[5]H)(C[5a]N[5a]H){1|H<1>,1|O<2>,2|C<3>,2|C<4>}`.

In `allAtomTypesFromMolsCoded.list`, these labels are paired with coded
keys (for example `240_652_0    C(C)(H)3`).

How to read complex signatures
------------------------------

A full AceDRG signature can contain several layers of context:

* center token (for example `C[5a,6a]`) describes the central atom and
  ring/aromatic context,
* parenthesized groups encode key bonded-neighbor environments,
* optional brace blocks (for example `{1|C<4>,1|N<2>,...}`) summarize
  additional counted local features used to disambiguate similar motifs.

This expressiveness is one reason AceDRG signatures can separate subtle
chemical environments better than simpler type systems.
