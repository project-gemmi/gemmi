
gemmi/addends.hpp
    Addends to scattering form factors used in DensityCalculator
    and StructureFactorCalculator.

gemmi/align.hpp
    Sequence alignment, label_seq_id assignment, structure superposition.

gemmi/assembly.hpp
    Generating biological assemblies by applying operations
    from struct Assembly to a Model.
    Includes chain (re)naming utilities.

gemmi/asudata.hpp
    AsuData for storing reflection data.

gemmi/asumask.hpp
    AsuBrick and MaskedGrid that is used primarily as direct-space asu mask.

gemmi/atof.hpp
    Functions that convert strings to floating-point numbers ignoring locale.
    Simple wrappers around fastfloat::from_chars().

gemmi/atox.hpp
    Locale-independent functions that convert strings to integers,
    equivalents of standard isspace and isdigit, and a few helper functions.

gemmi/bessel.hpp
    Functions derived from modified Bessel functions I1(x) and I0(x).

gemmi/binner.hpp
    Binning - resolution shells for reflections.

gemmi/blob.hpp
    Finding maxima or "blobs" in a Grid (map).
    Similar to CCP4 PEAKMAX and COOT's "Unmodelled blobs".

gemmi/bond_idx.hpp
    BondIndex: for checking which atoms are bonded, calculating graph distance.

gemmi/c4322.hpp
    Electron scattering factor coefficients from the International Tables.

gemmi/calculate.hpp
    Calculate various properties of the model.

gemmi/ccp4.hpp
    CCP4 format for maps and masks.

gemmi/cellred.hpp
    Unit cell reductions: Buerger, Niggli, Selling-Delaunay.

gemmi/chemcomp.hpp
    ChemComp - chemical component that represents a monomer from Refmac
    monomer library, or from PDB CCD.

gemmi/cif.hpp
    CIF parser (based on PEGTL) with pluggable actions,
    and a set of actions that prepare Document.
    To just read the CIF format, include read_cif.hpp instead.

gemmi/cif2mtz.hpp
    A class for converting SF-mmCIF to MTZ (merged or unmerged).

gemmi/cifdoc.hpp
    struct Document that represents the CIF file (but can also be
    read from a different representation, such as CIF-JSON or mmJSON).

gemmi/contact.hpp
    Contact search, based on NeighborSearch from neighbor.hpp.

gemmi/crd.hpp
    Generate Refmac intermediate (prepared) files crd and rst

gemmi/ddl.hpp
    Using DDL1/DDL2 dictionaries to validate CIF/mmCIF files.

gemmi/dencalc.hpp
    Tools to prepare a grid with values of electron density of a model.

gemmi/dirwalk.hpp
    Classes for iterating over files in a directory tree, top-down,
    in alphabetical order. Wraps the tinydir library (as we cannot yet
    depend on C++17 <filesystem>).

gemmi/dssp.hpp
    DSSP (Define Secondary Structure of Proteins) implementation.

gemmi/ecalc.hpp
    Normalization of amplitudes F->E ("Karle" approach, similar to CCP4 ECALC).

gemmi/eig3.hpp
    Eigen decomposition code for symmetric 3x3 matrices.

gemmi/elem.hpp
    Elements from the periodic table.

gemmi/enumstr.hpp
    Converts between enums (EntityType, PolymerType, Connection::Type,
    SoftwareItem::Classification) and mmCIF strings.

gemmi/fail.hpp
    fail(), unreachable() and __declspec/__attribute__ macros

gemmi/fileutil.hpp
    File-related utilities.

gemmi/flat.hpp
    FlatStructure, FlatAtom

gemmi/floodfill.hpp
    The flood fill (scanline fill) algorithm for Grid.
    Assumes periodic boundary conditions in the grid and 6-way connectivity.

gemmi/formfact.hpp
    Calculation of atomic form factors approximated by a sum of Gaussians.
    Tables with numerical coefficients are in it92.hpp and c4322.hpp.

gemmi/fourier.hpp
    Fourier transform applied to map coefficients.

gemmi/fprime.hpp
    Cromer-Liberman calculation of anomalous scattering factors,
    with corrections from Kissel & Pratt.

gemmi/fstream.hpp
    Ofstream and Ifstream: wrappers around std::ofstream and std::ifstream.

gemmi/glob.hpp
    Glob pattern matching

gemmi/grid.hpp
    3d grids used by CCP4 maps, cell-method search and hkl data.

gemmi/gz.hpp
    Functions for transparent reading of gzipped files. Uses zlib.

gemmi/input.hpp
    Input abstraction.
    Used to decouple file reading and decompression.

gemmi/intensit.hpp
    Class Intensities that reads multi-record data from MTZ, mmCIF or XDS_ASCII
    and merges it into mean or anomalous intensities.
    It can also read merged data.

gemmi/interop.hpp
    Interoperability between Model (MX) and SmallStructure (SX).

gemmi/it92.hpp
    X-ray scattering factor coefficients from International Tables
    for Crystallography Volume C, edition from 1992 or later.

gemmi/iterator.hpp
    Bidirectional iterators (over elements of any container) that can filter,
    uniquify, group, or iterate with a stride.

gemmi/json.hpp
    Reading CIF-JSON (COMCIFS) and mmJSON (PDBj) formats into cif::Document.

gemmi/levmar.hpp
    Least-squares fitting - Levenberg-Marquardt method.

gemmi/linkhunt.hpp
    Searching for links based on the _chem_link table from monomer dictionary.

gemmi/logger.hpp
    Logger - a tiny utility for passing messages through a callback.

gemmi/math.hpp
    Math utilities. 3D linear algebra.

gemmi/metadata.hpp
    Metadata from coordinate files.

gemmi/mmcif.hpp
    Read mmCIF (PDBx/mmCIF) file into a Structure from model.hpp.

gemmi/mmcif_impl.hpp
    Functions used in both mmcif.hpp and refln.hpp (for coordinate and
    reflection mmCIF files).

gemmi/mmdb.hpp
    Converts between gemmi::Structure and mmdb::Manager.

gemmi/mmread.hpp
    Read any supported coordinate file. Usually, mmread_gz.hpp is preferred.

gemmi/mmread_gz.hpp
    Functions for reading possibly gzipped coordinate files.

gemmi/model.hpp
    Data structures to store macromolecular structure models.

gemmi/modify.hpp
    Modify various properties of the model.

gemmi/monlib.hpp
    Monomer library - (Refmac) restraints dictionary,
    which consists of monomers (chemical components), links, and modifications.

gemmi/mtz.hpp
    MTZ reflection file format.

gemmi/mtz2cif.hpp
    A class for converting MTZ (merged or unmerged) to SF-mmCIF

gemmi/neighbor.hpp
    Cell-linked lists method for atom searching (a.k.a. grid search, binning,
    bucketing, cell technique for neighbor search, etc).

gemmi/neutron92.hpp
    Neutron coherent scattering lengths of the elements,
    from Neutron News, Vol. 3, No. 3, 1992.

gemmi/numb.hpp
    Utilities for parsing CIF numbers (the CIF spec calls them 'numb').

gemmi/pdb.hpp
    Read the PDB file format and store it in Structure.

gemmi/pdb_id.hpp
    Handling PDB ID and $PDB_DIR: is_pdb_code(), expand_pdb_code_to_path(), ...

gemmi/pirfasta.hpp
    Read sequences from PIR or (multi-)FASTA formats.

gemmi/polyheur.hpp
    Heuristic methods for working with chains and polymers.
    Also includes a few well-defined functions, such as removal of waters.

gemmi/pymol_select.hpp
    PyMOL Selection syntax

gemmi/qcp.hpp
    Structural superposition, the QCP method.

gemmi/read_cif.hpp
    Functions for reading possibly gzipped CIF files.

gemmi/recgrid.hpp
    ReciprocalGrid -- grid for reciprocal space data.

gemmi/reciproc.hpp
    Reciprocal space helper functions.

gemmi/refln.hpp
    Reads reflection data from the mmCIF format.

gemmi/resinfo.hpp
    List of common residues with basic data.

gemmi/riding_h.hpp
    Place hydrogens according to bond lengths and angles from monomer library.

gemmi/scaling.hpp
    Anisotropic scaling of data (includes scaling of bulk solvent parameters).

gemmi/select.hpp
    Selections.

gemmi/seqalign.hpp
    Simple pairwise sequence alignment.

gemmi/seqid.hpp
    SeqId -- residue number and insertion code together.

gemmi/seqtools.hpp
    Functions for working with sequences (other than alignment).

gemmi/serialize.hpp
    Binary serialization for Structure (as well as Model, UnitCell, etc).

gemmi/sfcalc.hpp
    Direct calculation of structure factors.

gemmi/small.hpp
    Representation of a small molecule or inorganic crystal.
    Flat list of atom sites. Minimal functionality.

gemmi/smcif.hpp
    Read small molecule CIF file into SmallStructure (from small.hpp).

gemmi/solmask.hpp
    Flat bulk solvent mask. With helper tools that modify data on grid.

gemmi/span.hpp
    Span - span of array or std::vector.
    MutableVectorSpan - span of std::vector with insert() and erase()

gemmi/sprintf.hpp
    interface to stb_sprintf: snprintf_z, to_str(float|double)

gemmi/stats.hpp
    Statistics utilities: classes Covariance, Correlation, DataStats

gemmi/symmetry.hpp
    Crystallographic Symmetry. Space Groups. Coordinate Triplets.

gemmi/to_chemcomp.hpp
    Create cif::Block with monomer library _chem_comp* categories
    from struct ChemComp.

gemmi/to_cif.hpp
    Writing cif::Document or its parts to std::ostream.

gemmi/to_json.hpp
    Writing cif::Document or its parts as JSON (mmJSON, CIF-JSON, etc).

gemmi/to_mmcif.hpp
    Create cif::Document (for PDBx/mmCIF file) from Structure.

gemmi/to_pdb.hpp
    Writing PDB file format (Structure -> pdb file).

gemmi/topo.hpp
    Topo(logy) - restraints (from a monomer library) applied to a model.

gemmi/twin.hpp
    Twinning laws.

gemmi/unitcell.hpp
    Unit cell.

gemmi/utf.hpp
    Conversion between UTF-8 and wchar. Used only for file names on Windows.

gemmi/util.hpp
    Utilities. Mostly for working with strings and vectors.

gemmi/version.hpp
    Version number.

gemmi/xds2mtz.hpp
    Convert XDS_ASCII to MTZ.

gemmi/xds_ascii.hpp
    Read XDS files: XDS_ASCII.HKL and INTEGRATE.HKL.
