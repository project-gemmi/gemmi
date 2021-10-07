
gemmi/addends.hpp
    Addends to scattering form factors used in DensityCalculator
    and in StructureFactorCalculator.

gemmi/align.hpp
    Sequence alignment, label_seq_id assignment, structure superposition.

gemmi/assembly.hpp
    Generating biological assemblies by applying operations
    from struct Assembly to a Model.
    Includes chain (re)naming utilities.

gemmi/asudata.hpp
    AsuData for storing reflection data.

gemmi/atof.hpp
    Functions that convert string to floating-point number ignoring locale.
    Simple wrappers around fastfloat::from_chars().

gemmi/atox.hpp
    Locale-independent functions that convert string to integer,
    equivalents of standard isspace and isdigit, and a few helper functions.

gemmi/bessel.hpp
    Ratio of modified Bessel functions I1(x)/I0(x).

gemmi/blob.hpp
    Finding maxima or "blobs" in a Grid (map).
    Similar to CCP4 PEAKMAX and COOT's "Unmodelled blobs".

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

gemmi/chemcomp_xyz.hpp
    Reading coordinates from chemical component or Refmac monomer library files.

gemmi/cif.hpp
    CIF parser (based on PEGTL) with pluggable actions,
    and a set of actions that prepare Document.

gemmi/cif2mtz.hpp
    A class for converting SF-mmCIF to MTZ (merged or unmerged).

gemmi/cifdoc.hpp
    struct Document that represents the CIF file (but can be also
    read from JSON file, such as CIF-JSON or mmJSON).

gemmi/contact.hpp
    Contact search, based on NeighborSearch from neighbor.hpp.

gemmi/dencalc.hpp
    Tools to prepare a grid with values of electron density of a model.

gemmi/dirwalk.hpp
    Classes for iterating files in a directory tree, top-down,
    in an alphabetical order.  It wraps the tinydir library (as we cannot
    depend on C++17 <filesystem> yet).

gemmi/elem.hpp
    Elements from the periodic table.

gemmi/enumstr.hpp
    Converts between enums (EntityType, PolymerType, Connection::Type,
    SoftwareItem::Classification) and mmCIF strings.

gemmi/fail.hpp
    fail() and unreachable()

gemmi/fileutil.hpp
    File-related utilities.

gemmi/floodfill.hpp
    The flood fill (scanline fill) algorithm for Grid.
    Assumes periodic boundary conditions in the grid and 6-way connectivity.

gemmi/formfact.hpp
    Calculation of atomic form factors approximated by a sum of Gaussians.
    Tables with numeric coefficient are in it92.hpp and c4322.hpp.

gemmi/fourier.hpp
    Fourier transform applied to map coefficients.

gemmi/fprime.hpp
    C++ implementation of Cromer-Liberman calculation of anomalous scattering
    factors, with corrections from Kissel & Pratt, Acta Cryst. A46, 170 (1990).
    Single header. No dependencies.

gemmi/fstream.hpp
    Ofstream and Ifstream: wrappers around std::ofstream and std::ifstream.

gemmi/grid.hpp
    3d grids used by CCP4 maps, cell-method search and hkl data.

gemmi/gz.hpp
    Functions for transparent reading of gzipped files. Uses zlib.

gemmi/input.hpp
    Input abstraction.
    Used to decouple file reading and uncompression.

gemmi/interop.hpp
    Interoperability between Model (MX) and SmallStructure (SX).

gemmi/it92.hpp
    X-ray scattering factor coefficients from International Tables
    for Crystallography Volume C, edition from 1992 or later.

gemmi/iterator.hpp
    Iterators. Currently each of them is a BidirectionalIterator.

gemmi/json.hpp
    Reading CIF-JSON (COMCIFS) and mmJSON (PDBj) formats into cif::Document.

gemmi/levmar.hpp
    Least-squares fitting - Levenberg-Marquardt method.

gemmi/linkhunt.hpp
    Searching for links based on the _chem_link table from monomer dictionary.

gemmi/math.hpp
    Math utilities. 3D linear algebra.

gemmi/merge.hpp
    Class Intensities that reads multi-record data from MTZ, mmCIF or XDS_ASCII
    and merges it into mean or anomalous intensities.
    It can also read merged data.

gemmi/metadata.hpp
    Metadata from coordinate files.

gemmi/mmcif.hpp
    Read mmcif (PDBx/mmCIF) file into a Structure from model.hpp.

gemmi/mmcif_impl.hpp
    Function used in both mmcif.hpp and refln.hpp (for coordinate and
    reflection mmCIF files).

gemmi/mmread.hpp
    Read any supported coordinate file.

gemmi/model.hpp
    Data structures to keep macromolecular structure model.

gemmi/monlib.hpp
    Monomer library - (Refmac) restraints dictionary,
    which is made of monomers (chemical components), links and modifications.

gemmi/mtz.hpp
    MTZ reflection file format.

gemmi/mtz2cif.hpp
    A class for converting MTZ (merged or unmerged) to SF-mmCIF

gemmi/neighbor.hpp
    Cell-linked lists method for atom searching (a.k.a. grid search, binning,
    bucketing, cell technique for neighbor search, etc).

gemmi/numb.hpp
    Utilities for parsing CIF numbers (the CIF spec calls it 'numb').

gemmi/pdb.hpp
    Read PDB file format and store it in Structure.

gemmi/pirfasta.hpp
    Read sequence from PIR or FASTA format.

gemmi/placeh.hpp
    Place hydrogens according to bond lengths and angles from monomer library.

gemmi/polyheur.hpp
    Heuristic methods for working with chains and polymers.
    Includes also a few well-defined functions, such as removal of hydrogens.

gemmi/qcp.hpp
    Structural superposition, the QCP method.

gemmi/read_cif.hpp
    Functions for reading possibly gzipped CIF files.
    Trivial wrappers that can make compilation faster.

gemmi/read_coor.hpp
    Functions for reading possibly gzipped coordinate files.
    Trivial wrappers that can make compilation faster.

gemmi/read_map.hpp
    Functions for reading possibly gzipped CCP4 map files.
    Trivial wrappers that can make compilation faster.

gemmi/recgrid.hpp
    ReciprocalGrid -- grid for reciprocal space data.

gemmi/reciproc.hpp
    Reciprocal space utilities.

gemmi/refln.hpp
    Reads reflection data from the mmCIF format.

gemmi/reindex.hpp
    Reindex merged or unmerged MTZ file.

gemmi/remarks.hpp
    Function read_metadata_from_remarks() that interprets REMARK 3
    and REMARK 200/230/240 filling in Metadata.

gemmi/resinfo.hpp
    List of common residues with basic data.

gemmi/scaling.hpp
    Anisotropic scaling of data (includes scaling of bulk solvent parameters)

gemmi/select.hpp
    Selections.

gemmi/seqalign.hpp
    Simple pairwise sequence alignment.

gemmi/seqid.hpp
    SeqId -- residue number and insertion code together.

gemmi/sfcalc.hpp
    Direct calculation of structure factors.

gemmi/small.hpp
    Representation of small molecule or inorganic crystal.
    Flat list of atom sites. Minimal functionality.

gemmi/smcif.hpp
    Read small molecule CIF file into SmallStructure (from small.hpp).

gemmi/solmask.hpp
    Flat bulk solvent mask. With helper tools that modify data on grid.

gemmi/span.hpp
    Span - span of array or std::vector.
    MutableVectorSpan - span of std::vector with insert() and erase()

gemmi/sprintf.hpp
    to_str(float|double), gf_snprintf - wrappers around stb_sprintf.

gemmi/symmetry.hpp
    Crystallographic Symmetry. Space Groups. Coordinate Triplets.

gemmi/to_cif.hpp
    Writing cif::Document or its parts to std::ostream.

gemmi/to_json.hpp
    Writing cif::Document or its parts as JSON (mmJSON, CIF-JSON, etc).

gemmi/to_mmcif.hpp
    Create cif::Document (for PDBx/mmCIF file) from Structure.

gemmi/to_mmdb.hpp
    copy_to_mmdb(): converts gemmi::Structure to mmdb::Manager.

gemmi/to_pdb.hpp
    Writing PDB file format (Structure -> pdb file).

gemmi/topo.hpp
    Topo(logy) - restraints (from a monomer library) applied to a model.

gemmi/tostr.hpp
    gemmi::tostr() - converts a list of arguments to string (uses ostringstream).

gemmi/unitcell.hpp
    Unit cell.

gemmi/utf.hpp
    Conversion between UTF-8 and wchar. Used only for file names on Windows.

gemmi/util.hpp
    Utilities. Mostly for working with strings and vectors.

gemmi/version.hpp
    Version number.

gemmi/xds_ascii.hpp
    Read XDS_ASCII.HKL. For now, only unmerged files are read.
