
Reflection files
################

This section is about reciprocal space and is primarily for people
working with crystallographic data.

Reflection files list reflections, identified by Miller indices (h,k,l),
with various numbers assigned to them.
Initially, the numbers represent observations derived from
diffraction images (estimated intensities and their errors).
Then other quantities derived from these are added,
as well as quantities derived from a macromolecular model
or from both the model and experimental data.

Electron density maps, which are used in real space, are also kept as
numbers assigned to reflections. They must be Fourier-transformed
to the real space before use, but this operation is cheap enough
to be performed on the fly.

Gemmi can read two formats of reflection files:

* MTZ files -- the most popular format in macromolecular crystallography,
* and mmCIF files -- used for data archiving in the Protein Data Bank.

MTZ
===

MTZ format has textual headers and a binary data table, where all numbers
are stored in a 32-bit floating point format.
The headers, as well as data, are read into ``Mtz`` class.

In C++, the MTZ file can be read using either stand-alone functions::

  Mtz read_mtz_file(const std::string& path)
  Mtz read_mtz_stream(std::FILE* stream, bool with_data)

or member functions of the Mtz class, when more control over the reading
process is needed.

In Python, we have a single function for reading MTZ files:

.. doctest::
  :skipif: no_mtz_file

  >>> import gemmi
  >>> mtz = gemmi.read_mtz_file('example.mtz')

The Mtz class has a number of properties that are based on the header
records in the MTZ format, such as ``spacegroup``, ``cell``,
``title``, ``history``. (TODO: document them all.)
Importantly, it also has a list of datasets and a list of columns.

Datasets are stored in variable ``datasets``::

  std::vector<Mtz::Dataset> Mtz::datasets

.. doctest::
  :skipif: no_mtz_file

  >>> mtz.datasets
  MtzDatasets[<gemmi.Mtz.Dataset 0 HKL_base/HKL_base/HKL_base>, <gemmi.Mtz.Dataset 1 5dei/5dei/1>]

In the MTZ file, each dataset is identified internally by an integer
"dataset ID". To get dataset with specified ID use function::

  Dataset& Mtz::dataset(int id)

.. doctest::
  :skipif: no_mtz_file

  >>> mtz.dataset(0)
  <gemmi.Mtz.Dataset 0 HKL_base/HKL_base/HKL_base>

Dataset has a few properties that can be accessed directly::

  struct Dataset {
    int number;
    std::string project_name;
    std::string crystal_name;
    std::string dataset_name;
    UnitCell cell;
    double wavelength = 0.;
  };

Python bindings provide the same properties:

.. doctest::
  :skipif: no_mtz_file

  >>> mtz.dataset(0).project_name
  'HKL_base'
  >>> mtz.dataset(0).crystal_name
  'HKL_base'
  >>> mtz.dataset(0).dataset_name
  'HKL_base'
  >>> mtz.dataset(0).cell
  <gemmi.UnitCell(70.166, 92.451, 93.715, 63.586, 72.425, 72.884)>
  >>> mtz.dataset(0).wavelength
  0.0


Columns are stored in variable ``columns``::

  std::vector<Mtz::Column> Mtz::columns

.. doctest::
  :skipif: no_mtz_file

  >>> mtz.columns[0]
  <gemmi.Mtz.Column H type H>
  >>> len(mtz.columns)
  6

To get the first column with the specified label or type use functions:

.. doctest::
  :skipif: no_mtz_file

  >>> mtz.column_with_label('FREE')
  <gemmi.Mtz.Column FREE type I>
  >>> mtz.column_with_type('Q')
  <gemmi.Mtz.Column SIGI type Q>

Column has properties read from MTZ headers COLUMN and COLSRC::

  struct Column {
    int dataset_id;
    char type;
    std::string label;
    float min_value = NAN;
    float max_value = NAN;
    std::string source;  // from COLSRC
    // (and functions and variables that help in accessing data)
  };

Python bindings provide the same properties:

.. doctest::
  :skipif: no_mtz_file

  >>> intensity = mtz.column_with_label('I')
  >>> intensity.dataset_id
  1
  >>> intensity.type
  'J'
  >>> intensity.label
  'I'
  >>> intensity.min_value, intensity.max_value
  (-513.5999755859375, 138805.0)
  >>> intensity.source
  'CREATED_31/01/2019_20:31:27'

In both C++ and Python ``Column`` supports the iteration protocol:

.. doctest::
  :skipif: no_mtz_file

  >>> sum(intensity) / len(intensity)
  2665.906147073007

In Python we also have a read-only ``array`` property that provides
a view of the data compatible with NumPy. It does not copy the data,
so the data is not contiguous (because it's stored row-wise in MTZ):

.. doctest::
  :skipif: no_mtz_file

  >>> intensity.array.strides
  (24,)

Here is an example that uses the array property
to make a plot similar to `AUSPEX <http://www.auspex.de/>`_:

.. literalinclude:: ../examples/mtz_i_sigi.py
  :language: python
  :lines: 4-

To test this program we run it on data from 5DEI and we get result
similar to `#3 <http://www.auspex.de/pathol/#3>`_:

.. image:: img/5dei-Ioversigma.png
  :align: center
  :scale: 100

Another way to access all the MTZ data in Python is through the
`buffer protocol <https://docs.python.org/3/c-api/buffer.html>`_:
For example, to view the data as 2D NumPy array (C-style contiguous)
without copying do:

.. doctest::
  :skipif: no_mtz_file

  >>> import numpy
  >>> all_data = numpy.array(mtz, copy=False)

It helps to have labels on the columns. A good data structure for this
is Pandas DataFrame:

.. doctest::
  :skipif: no_mtz_file

  >>> import pandas
  >>> df = pandas.DataFrame(data=all_data, columns=mtz.column_labels())
  >>> # now we can handle columns using their labels:
  >>> I_over_sigma = df['I'] / df['SIGI']

In C++, the MTZ file can be written to a file using one of the functions::

  void Mtz::write_to_stream(std::FILE* stream) const
  void Mtz::write_to_file(const std::string& path) const

and in Python using:

.. doctest::
  :skipif: no_mtz_file

  >>> mtz.write_to_file('output.mtz')


SF mmCIF
========

Reflection data mmCIF files usually contain one block, but may have
multiple blocks, for example merged and unmerged data in separate blocks.
Usually blocks contain information about the unit cell and space group,
sometimes also about the radiation wavelength.
Usually, the reflections are listed as mmCIF category ``_refln``,
but in some blocks the ``_diffrn_refln`` category is used
to provide intensities.

The script below renders the same colorful *I*/*σ* image as in the previous
section, but it can take as an argument a file downloaded directly from
the wwPDB (for example,
``$PDB_DIR/structures/divided/structure_factors/de/r5deisf.ent.gz``).

.. literalinclude:: ../examples/cif_i_sigi.py
  :language: python
  :lines: 4-
