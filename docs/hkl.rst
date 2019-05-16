
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
The headers, as well as data, are stored in class ``Mtz``.
The columns of data are grouped hierarchically (Project -> Crystal -> Dataset
-> Column), but normally the column name is all that is needed,
and the hierarchy can be ignored. In Gemmi, the hierarchy is flattened:
we have a list of columns and a list of datasets.
Each columns is associated with one dataset, and each dataset have properties
``project_name`` and ``crystal_name`` (in addition to ``dataset_name``),
which is enough to construct tree-like hierarchy if needed.


Reading
-------

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

Modifying
---------

To show how the data and metadata in an ``Mtz`` object can be modified,
we will create a complete ``Mtz`` object from scratch.
The code below is in Python, but all the functions and properties have
equivalents with the same names in C++.

We start from creating an object and setting its space group and unit cell:

.. doctest::

  >>> mtz = gemmi.Mtz()
  >>> mtz.spacegroup = gemmi.find_spacegroup_by_name('P 21 21 2')
  >>> mtz.cell.set(77.7, 149.5, 62.4, 90, 90, 90)

Now we need to add datasets and columns.

.. doctest::

  >>> mtz.add_dataset('HKL_base')
  <gemmi.Mtz.Dataset 0 HKL_base/HKL_base/HKL_base>

The name passed to ``add_dataset()`` is used to set ``dataset_name``,
``crystal_name`` and ``project_name``. These three names are often
kept the same, but if needed, they can be changed afterwards.

The MTZ format stores general, per file cell dimensions (keyword CELL),
as well per dataset ones (keyword DCELL).
Function ``add_dataset`` initializes the latter from the former.
This also can be changed afterwards.

The first three columns must always be Miller indices:

.. doctest::

  >>> for label in ['H', 'K', 'L']: mtz.add_column(label, 'H')
  <gemmi.Mtz.Column H type H>
  <gemmi.Mtz.Column K type H>
  <gemmi.Mtz.Column L type H>

Let's add two more columns in a separate dataset:

.. doctest::

  >>> mtz.add_dataset('synthetic')
  <gemmi.Mtz.Dataset 1 synthetic/synthetic/synthetic>
  >>> mtz.add_column('F', 'F')
  <gemmi.Mtz.Column F type F>
  >>> mtz.add_column('SIGF', 'Q')
  <gemmi.Mtz.Column SIGF type Q>

By default, the new column is assigned to the last dataset
and is placed at the end of the column list.
But we can choose any dataset and position:

.. doctest::

  >>> mtz.add_column('FREE', 'I', dataset_id=0, pos=3)
  <gemmi.Mtz.Column FREE type I>
  >>> mtz.column_labels()
  ['H', 'K', 'L', 'FREE', 'F', 'SIGF']

Now it is time to add data.
In Python we have function ``set_data()`` that expects 2D NumPy array
of floating point numbers (even indices are converted to floats,
but they need to be converted at some point anyway -- the MTZ format
stores all numbers as 32-bit floats).

.. doctest::

  >>> data = numpy.array([[4, 13, 8, 1, 453.9, 19.12],
  ...                     [4, 13, 9, 0, 102.0, 27.31]], numpy.float)
  >>> mtz.set_data(data)
  >>> mtz
  <gemmi.Mtz with 6 columns, 2 reflections>

In C++ the ``set_data`` function takes a pointer to row-wise ordered data
and its size (columns x rows)::

  void Mtz::set_data(const float* new_data, size_t n)

Writing
-------

In C++, the MTZ file can be written to a file using one of the functions::

  void Mtz::write_to_stream(std::FILE* stream) const
  void Mtz::write_to_file(const std::string& path) const

and in Python using:

.. doctest::
  :skipif: no_mtz_file

  >>> mtz.write_to_file('output.mtz')


SF mmCIF
========

The mmCIF format is also used to store reflection data from crystallographic
experiments. Reflections and coordinates are normally stored in separate
mmCIF files. The files with reflections are called structure factor mmCIF
or shortly SF mmCIF. Such file for the 1ABC PDB entry
is usually named either ``1ABC-sf.cif`` (if downloaded through RCSB website)
or ``r1abcsf.ent`` (if downloaded through PDBe website or through FTP).

SF mmCIF files usually contain one block, but may have
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
