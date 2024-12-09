.. _scattering:

Scattering
==========

On this page:

* First, we cover parameters describing how atoms scatter
  :ref:`X-rays <xray_scattering>`,
  :ref:`electrons <electron_scattering>` and
  :ref:`neutrons <neutron_scattering>`.
  This data can be viewed as properties of :ref:`Elements <elements>`.
* Next, we introduce functions that use these properties to calculate
  structure factors (:ref:`Direct summation <direct_summation>`),
  and functions that use the same properties to calculate
  :ref:`real-space density <density>`,
  which can then be used to calculate structure factors more efficiently.
* Finally, we cover complementary calculations, such as approximating
  the contribution of bulk solvent to the total scattering.

.. _xray_scattering:

X-ray form factors
------------------

The X-ray scattering contribution from an isolated atom is described by the
atomic form factor. Nowadays, two form factor parametrizations are commonly
used: one from the International Tables of Crystallography Vol. C
(4 Gaussians + constant factor), and the other from Waasmaier & Kirfel (1995),
`Acta Cryst. A51, 416 <https://doi.org/10.1107/S0108767394013292>`_
(5 Gaussians).
Other parametrizations exist (for example, with only two Gaussians),
but are not widely used.

Currently, Gemmi includes only the ITC parametrization.

In C++, the form factor coefficients are listed in the :file:`it92.hpp` header.
In Python, these coefficients can be accessed as a property of an element
(only coefficients for neutral atoms):

.. doctest::

  >>> gemmi.Element('Fe').it92  #doctest: +ELLIPSIS
  <gemmi.IT92Coef object at 0x...>
  >>> gemmi.Element('Es').it92  # -> None (no parametrization for Einsteinium)

or by using the function `IT92_get_exact()` that takes an element and a charge
as arguments and returns None if this exact atom or ion is absent in the table.

The charges listed in a coordinates file might not be reliable for determining
the scattering of an atom, so by default, we use coefficients for neutral atoms
regardless of the formal charge. To use coefficients for ions, you first need to
change `IT92::ignore_charge`:

.. doctest::

  >>> gemmi.IT92_set_ignore_charge(False)
  >>> gemmi.IT92_get_exact(gemmi.Element('Mg'), +2)  # for Mg2+ #doctest: +ELLIPSIS
  <gemmi.IT92Coef object at 0x...>

You can get the coefficients as numbers:

.. doctest::

  >>> fe_coef = gemmi.Element('Fe').it92
  >>> fe_coef.a
  [11.7695, 7.3573, 3.5222, 2.3045]
  >>> fe_coef.b
  [4.7611, 0.3072, 15.3535, 76.8805]
  >>> fe_coef.c
  1.0369

The *a*'s and *c* should sum up to the number of electrons (*Z* - *charge*),
but the numbers from the Tables may differ slightly:

.. doctest::

  >>> gemmi.Element('Fe').atomic_number
  26
  >>> sum(fe_coef.a) + fe_coef.c  # doctest: +ELLIPSIS
  25.9904...

Some programs (Refmac) normalize the coefficients,
which means multiplying them by a factor between 0.99995 and 1.00113.

.. doctest::

  >>> # let's store the original values first
  >>> orig_coefs = {i : gemmi.Element(i).it92.get_coefs() for i in range(1, 99)}
  >>> # now we normalize the values
  >>> gemmi.IT92_normalize()
  >>> # and we can see that the values have changed
  >>> fe_coef.a  # doctest: +ELLIPSIS
  [11.7738..., 7.3600..., 3.5235..., 2.30535...]

Now let's use the function `set_coefs()` to change the coefficients back
to the original values (for neutral atoms):

.. doctest::

  >>> for i in range(1, 99):
  ...     gemmi.Element(i).it92.set_coefs(orig_coefs[i])


Macromolecular models may have unknown atoms with the element specified as X.
By default, we use oxygen's coefficients for X,
but you may change it (as well as the coefficients of any other atom):

.. doctest::

  >>> c_coefs = gemmi.Element('C').it92.get_coefs()
  >>> gemmi.Element('X').it92.set_coefs(c_coefs)

The coefficients can be used to directly calculate the sum of Gaussians --
the structure factor contribution:

.. doctest::
  :skipif: numpy is None

  >>> fe_coef.calculate_sf(stol2=0.4)  # argument: (sin(theta)/lambda)^2
  9.303603172302246

The large number of reflections in macromolecular crystallography makes direct
calculation of structure factors inefficient. Instead, we can calculate electron
density on a grid and use FFT to obtain structure factors.
In the simplest case, the atom's contribution to the electron density at a grid
point can be calculated as:

.. doctest::
  :skipif: numpy is None

  >>> # arguments are distance^2 and isotropic ADP
  >>> fe_coef.calculate_density_iso(r2=2.3, B=50)
  0.5279340744018555

The C++ interface provides more functions to calculate the electron density.
We have separate functions to work with isotropic and anisotropic ADPs.
For efficiency, the computations are split into two stages: the first one is
performed once per atom, and the second one for each nearby grid point.

In the usual scenario, we add *f'* (the real component of anomalous
scattering -- see the next section) to the constant coefficient *c*.

.. _anomalous:

Anomalous scattering
--------------------

The anomalous dispersion is wavelength-dependent.
Gemmi provides the function `cromer_liberman`, which calculates
the real and imaginary components *f'* and *f"*
for isolated atoms from Z=3 to Z=92.

As the name suggests, this function uses the Cromer-Liberman algorithm.
This algorithm, as noted on the
`pyFprime website <https://subversion.xray.aps.anl.gov/trac/pyFprime/>`_,
"fails in computing *f'* for wavelengths < 0.16 Å (> 77.48 keV)
for the heaviest elements (Au-Cf)
and fails to correctly compute *f'*, *f"* and *μ*
for wavelengths > 2.67 Å (< 4.64 keV) for very heavy elements (Am-Cf)."

The implementation is contained in a single C++ header file
:file:`fprime.hpp` with no dependencies.
All the data is embedded in the code.
The binary size after compilation is about 100kB.

Reportedly, the data tables synthesized by C.T. Chantler are more accurate.
Consider using them instead. They are available from the
`NIST website <https://www.nist.gov/pml/x-ray-form-factor-attenuation-and-scattering-tables>`_
and the `XrayDB <https://github.com/xraypy/XrayDB>`_ project.

The code in :file:`fprime.hpp` is based on the X-ray spectroscopy project
`Larch <https://xraypy.github.io/xraylarch/>`_
and should give the same results as the
`f1f2_cl <https://xraypy.github.io/xraylarch/xray/index.html#_xray.f1f2_cl>`_
function there. The Fortran code in Larch is, in turn, based on the
`Brennan and Cowan <https://aip.scitation.org/doi/10.1063/1.1142625>`_
routines. Which, in turn, were based on the
`original program <https://doi.org/10.1107/S0021889883010791>`_
from Don Cromer. Along the way, the code was extensively modified.
Importantly, the Jensen correction has been removed (as is recommended
in the chapter 4.2.6 of `ITvC <https://it.iucr.org/Cb/contents/>`_)
and the `Kissel and Pratt (1990) <https://doi.org/10.1107/S0108767389010718>`_
correction has been added.
Therefore, it gives different results than the
`crossec <http://legacy.ccp4.ac.uk/html/crossec.html>`_ program,
which was contributed to CCP4 directly by Don Cromer in the 1990's.

The `cromer_liberman` function is available in both C++ and Python:

.. doctest::

  >>> gemmi.Element('Se').atomic_number
  34
  >>> gemmi.cromer_liberman(z=_, energy=10332.0) # energy in eV  #doctest: +ELLIPSIS
  (-1.41862..., 0.72389...)
  >>> # use gemmi.hc to convert wavelength [A] to energy [eV]
  >>> gemmi.cromer_liberman(z=34, energy=gemmi.hc/0.71073)  #doctest: +ELLIPSIS
  (-0.09201..., 2.23336...)

The same values can be printed from the command line program
:ref:`gemmi-fprime <fprime>`.

.. _electron_scattering:

Electron form factors
---------------------

The electron form factors are parametrized as 5 Gaussians,
using coefficients from International Tables for Crystallography
Volume C, edition 2011, table 4.3.2.2 (pp. 282-283),
"Elastic atomic scattering factors of electrons for neutral atoms
and *s* up to 2.0 A\ :sup:`--1`".
The same parametrization is used in cctbx and CCP4.

In gemmi, we refer to these coefficients as C4322 (after the volume and
table number in the International Tables).
For C++, the form factor coefficients are tabulated
in the :file:`c4322.hpp` header.
In Python, they can be accessed as a property of an element:

.. doctest::

  >>> fe_coef = gemmi.Element('Fe').c4322
  >>> fe_coef.a
  [0.3946, 1.2725, 1.7031, 2.314, 1.4795]
  >>> fe_coef.b
  [0.2717, 2.0443, 7.6007, 29.9714, 86.2265]

Similarly to the X-ray coefficients, they can be used to calculate
structure factors and density:

.. doctest::

  >>> fe_coef.calculate_sf(stol2=0.4)  # argument: (sin(theta)/lambda)^2
  0.9971507787704468
  >>> fe_coef.calculate_density_iso(r2=2.3, B=50)
  0.13794779777526855

Unlike for X-ray form factors, we do not add anomalous scattering here.

.. _neutron_scattering:

Neutron scattering
------------------

For neutrons, we use bound coherent scattering lengths
`from NCNR <https://www.ncnr.nist.gov/resources/n-lengths/list.html>`_,
which are based on Neutron News, Vol. 3, No. 3, 1992.

In C++, this data is in the :file:`neutron92.hpp` header.
In Python, it can be accessed as a property of an element:

.. doctest::

  >>> h_neut = gemmi.Element('H').neutron92
  >>> h_neut.get_coefs()  # bound coherent scat. length in fm (10^-15 m)
  [-3.739]
  >>> h_neut.calculate_sf(0)  # the argument is ignored
  -3.739
  >>> h_neut.calculate_density_iso(r2=1.3, B=25)
  -0.17104358254308388

We don't store data for individual isotopes, except for deuterium.

.. doctest::

  >>> gemmi.Element('D').neutron92.get_coefs()
  [6.671]

.. _direct_summation:

Direct summation
----------------

Gemmi has functions to calculate structure factor by direct summation
of structure factors from individual atoms. This route is not commonly used
in macromolecular crystallography and it was implemented primarily to check
the accuracy of FFT-based computations. Therefore, the efficiency was not
a priority here. In particular, space group specific optimizations, described
by `Bourhis et al (2014) <https://doi.org/10.1107/S2053273314022207>`_,
are not included.

In Python classes `StructureFactorCalculatorX` and `StructureFactorCalculatorE`
perform direct summation using X-ray and electron form factors, respectively.
Class `StructureFactorCalculatorN` does calculations for neutrons.
The C++ interface is similar, although it uses a single templated class
`StructureFactorCalculator`.
These classes are constructed with `UnitCell` as a parameter (we don't pass
`SpaceGroup` because `UnitCell` already contains a list of symmetry operations).

.. doctest::

  >>> st = gemmi.read_structure('../tests/4oz7.pdb')
  >>> calc_e = gemmi.StructureFactorCalculatorE(st.cell)

Now we can compute structure factors from Model for any (hkl):

.. doctest::

  >>> calc_e.calculate_sf_from_model(st[0], (3,4,5))  #doctest: +ELLIPSIS
  (54.50873...+53.39498...j)

Similarly, for small molecules:

.. doctest::
  :skipif: sys.platform == 'win32'

  >>> small = gemmi.read_small_structure('../tests/1011031.cif')
  >>> small.change_occupancies_to_crystallographic()
  >>> calc_x = gemmi.StructureFactorCalculatorX(small.cell)
  >>> calc_x.calculate_sf_from_small_structure(small, (0,2,4))
  (17.849694728851315-6.557871454633539e-15j)

For each atom, the Debye-Waller factor (used in the structure factor
calculation) is obtained using either isotropic or anisotropic ADPs
(B-factors). If anisotropic ADPs are non-zero, isotropic ADP is ignored.

.. _addends:

Addends
~~~~~~~

`StructureFactorCalculator` contains also addends -- per-element values
that are to be added to the form factor coefficient *c*
(which is the same as adding it to the total atomic contribution
to the structure factor):

.. doctest::

  >>> calc_x = gemmi.StructureFactorCalculatorX(st.cell)
  >>> calc_x.addends  #doctest: +ELLIPSIS
  <gemmi.Addends object at 0x...>

When calculating X-ray structure factors, one may want to include *f'*
(real part of the :ref:`anomalous scattering <anomalous>`).
To do this, *f'* needs to be set for each element present in the system:
This is done with function `set_addend()`, which sets angle-independent
value that will be added to the value calculated from the form factors.

.. doctest::

  >>> energy = gemmi.hc / 0.8  # for wavelength 0.8A
  >>> for symbol in ['C', 'N', 'O', 'S', 'Cu']:
  ...     el = gemmi.Element(symbol)
  ...     fp, _ = gemmi.cromer_liberman(z=el.atomic_number, energy=energy)
  ...     calc_x.addends.set(el, fp)
  ...
  >>> calc_x.addends.get(gemmi.Element('Cu'))
  0.265503853559494

Alternatively, you could call a single function (here we first need
to remove the previous values, so it makes two functions):

.. doctest::

  >>> calc_x.addends.clear()
  >>> calc_x.addends.add_cl_fprime(energy)

which calculates *f'* for all elements handled by the Cromer-Liberman
algorithm (*Z* from 3 to 92). Although it seems wasteful, it takes
well below 1ms.

Structure factors calculated at this point incorporate the addends:

.. doctest::

  >>> calc_x.calculate_sf_from_model(st[0], (3,4,5))  #doctest: +ELLIPSIS
  (182.3655...+269.0002...j)

Addends can also be employed to calculate the electron scattering
from X-ray form factors, according to the Mott–Bethe formula:

.. doctest::

  >>> calc_x.addends.clear()
  >>> calc_x.addends.subtract_z()
  >>> calc_x.mott_bethe_factor() * calc_x.calculate_sf_from_model(st[0], (3,4,5))  #doctest: +ELLIPSIS
  (54.065...+52.968...j)

The next section gives slightly more details on the Mott-Bethe formula.

.. _density:

Density for FFT
---------------

To use FFT to calculate structure factors, we first need to calculate
density of the scatterer (usually electrons) on a grid. For this we use

* in C++ -- class DensityCalculator templated with a form factor table,
* in Python -- classes DensityCalculatorX (corresponding to X-ray
  form factors), DensityCalculatorE (electron form factors)
  and DensityCalculatorN (neutrons).

DensityCalculator contains a grid. The size of the grid is determined
from two parameters that we need to set: `d_min` which corresponds to
our resolution limit, and `rate` -- the oversampling rate (1.5 by default).

.. doctest::

  >>> dencalc = gemmi.DensityCalculatorE()
  >>> dencalc.d_min = 2.5 # 2.5A
  >>> dencalc.rate = 1.5  # we could skip it, this is the default value

These two parameters are only used to determine the spacing of the grid.
In this case, about 2.5Å / (2 \* 1.5) = 0.83Å.

As with StructureFactorCalculator, here we also have :ref:`addends <addends>`:

.. doctest::

  >>> dencalc.addends  #doctest: +ELLIPSIS
  <gemmi.Addends object at 0x...>

To create a grid and calculate the density we use two function calls.
Almost all the work is in the latter:

.. doctest::

  >>> dencalc.grid.setup_from(st)
  >>> dencalc.put_model_density_on_grid(st[0])

Calling `put_model_density_on_grid` is equivalent to these three functions:

.. doctest::

  >>> dencalc.initialize_grid()
  >>> dencalc.add_model_density_to_grid(st[0])
  >>> dencalc.grid.symmetrize_sum()

`initialize_grid()`, in this case, uses `d_min` and `rate`
to determine required grid spacing and uses this spacing to setup the grid.
If `d_min` is not set, but the grid size is already set,
initialize_grid() only zeros all the grid values.

At this point, we have a grid with density:

.. doctest::

  >>> dencalc.grid
  <gemmi.FloatGrid(48, 48, 50)>

which can be transformed (as described :ref:`above <fft>`)
into a structure factor grid:

.. doctest::

  >>> sf_grid = gemmi.transform_map_to_f_phi(dencalc.grid)
  >>> sf_grid
  <gemmi.ReciprocalComplexGrid(48, 48, 50)>
  >>> sf_grid.get_value(3, 4, 5)  #doctest: +ELLIPSIS
  (54.5276...+53.4189...j)

In addition to `d_min` and `rate`, which govern the grid density,
DensityCalculator has two more parameters that affect the accuracy
of the calculated structure factors:

* `cutoff` (default: 1e-5) -- density cut-off in the same unit as the map.
  It is used to determine the atomic radius in which the density is calculated
  (density at the radius distance should be approximately `cutoff`).
  A smaller cutoff means more accurate but slower calculations.
* `blur` (default: 0) -- Gaussian dampening (blurring) factor --
  artificial temperature factor *B*\ :sub:`extra` added to all atomic B-factors
  (the structure factors must be later corrected to cancel it out).

.. _blur:

Choosing these parameters is a trade-off between efficiency and accuracy.
*B*\ :sub:`extra` is the most interesting one.
It is discussed in the `ITfC vol B <https://it.iucr.org/Bb/contents/>`_,
section 1.3.4.4.5 by G. Bricogne, and further in papers by
`J. Navaza (2002) <https://doi.org/10.1107/S0108767302016318>`_ and by
`P. Afonine and A. Urzhumtsev (2003) <https://doi.org/10.1107/S0108767303022062>`_,
but no formula for the optimal value exists.
The value of *B*\ :sub:`extra` that
gives the most accurate results depends on the resolution, oversampling,
atomic radius cut-off, and on the distribution of B-factors
(in particular, on the minimal B-factor in the model).
Additionally, increasing the dampening makes the computations slower
(because it increases the atomic radius).

*B*\ :sub:`extra` can be set explicitly (it can be negative):

.. doctest::

  >>> dencalc.blur = 10

or using the formula from Refmac (which is a function of the grid spacing
and *B*\ :sub:`min`). If the formula results in a negative number,
Refmac sets *B*\ :sub:`extra`\ to 0, and Servalcat, which employs the same
formula, uses the negative value:

.. doctest::

  >>> dencalc.set_refmac_compatible_blur(st[0], allow_negative=False)
  >>> dencalc.blur
  31.01648695048263

The :ref:`sfcalc <sfcalc>` program can be used to test different choices
of *B*\ :sub:`extra`.

If the density was blurred, you need to calculate it by either
adding option `unblur` to `prepare_asu_data()`:

.. doctest::

  >>> asu_data = grid.prepare_asu_data(dmin=dencalc.d_min, unblur=dencalc.blur)

or multiplying individual structure factor by
`dencalc.reciprocal_space_multiplier(inv_d2)`.

Mott-Bethe formula
~~~~~~~~~~~~~~~~~~

(It's a niche topic that most of the readers should skip.)

To calculate *f*\ :sub:`e` according to the Mott-Bethe formula
we first employ addends to calculate *f*\ :sub:`x`\ --\ *Z*:

.. doctest::

  >>> dc = gemmi.DensityCalculatorX()
  >>> dc.d_min = 2.5
  >>> dc.addends.subtract_z()
  >>> dc.grid.setup_from(st)
  >>> dc.set_refmac_compatible_blur(st[0])
  >>> dc.put_model_density_on_grid(st[0])
  >>> grid = gemmi.transform_map_to_f_phi(dc.grid)

Then we multiply it by –1/(2\ *π*:sup:`2`\ *a*:sub:`0`\ *d*:sup:`2`)
to get *f*\ :sub:`e`.

We either multiply individual values by `mott_bethe_factor()`
(which includes `reciprocal_space_multiplier()`):

.. doctest::

  >>> dc.mott_bethe_factor([3,4,5]) * grid.get_value(3,4,5)  #doctest: +ELLIPSIS
  (54.063...+52.970...j)

or we call `prepare_asu_data()` with `mott_bethe=True`:

.. doctest::
  :skipif: numpy is None

  >>> asu_data = grid.prepare_asu_data(dmin=2.5, mott_bethe=True, unblur=dencalc.blur)
  >>> asu_data.value_array[numpy.all(asu_data.miller_array == [3,4,5], axis=1)]  #doctest: +ELLIPSIS
  array([54.063...+52.970...j], dtype=complex64)

That is all.
If you would like to separate positions of hydrogen nuclei
and electron clouds then, assuming that the model
has positions for electrons, call `subtract_z()` as:

.. doctest::

  >>> dc.addends.subtract_z(except_hydrogen=True)

change hydrogen coordinates to proton positions,
and subtract Z=1 by adding c=-1:

.. doctest::

  >>> for cra in st[0].all():
  ...     if cra.atom.is_hydrogen():
  ...         dc.add_c_contribution_to_grid(cra.atom, -1)


.. _scaling:

Scaling and bulk solvent correction
-----------------------------------

Anisotropic scaling of calculated structure factors **F** and bulk solvent
correction are usually handled together because the bulk solvent parameters
are optimized along with the scaling parameters.
In Gemmi, both are implemented in one class: Scaling.

The bulk solvent correction is optional and used only with crystallographic
data. The solvent (not modeled as individual atoms) occupies a significant
volume of a macromolecular crystal, so accounting for it makes a difference.
If we have a model with atomic coordinates, we can create a mask,
as described in the section about :ref:`solvent masking <solventmask>`,
to model bulk solvent as a flat scatterer.
Then, we Fourier-transform the mask and obtain **F**\ :sub:`mask`.

The solvent correction is parametrized as:

  **F**\ :sub:`sol` = *k*\ :sub:`sol` exp(–\ *B*\ :sub:`sol` *s*\ :sup:`2`/4) **F**\ :sub:`mask`

*k*\ :sub:`sol` and *B*\ :sub:`sol` are refined together with the overall
anisotropic scaling parameters, *k*\ :sub:`ov` and **B**\ :sub:`ov`:

  **F**\ :sub:`tot` = *k*\ :sub:`ov` exp(–\ **s**\ :sup:`T` **B**\ :sub:`ov` **s** / 4) (**F**\ :sub:`cryst` + **F**\ :sub:`sol`)

**B**\ :sub:`ov` is a crystallographic anisotropic tensor that obeys
the symmetry; effectively, from a single fittable parameter for a cubic
crystal to 6 parameters for a primitive crystal.
All the parameters (*k*\ :sub:`ov`, **B**\ :sub:`ov`,
*k*\ :sub:`sol` and *B*\ :sub:`sol`), or only selected ones,
are optimized to make the total calculated amplitude
*F*\ :sub:`tot` = \|\ **F**\ :sub:`tot`\|
match the diffraction data *F*\ :sub:`obs`.
But according to what criterion?

First, there is least-squares scaling, which minimizes
∑(*F*\ :sub:`tot` – *F*\ :sub:`obs`)\ :sup:`2`.
While this target function has fallen out of favor in macromolecular
crystallography (because the errors are not really Gaussian-distributed),
it is still routinely used for the calculation of statistics such as R-factors.

The workhorse algorithm for non-linear least squares optimization
is the Levenberg-Marquardt method. Unlike other popular algorithms,
such as BFGS and conjugate gradient methods, L-M works only for least squares.
This is because while the other methods estimate the Hessian from the result of
the previous optimization step (BFGS) or ignore the second derivatives (CG),
L-M approximates the Hessian as a squared Jacobian (**J**\ :sup:`T` **J**),
which holds only for quadratic forms. However, for these forms L-M tends
to converge faster than the other methods.

Here is a full example. We read *F*\ :sub:`obs` from an mtz file,
then we read a corresponding model and use it to calculate
**F**\ :sub:`cryst` and **F**\ :sub:`mask`.
These functions were described in the previous sections.
Then we instantiate the Scaling class, pass the data points to it,
and finally we perform the least squares scaling.

**Note: functions in class Scaling are likely to change soon**.

.. doctest::

  >>> # read Fobs
  >>> mtz = gemmi.read_mtz_file('../tests/5e5z.mtz')
  >>> mtz.update_reso()
  >>> f_obs = mtz.get_value_sigma('FP', 'SIGFP')

  >>> # calculate Fcryst
  >>> st = gemmi.read_structure('../tests/5e5z.pdb')
  >>> d_min = mtz.resolution_high() - 1e-9
  >>> dc = gemmi.DensityCalculatorX()
  >>> dc.d_min = d_min
  >>> dc.set_refmac_compatible_blur(st[0])
  >>> dc.grid.setup_from(st)
  >>> dc.put_model_density_on_grid(st[0])
  >>> grid = gemmi.transform_map_to_f_phi(dc.grid, half_l=True)
  >>> f_cryst = grid.prepare_asu_data(dmin=d_min, unblur=dc.blur)

  >>> # calculate Fmask
  >>> mask_grid = gemmi.FloatGrid()
  >>> mask_grid.setup_from(st, spacing=min(0.6, d_min / 2 - 1e-9))
  >>> masker = gemmi.SolventMasker(gemmi.AtomicRadiiSet.Refmac)
  >>> masker.put_mask_on_float_grid(mask_grid, st[0])
  >>> fmask_gr = gemmi.transform_map_to_f_phi(mask_grid, half_l=True)
  >>> f_mask = fmask_gr.prepare_asu_data(dmin=d_min)

  >>> # now we can start with Scaling
  >>> scaling = gemmi.Scaling(st.cell, st.find_spacegroup())
  >>> scaling.use_solvent = True
  >>> scaling.prepare_points(f_cryst, f_obs, f_mask)
  >>> scaling.fit_isotropic_b_approximately()
  >>> wssr = scaling.fit_parameters()
  >>> print(f'RMSE: {(wssr / len(f_obs)) ** 0.5:.4f}')
  RMSE: 6.4176
  >>> print(f'R-factor: {scaling.calculate_r_factor():.2%}')
  R-factor: 17.73%

  >>> # get scaled Ftotal
  >>> f_tot = f_cryst.copy()
  >>> scaling.scale_data(f_tot, f_mask)

Least squares are sensitive to outliers. To make the scaling less sensitive,
we must change the target function. The absolute differences in R-factor are
less affected by outliers than the squared differences.
If we choose a target that's more similar to the R-factor,
we'll get more robust scaling, and also a lower R-factor.

This requires a different optimization algorithm. Actually, it's also possible
to *robustify* least squares and keep using (modified) Levenberg-Marquardt,
as was `reviewed <https://link.springer.com/chapter/10.1007/978-3-319-10602-1_50>`_
for the bundle adjustment procedure, but I haven't found any research indicating
these methods are preferable to simply using another algorithm.

TBC
