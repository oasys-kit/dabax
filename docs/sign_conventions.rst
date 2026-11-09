==========================================================
Sign conventions of structure factors and susceptibilities
==========================================================

Abstract
--------

The complex quantities that describe X-ray scattering by crystals (atomic scattering factors,
structure factors, susceptibilities, refractive index) change sign in their imaginary parts with the
time convention, :math:`e^{+i\omega t}` or :math:`e^{-i\omega t}`. We review the conventions used by
the atomic-data libraries DABAX and xraylib, which are the basis of many OASYS codes.

The structure factors and susceptibilities follow the :math:`e^{+i\omega t}` convention of the
International Tables, whereas the anomalous factor ``Fii`` and the refractive index follow
:math:`e^{-i\omega t}`.

1. Introduction
---------------

Many of the OASYS codes that compute X-ray interactions with matter and diffraction by crystals rely
on two atomic-data libraries: DABAX [1] and xraylib [2]. This note restricts itself to these two
libraries.

The literature uses two conventions for the sign of the imaginary unit, which appear in the phase
factor of the structure factor, :math:`e^{\pm2\pi i\mathbf{H}\cdot\mathbf{r}}`, and in the atomic
scattering factor, :math:`f_0+f'\pm i f''`. Mixing them gives wrong results that may pass unnoticed
in simple cases.

**Table 1.** The two time conventions. :math:`f''>0`, :math:`\delta` and :math:`\beta>0` are the
same in both.

+--------------------------------+--------------------------------+--------------------------------+
|                                | :math:`e^{+i\omega t}`         | :math:`e^{-i\omega t}`         |
|                                | (International Tables)         | (optics)                       |
+================================+================================+================================+
| plane wave                     | :math:`e^{i(\omega             | :math:`e^{i(\mathbf            |
|                                | t-\mathbf{k}\cdot\mathbf{r})}` | {k}\cdot\mathbf{r}-\omega t)}` |
+--------------------------------+--------------------------------+--------------------------------+
| :math:`f`                      | :math:`f_0+f'+if''`            | :math:`f_0+f'-if''`            |
+--------------------------------+--------------------------------+--------------------------------+
| phase factor in :math:`F`      | :math:`e^{+2\p                 | :math:`e^{-2\p                 |
|                                | i i\mathbf{H}\cdot\mathbf{r}}` | i i\mathbf{H}\cdot\mathbf{r}}` |
+--------------------------------+--------------------------------+--------------------------------+
| :math:`n`                      | :math:`1-\delta-i\beta`        | :math:`1-\delta+i\beta`        |
+--------------------------------+--------------------------------+--------------------------------+
| :math:`\chi_0`,                | :math:`\mathrm{Im}\,\chi_0<0`, | :math:`\mathrm{Im}\,\chi_0>0`, |
| :math:`\chi_H\chi_{\bar H}`    | :math:`\math                   | :math:`\math                   |
| (Si)                           | rm{Im}\,\chi_H\chi_{\bar H}>0` | rm{Im}\,\chi_H\chi_{\bar H}<0` |
+--------------------------------+--------------------------------+--------------------------------+

A monochromatic field can be written with a time dependence :math:`e^{+i\omega t}` or
:math:`e^{-i\omega t}`. Both describe the same physics, and the complex quantities of one convention
are the complex conjugates of those of the other (Table 1 summarizes the results that are discussed
in the text).

This note fixes the definitions and checks the conventions used for the atomic scattering factors
and the refractive index (Section 2) and for the structure factors and susceptibilities (Section 3).
All the numbers in the tables are reproduced by ``scripts/reproduce_tables.py``.

2. Atomic scattering factors and refractive index
-------------------------------------------------

2.1 The :math:`e^{+i\omega t}` convention (International Tables)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A plane wave is :math:`e^{i(\omega t - \mathbf{k}\cdot\mathbf{r})}`, and the atomic scattering
factor is [3]

.. math:: f(\mathbf{Q},E) = f_0(\mathbf{Q}) + f'(E) + i\,f''(E), \tag{1}

with :math:`f_0` the non-resonant form factor and :math:`f'`, :math:`f''` the anomalous corrections.
:math:`f''>0`, since it is proportional to the photoabsorption cross section,
:math:`\sigma_a = 2 r_e \lambda f''`; :math:`f'` is small and positive far above an absorption edge
and negative near it.

The refractive index of a medium with :math:`N_i = N_A\rho w_i/A_i` atoms of element :math:`i` per
unit volume (:math:`\rho` density, :math:`w_i` mass fraction, :math:`A_i` atomic weight) is

.. math:: n = 1 - \frac{r_e \lambda^2}{2\pi} \sum_i N_i\, f_i(0,E), \tag{2}

with :math:`f_i(0,E)=Z_i+f'_i+if''_i`. Writing, as usual for X-rays,

.. math:: n = 1 - \delta - i\beta, \tag{3}

we obtain

.. math:: \delta = \frac{r_e \lambda^2}{2\pi} \sum_i N_i\,(Z_i + f'_i), \qquad \beta = \frac{r_e \lambda^2}{2\pi} \sum_i N_i\,f''_i = \frac{\lambda}{4\pi}\,\mu_a, \tag{4}

where :math:`\mu_a` is the linear photoabsorption coefficient. The wave in the medium,
:math:`e^{i(\omega t - nkz)}`, decays as :math:`e^{-\beta k z}`.

An atom at :math:`\mathbf{r}_j` scatters with the phase factor
:math:`e^{+i\mathbf{Q}\cdot\mathbf{r}_j}`, :math:`\mathbf{Q}=\mathbf{k}'-\mathbf{k}`, which gives
:math:`e^{+2\pi i\mathbf{H}\cdot\mathbf{r}_j}` in the structure factor.

2.2 The :math:`e^{-i\omega t}` convention
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A plane wave is :math:`e^{i(\mathbf{k}\cdot\mathbf{r} - \omega t)}`, and all complex quantities are
conjugated:

.. math:: f^{(-)} = f_0 + f' - i\,f'', \qquad n^{(-)} = 1 - \delta + i\beta, \tag{5}

with the same :math:`f'`, :math:`f''>0`, :math:`\delta` and :math:`\beta`, and the phase factor
:math:`e^{-2\pi i\mathbf{H}\cdot\mathbf{r}_j}` in the structure factor.

2.3 Tabulations in DABAX
~~~~~~~~~~~~~~~~~~~~~~~~

DABAX provides several tabulations of :math:`f'` and :math:`f''` (Table 2). All store
:math:`f_2=f''>0`; :math:`f_1` is either :math:`Z+f'` (Windt [4], Henke [5], Chantler [6,7], EPDL97
[8], Kissel [9]) or :math:`f'` (Brennan–Cowan [10], Sasaki [11]). The methods ``Fi(Z,E)`` and
``Fii(Z,E)`` of ``DabaxXraylib`` reproduce the interface and the signs of xraylib: they subtract
:math:`Z` when needed and change the sign of :math:`f_2`,

.. math:: \texttt{Fi} = f', \qquad \texttt{Fii} = -f''.

Therefore :math:`f_0+\texttt{Fi}+i\,\texttt{Fii}` is the scattering factor in the
:math:`e^{-i\omega t}` convention, and

.. math:: f = f_0 + \texttt{Fi} - i\,\texttt{Fii} = f_0 + f' + i\,f'' \tag{6}

the one in the :math:`e^{+i\omega t}` convention, which is used by the DABAX structure-factor
routines. Code that uses :math:`+i\,\texttt{Fii}` together with
:math:`e^{+2\pi i\mathbf{H}\cdot\mathbf{r}}` changes the sign of the absorption. For Si at 8 keV the
tabulations differ by about 20% in :math:`f'` and 4% in :math:`f''`.

**Table 2.** Tabulated :math:`f_1,f_2` and values returned by ``Fi`` and ``Fii`` for Si at 8 keV.

+----------------+----------------+-------------+-------------+----------------+----------------+
| DABAX file     | :math:`f_1`    | :math:`f_1` | :math:`f_2` | ``Fi``         | ``Fii``        |
|                | column         |             |             | :math:`=f'`    | :math:`=-f''`  |
+================+================+=============+=============+================+================+
| ``f1           | :math:`Z+f'`   | 14.2983     | 0.3287      | +0.2983        | −0.3287        |
| f2_Windt.dat`` |                |             |             |                |                |
+----------------+----------------+-------------+-------------+----------------+----------------+
| ``f1           | :math:`Z+f'`   | 14.2983     | 0.3287      | +0.2983        | −0.3287        |
| f2_Henke.dat`` |                |             |             |                |                |
+----------------+----------------+-------------+-------------+----------------+----------------+
| ``f1f2_        | :math:`Z+f'`   | 14.2677     | 0.3338      | +0.2677        | −0.3338        |
| Chantler.dat`` |                |             |             |                |                |
+----------------+----------------+-------------+-------------+----------------+----------------+
| ``f1f          | :math:`Z+f'`   | 14.2582     | 0.3433      | +0.2582        | −0.3433        |
| 2_EPDL97.dat`` |                |             |             |                |                |
+----------------+----------------+-------------+-------------+----------------+----------------+
| ``f1f2_as      | :math:`Z+f'`   | 14.2605     | 0.3398      | +0.2605        | −0.3398        |
| f_Kissel.dat`` |                |             |             |                |                |
+----------------+----------------+-------------+-------------+----------------+----------------+
| ``f1f2_Bren    | :math:`f'`     | 0.2561      | 0.3340      | +0.2561        | −0.3340        |
| nanCowan.dat`` |                |             |             |                |                |
+----------------+----------------+-------------+-------------+----------------+----------------+
| ``f1f          | :math:`f'`     | 0.2455      | 0.3340      | +0.2455        | −0.3340        |
| 2_Sasaki.dat`` |                |             |             |                |                |
+----------------+----------------+-------------+-------------+----------------+----------------+
| xraylib        |                |             |             | +0.2584        | −0.3401        |
+----------------+----------------+-------------+-------------+----------------+----------------+

..

   ``Fi`` returns :math:`f'`, valid for both conventions (positive far above an absorption edge,
   negative near it, e.g., −9.1 for Si at 1.83 keV).

   ``Fii`` returns :math:`-f''`, which is negative, in agreement with the :math:`e^{-i\omega t}`
   convention.

2.4 Refractive index in DABAX
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The DABAX methods ``Refractive_Index`` (and ``_Re``, ``_Im``), like those of xraylib, return
:math:`n=1-\delta+i\beta`, with

.. math:: \delta = K_D \,\frac{\rho}{E^2} \sum_i \frac{w_i}{A_i}\,\bigl(Z_i + \texttt{Fi}_i\bigr), \qquad \beta = \frac{hc}{4\pi E}\, \rho \sum_i w_i \left(\frac{\mu}{\rho}\right)_{\mathrm{tot},i}, \tag{7}

:math:`K_D=r_e(hc)^2N_A/2\pi`, :math:`E` in keV. :math:`\delta` coincides with Eq. (4), but
:math:`\beta` uses the *total* attenuation coefficient (including Rayleigh and Compton scattering,
file ``CrossSec_EPDL97.dat``) instead of :math:`f''`. Table 3 shows that :math:`\beta` then exceeds
the value from :math:`f''` by 5% for Si at 8 keV and by about 25% for Si at 30 keV and Be at 8 keV,
where scattering is a large part of the attenuation.

**Table 3.** :math:`\delta` and :math:`\beta` computed by DABAX (``f1f2_Windt.dat``,
``CrossSec_EPDL97.dat``) and xraylib, and :math:`\beta` from :math:`f''` [Eq. (4)] and from the
photoelectric attenuation coefficient :math:`\mu_\mathrm{ph}`. Densities (g cm⁻³): Si 2.33, Be
1.848, SiO₂ 2.2, Pt 21.45.

+-----------+---------------+-----------+-----------+-----------+-----------+-----------+-----------+
|           | :math:`E`     | :math     | :math     | :mat      | :mat      | :mat      | :mat      |
|           | (keV)         | :`\delta` | :`\delta` | h:`\beta` | h:`\beta` | h:`\beta` | h:`\beta` |
|           |               | DABAX     | xraylib   | DABAX     | from      | :math     | xraylib   |
|           |               |           |           |           | :m        | :`\mu_\ma |           |
|           |               |           |           |           | ath:`f''` | thrm{ph}` |           |
+===========+===============+===========+===========+===========+===========+===========+===========+
| Si        | 8.0           | 7.695e-6  | 7.672e-6  | 1.858e-7  | 1.769e-7  | 1.831e-7  | 1.859e-7  |
+-----------+---------------+-----------+-----------+-----------+-----------+-----------+-----------+
| Si        | 1.8           | 1.046e-4  | 9.893e-5  | 4.165e-6  | 4.062e-6  | 4.153e-6  | 4.184e-6  |
+-----------+---------------+-----------+-----------+-----------+-----------+-----------+-----------+
| Si        | 1.9           | 1.024e-4  | 1.031e-4  | 3.784e-5  | 3.790e-5  | 3.782e-5  | 3.668e-5  |
+-----------+---------------+-----------+-----------+-----------+-----------+-----------+-----------+
| Si        | 30.0          | 5.365e-7  | 5.367e-7  | 1.102e-9  | 8.743e-10 | 8.889e-10 | 1.101e-9  |
+-----------+---------------+-----------+-----------+-----------+-----------+-----------+-----------+
| Be        | 8.0           | 5.326e-6  | 5.328e-6  | 2.568e-9  | 2.074e-9  | 2.012e-9  | 2.561e-9  |
+-----------+---------------+-----------+-----------+-----------+-----------+-----------+-----------+
| SiO₂      | 8.0           | 7.222e-6  | 7.211e-6  | 9.880e-8  | 9.429e-8  | 9.692e-8  | 9.885e-8  |
+-----------+---------------+-----------+-----------+-----------+-----------+-----------+-----------+
| Pt        | 8.0           | 5.257e-5  | 5.237e-5  | 5.238e-6  | 5.211e-6  | 5.103e-6  | 5.257e-6  |
+-----------+---------------+-----------+-----------+-----------+-----------+-----------+-----------+

..

   The refractive index returned by DABAX (``Refractive_Index``) and by xraylib has
   :math:`\mathrm{Im}\,n = \beta > 0`: it complies with the :math:`e^{-i\omega t}` convention, Eq.
   (5).

3. Crystal structure factor and susceptibility
----------------------------------------------

3.1 Structure factor
~~~~~~~~~~~~~~~~~~~~

In the :math:`e^{+i\omega t}` convention,

.. math:: F(\mathbf{H}) = \sum_j f_j\, T_j\, e^{+2\pi i\,\mathbf{H}\cdot\mathbf{r}_j}, \tag{8}

with :math:`f_j` evaluated at :math:`s=\sin\theta_B/\lambda=1/(2d_H)` and the Debye–Waller factor
:math:`T_j = e^{-B_j s^2}`. Dynamical diffraction needs :math:`F_0=F(\mathbf{0})`, with
:math:`T_j=1` and :math:`\mathrm{Im}\,F_0>0`, :math:`F_H=F(\mathbf{H})` and
:math:`F_{\bar H}=F(-\mathbf{H})`.

In the :math:`e^{-i\omega t}` convention :math:`F^{(-)}(\mathbf{H}) = F(\mathbf{H})^*`, with the
same moduli.

The phases of :math:`F_H` and :math:`F_{\bar H}` depend on the origin of the unit cell; :math:`F_0`,
:math:`\lvert F_H\rvert`, :math:`\lvert F_{\bar H}\rvert` and :math:`F_HF_{\bar H}` do not, and only
these can be compared between codes. For the diamond structure with origin at an atom,
:math:`F(hkl) = f\, G_{fcc}\,[1 + e^{i\pi(h+k+l)/2}]`, with :math:`G_{fcc}=4` for unmixed
:math:`hkl` and 0 otherwise: :math:`F_{111}=4f(1-i)`, :math:`F_{\bar1\bar1\bar1}=4f(1+i)`, and
:math:`F_{111}=F_{\bar1\bar1\bar1}=-4\sqrt2 f` with the origin at the inversion centre.

DABAX (``Crystal_F_0_F_H_F_H_bar_StructureFactor``) and xraylib (``Crystal_F_H_StructureFactor``)
evaluate Eq. (8) with Eq. (6). Their results for Si (Table 4) agree with the analytical expression
to :math:`10^{-13}`. The Debye–Waller factor is passed as a single number for the reflection
:math:`\mathbf{H}`; xraylib applies it also when called with :math:`hkl=000`, so the caller must
pass 1 for :math:`F_0`. DABAX applied it to :math:`F_0` as well; this error has been corrected
(October 2026).

**Table 4.** Structure factors of Si at 8 keV (DABAX, ``f1f2_Windt.dat``, ``f0_InterTables.dat``,
:math:`T=1`) for the cells ``Si`` (origin at an atom) and ``Si2`` (origin at the inversion centre).
:math:`F_0 = 114.3669+2.6295i`.

+---------+-------------+-----------------+-----------------+-----------------+-----------------+
| cell    | :math:`hkl` | :math:`F_H`     | :ma             | :math:`\l       | :math:`         |
|         |             |                 | th:`F_{\bar H}` | vert F_H\rvert` | F_H F_{\bar H}` |
+=========+=============+=================+=================+=================+=================+
| ``Si``  | 111         | :math:`44       | :math:`42       | 61.3149         | :math:`3        |
|         |             | .6510-42.0215i` | .0215+44.6510i` |                 | 752.61+227.91i` |
+---------+-------------+-----------------+-----------------+-----------------+-----------------+
| ``Si2`` | 111         | :math:`-6       | :math:`-6       | 61.3149         | :math:`3        |
|         |             | 1.2867-1.8594i` | 1.2867-1.8594i` |                 | 752.61+227.91i` |
+---------+-------------+-----------------+-----------------+-----------------+-----------------+
| ``Si``  | 220         | :math:`7        | :math:`7        | 72.1760         | :math:`5        |
|         |             | 2.1281+2.6295i` | 2.1281+2.6295i` |                 | 195.54+379.33i` |
+---------+-------------+-----------------+-----------------+-----------------+-----------------+
| ``Si2`` | 220         | :math:`-7       | :math:`-7       | 72.1760         | :math:`5        |
|         |             | 2.1281-2.6295i` | 2.1281-2.6295i` |                 | 195.54+379.33i` |
+---------+-------------+-----------------+-----------------+-----------------+-----------------+

..

   The structure factors returned by DABAX (``Crystal_F_0_F_H_F_H_bar_StructureFactor``,
   ``Crystal_F_H_StructureFactor``) and xraylib (``Crystal_F_H_StructureFactor`` with :math:`hkl`,
   :math:`\bar h\bar k\bar l`, :math:`000`) use :math:`f_0+\texttt{Fi}-i\,\texttt{Fii}` and
   :math:`e^{+2\pi i\mathbf{H}\cdot\mathbf{r}}`, and have :math:`\mathrm{Im}\,F_0>0`: they comply
   with the :math:`e^{+i\omega t}` convention, Eq. (8).

3.2 Friedel’s law and Bijvoet pairs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Neglecting the imaginary part :math:`f''` (no absorption), the atomic scattering factors are real
and Eq. (8) gives :math:`F(-\mathbf{H}) = F(\mathbf{H})^*` for any structure. This is Friedel’s law
in its strict form; it implies equal intensities,
:math:`\lvert F(\mathbf{H})\rvert=\lvert F(-\mathbf{H})\rvert`, so that diffraction cannot
distinguish a crystal from its inverted image. Since :math:`f''>0` for all atoms, the strict form
never holds exactly: in general :math:`F(-\mathbf{H})\neq F(\mathbf{H})^*`.

The equality of intensities, the form in which Friedel’s law is usually stated, still holds with
absorption in two cases. In a centrosymmetric crystal every atom at :math:`\mathbf{r}` has an
identical partner at :math:`-\mathbf{r}` (origin at the inversion centre), so that, summing over the
pairs,
:math:`F(\mathbf{H}) = \sum_j 2 f_j T_j \cos(2\pi\mathbf{H}\cdot\mathbf{r}_j) = F(-\mathbf{H})`, and
:math:`\lvert F(\mathbf{H})\rvert=\lvert F(-\mathbf{H})\rvert` for any origin. In a crystal of a
single element, centrosymmetric or not, :math:`F(\pm\mathbf{H}) = f\,G(\pm\mathbf{H})` with
:math:`G(-\mathbf{H})=G(\mathbf{H})^*`, so that
:math:`\lvert F(\mathbf{H})\rvert=\lvert F(-\mathbf{H})\rvert`; more generally, this holds whenever
all atoms have the same ratio :math:`f''/(f_0+f')`.

In a non-centrosymmetric crystal with atoms of different kinds, instead,
:math:`\lvert F(\mathbf{H})\rvert\neq\lvert F(-\mathbf{H})\rvert`: the reflections
:math:`\mathbf{H}` and :math:`-\mathbf{H}` form a Bijvoet pair with different intensities, a
difference used to determine the absolute structure of polar crystals. For GaAs (zinc-blende
structure, Ga at the origin, As at :math:`(\tfrac14,\tfrac14,\tfrac14)`),
:math:`F(111)=4(f_\mathrm{Ga}-i f_\mathrm{As})` and
:math:`F(\bar1\bar1\bar1)=4(f_\mathrm{Ga}+i f_\mathrm{As})`, so that

.. math:: \lvert F(111)\rvert^2-\lvert F(\bar1\bar1\bar1)\rvert^2 = 64\,\bigl(A\,f''_\mathrm{As} - B\,f''_\mathrm{Ga}\bigr), \qquad A=f_{0,\mathrm{Ga}}+f'_\mathrm{Ga},\quad B=f_{0,\mathrm{As}}+f'_\mathrm{As}, \tag{9}

which vanishes if :math:`f''=0` or if the two atoms are identical (diamond structure). At 8 keV
(:math:`f''_\mathrm{Ga}=0.81`, :math:`f''_\mathrm{As}=1.06`) the two moduli differ by 0.7% (Table
5).

xraylib gives 149.396 and 148.416, with the same ordering.

**Table 5.** Friedel’s law and Bijvoet pairs for the 111 reflection of Si (centrosymmetric) and GaAs
(non-centrosymmetric) at 8 keV, computed with DABAX (``f1f2_Windt.dat``, ``f0_InterTables.dat``,
:math:`T=1`, origin at an atom).

+--------------------------------------------+------------------------+----------------------------+
|                                            | Si (centrosymmetric)   | GaAs (non-centrosymmetric) |
+============================================+========================+============================+
| :math:`F(111)`                             | :math:`44.651-42.021i` | :math:`105.733-105.674i`   |
+--------------------------------------------+------------------------+----------------------------+
| :math:`F(\bar1\bar1\bar1)`                 | :math:`42.021+44.651i` | :math:`97.274+112.160i`    |
+--------------------------------------------+------------------------+----------------------------+
| :math:`F(111)^*`                           | :math:`44.651+42.021i` | :math:`105.733+105.674i`   |
+--------------------------------------------+------------------------+----------------------------+
| :math:`F(\bar1\bar1\bar1)=F(111)^*`        | no                     | no                         |
| (strict form)                              |                        |                            |
+--------------------------------------------+------------------------+----------------------------+
| :math:`\lvert F(\                          | yes                    | no                         |
| bar1\bar1\bar1)\rvert=\lvert F(111)\rvert` |                        |                            |
| (equal intensities)                        |                        |                            |
+--------------------------------------------+------------------------+----------------------------+
| :math:`\lvert F(111)\rvert`,               | 61.315, 61.315         | 149.487, 148.466           |
| :math:`\lvert F(\bar1\bar1\bar1)\rvert`    |                        |                            |
+--------------------------------------------+------------------------+----------------------------+
| same, with :math:`f''=0`                   | 61.287, 61.287         | 148.882, 148.882           |
+--------------------------------------------+------------------------+----------------------------+

3.3 Susceptibilities
~~~~~~~~~~~~~~~~~~~~

The susceptibilities (Zachariasen’s :math:`\psi` [16]) are

.. math:: \chi_H = -\Gamma F_H , \qquad \Gamma = \frac{r_e\lambda^2}{\pi V}, \tag{10}

and similarly :math:`\chi_0` and :math:`\chi_{\bar H}`, with :math:`V` the cell volume. Since
:math:`F_0/V=\sum_iN_if_i(0,E)`, :math:`n=1+\chi_0/2` [Eq. (2)].

In the :math:`e^{+i\omega t}` convention

.. math:: \mathrm{Re}\,\chi_0 < 0, \qquad \mathrm{Im}\,\chi_0 < 0 . \tag{11}

In the :math:`e^{-i\omega t}` convention the signs of all the imaginary parts are reversed.

The two-beam diffraction profile depends on :math:`\chi_0` and :math:`\chi_H\chi_{\bar H}`, which do
not depend on the origin.

A simple test applies to a centrosymmetric crystal of one element (Si): with the origin at the
inversion centre, :math:`F_H=F_{\bar H}=GfT` with :math:`G` real, so
:math:`\mathrm{Im}(\chi_H\chi_{\bar H}) = 2\Gamma^2G^2T^2\,f''\,\mathrm{Re}f`. Provided that
:math:`\mathrm{Re}f=f_0+f'>0`, which holds for Si except very close to the K edge (1.839 keV), where
:math:`f'` reaches −11.2 (``f1f2_Windt.dat``) and :math:`\lvert f'\rvert` can exceed
:math:`f_0\geq10.8` for reflections near backscattering, this gives

.. math:: \mathrm{Im}\,\chi_0 < 0, \qquad \mathrm{Im}(\chi_H\chi_{\bar H}) > 0 , \tag{12}

both reversed in the :math:`e^{-i\omega t}` convention. If both have the same sign, the set mixes
the two conventions.

3.4 Crystal susceptibility and refractive index
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For Si at 8 keV and the density of the unit cell (2.32945 g cm⁻³), DABAX gives
:math:`1+\chi_0/2 = 0.999992308 - 1.7686\times10^{-7}\,i` and ``Refractive_Index``
:math:`= 0.999992307 + 1.8577\times10^{-7}\,i`. The real parts agree to :math:`10^{-9}`; the
imaginary parts do not:

-  they have opposite signs (the two conventions);
-  their absolute values are different: the first one uses the photoelectric attenuation only, from
   :math:`f''`, whereas the second one uses the total attenuation (photoelectric, Compton and
   Rayleigh; Section 2.4).

The two quantities must not be combined without conjugating one of them.

4. Conclusions
--------------

-  DABAX and xraylib compute structure factors in the :math:`e^{+i\omega t}` convention of the
   International Tables, with :math:`f=f_0+\texttt{Fi}-i\,\texttt{Fii}`, while
   :math:`\texttt{Fii}=-f''` and the refractive index (:math:`\mathrm{Im}\,n>0`) follow
   :math:`e^{-i\omega t}`. The refractive index uses the total attenuation for :math:`\beta`, up to
   25% larger than the value from :math:`f''`.
-  Equation (12) is a simple test of the consistency of a set of susceptibilities: both signs must
   be reversed together when changing convention.
-  Friedel’s law holds strictly only for :math:`f''=0`; with absorption, the equality of intensities
   survives in centrosymmetric and single-element crystals, while GaAs shows a Bijvoet difference of
   0.7% for the 111 reflection at 8 keV.

References
----------

1.  M. Sanchez del Rio, *DABAX: a database for X-ray applications*,
    https://github.com/oasys-kit/dabax.
2.  T. Schoonjans, A. Brunetti, B. Golosio, M. Sanchez del Rio, V. A. Solé, C. Ferrero, L. Vincze,
    The xraylib library for X-ray–matter interactions. Recent developments, *Spectrochim. Acta B*
    **66**, 776–784 (2011).
3.  E. Prince (ed.), *International Tables for Crystallography, Vol. C*, 3rd ed., IUCr (2004).
4.  D. L. Windt, IMD — software for modeling the optical properties of multilayer films, *Computers
    in Physics* **12**, 360–370 (1998).
5.  B. L. Henke, E. M. Gullikson, J. C. Davis, X-ray interactions: photoabsorption, scattering,
    transmission, and reflection at E = 50–30,000 eV, Z = 1–92, *At. Data Nucl. Data Tables* **54**,
    181–342 (1993).
6.  C. T. Chantler, Theoretical form factor, attenuation, and scattering tabulation for Z = 1–92
    from E = 1–10 eV to E = 0.4–1.0 MeV, *J. Phys. Chem. Ref. Data* **24**, 71–643 (1995).
7.  C. T. Chantler, Detailed tabulation of atomic form factors, photoelectric absorption and
    scattering cross section, and mass attenuation coefficients in the vicinity of absorption edges
    in the soft X-ray (Z = 30–36, Z = 60–89, E = 0.1 keV–10 keV), addressing convergence issues of
    earlier work, *J. Phys. Chem. Ref. Data* **29**, 597–1048 (2000).
8.  D. E. Cullen, J. H. Hubbell, L. Kissel, *EPDL97: the evaluated photon data library, ’97
    version*, LLNL report UCRL-50400, Vol. 6, Rev. 5 (1997).
9.  L. Kissel, B. Zhou, S. C. Roy, S. K. Sen Gupta, R. H. Pratt, The validity of form-factor,
    modified-form-factor and anomalous-scattering-factor approximations in elastic scattering
    calculations, *Acta Cryst. A* **51**, 271–288 (1995).
10. S. Brennan, P. L. Cowan, A suite of programs for calculating X-ray absorption, reflection, and
    diffraction performance for a variety of materials at arbitrary wavelengths, *Rev. Sci.
    Instrum.* **63**, 850–853 (1992).
11. S. Sasaki, *Numerical tables of anomalous scattering factors calculated by the Cromer and
    Liberman’s method*, KEK Report 88-14 (1989).
12. W. H. Zachariasen, *Theory of X-Ray Diffraction in Crystals*, Wiley, New York (1945).
