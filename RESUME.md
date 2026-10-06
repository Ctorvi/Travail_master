# Repository for my Master’s thesis: “[Error-field induced magnetic tangles in TCV (PDF)](V_Despland_master_project_compressed.pdf)”

This work was carried out by Victor Léonard Despland under the supervision of Dr. Christopher Smiet and Dr. Joaquim Loizu.

This project investigated error-field-induced homoclinic tangles in the divertor region of the TCV tokamak at EPFL, Lausanne. The aim was to investigate whether these magnetic structures could explain unusual light-emission patterns observed by the MANTIS diagnostic in configurations with a secondary X-point. Figure 1 shows three “jellyfish” discharges exhibiting tentacle-like emission structures. 

![MANTIS tomographic reconstructions](figures/MANTIS_MF_intro.png)

*Figure 1. Tomographic reconstructions of three “jellyfish” discharges showing tentacle-like emission structures in the divertor region.*

The working hypothesis was that these structures arise from magnetic tangles induced by error fields associated with non-axisymmetric displacements of the poloidal field (PF) coils.

## 1. Estimation of non-axisymmetric PF coil displacements

The coil displacements were estimated using the method described by [F. Piras et al.](https://www.sciencedirect.com/science/article/abs/pii/S0920379610001900). This approach uses a first-order model of coil displacements and a least-squares fit to magnetic measurements from the saddle-loop system.

For each of the 18 PF coils - 16 plasma shaping coil (E, F) and 2 Ohmic Heating system coils (OH) - four displacement parameters were considered: two centre shifts along the X and Y axes, and two tilts about the X and Y axes.

To isolate the contribution of each coil, measurements from plasma-free “backoff” shots were used, with current flowing through only one coil per shots. The estimated shifts and tilts are presented in Figure 2.

<p align="center">
  <img src="figures/Final_DX_78000.png" width="48%" alt="Polar plot of estimated PF coil centre shifts">
  <img src="figures/Final_tilts_78000.png" width="48.5%" alt="Polar plot of estimated PF coil tilts">
</p>

<p align="center">
  <em>Figure 2. Estimated PF coil displacements. Left: X and Y centre shifts. Right: tilts about the X and Y axes.</em>
</p>

## 2. Computation of a perturbed TCV magnetic configuration

Once the displacements are computed, the corresponding ('n = 1`) non-axisymmetric magnetic perturbations are added on top of the axisymmetric equilibrium, provided by [LIUQE](https://www.researchgate.net/publication/270825456_Tokamak_equilibrium_reconstruction_code_LIUQE_and_its_real_time_implementation) 
To analyze the 3-D magnetic topology, Poincaré plots of the magnetic field trajectories are computed. Existing tools from the [pyoculus](https://pypi.org/project/pyoculus/) package are adapted, for working on TCV experimental datas. 
An example of a non-pertubed and a pertubed poincaré plot of a TCV discharge are shown on figure 3.

<p align="center">
  <img src="script/jellyfisch/75979_12/figures/JF_75979_NP_n80_TRY.png" width="48%" alt="Polar plot of estimated PF coil centre shifts">
  <img src="script/jellyfisch/75979_12/figures/JF_75979_P_n80_TRY.png" width="35.6%" alt="Polar plot of estimated PF coil tilts">
</p>

<p align="center">
  <em>Figure 3. Poincaré plots of a Jellyfish discharge, using LIUQE and PyOculus. Left: Non-pertubed and axisymmetric magnetic configuration.  Right: Pertubed and NON-axisymmetric magnetic configuration.</em>
</p>


Poincaré plots and invariant manifolds are subsequently computed using the methods implemented in [pyoculus](https://pypi.org/project/pyoculus/). These manifolds form a homoclinic tangle around the divertor X-point.

The resulting magnetic structures can be overlaid on tomographic reconstructions of the MANTIS light emission, allowing the spatial correlation between the homoclinic tangle and the anomalous emission patterns to be investigated.

More detailed documentation of the numerical routine is available [here](https://gitlab.epfl.ch/spc/tcv/analysis/tokatangle).

![Poincaré](figures/MANTIS_MF_He.png)

*Figure 3 — Error-field tangles for 4 Jellyfish shots and one standard shot. On top,
poincaré plots with manifolds/ fixed-points. On the middle, MANTIS/ Manifold plot
zoomed in the divertor region. On bottom, simple MANTIS diagnostic plots. The
appearance and the spatial location of the tentacles seem to correlate with the manifolds*
