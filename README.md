# Repository for my Master’s thesis: “[Error-field induced magnetic tangles in TCV (PDF)](V_Despland_master_project_compressed.pdf)”

This project focused on computing homoclinic tangles in the divertor region of the TCV tokamak at EPFL in Lausanne and comparing these calculations with MANTIS observations of anomalous light emission in this region.

This work was carried out by Victor Léonard Despland.

The analysis consists of extracting the axisymmetric equilibrium magnetic field of a TCV discharge at a specific time using [LIUQE](https://www.researchgate.net/publication/270825456_Tokamak_equilibrium_reconstruction_code_LIUQE_and_its_real_time_implementation). The `n = 1` component of the non-axisymmetric magnetic perturbation caused by coil displacements is then calculated using a modified version of [Piras et al.’s method](https://www.sciencedirect.com/science/article/abs/pii/S0920379610001900).

Poincaré plots and invariant manifolds are subsequently computed using the methods implemented in [pyoculus](https://pypi.org/project/pyoculus/). These manifolds form a homoclinic tangle around the divertor X-point.

The resulting magnetic structures can be overlaid on tomographic reconstructions of the MANTIS light emission, allowing the spatial correlation between the homoclinic tangle and the anomalous emission patterns to be investigated.

More detailed documentation of the numerical routine is available [here](https://gitlab.epfl.ch/spc/tcv/analysis/tokatangle).

![Poincaré](figures/MANTIS_MF_He.png)

*Figure 1 — Error-field tangles for 4 Jellyfish shots and one standard shot. On top,
poincaré plots with manifolds/ fixed-points. On the middle, MANTIS/ Manifold plot
zoomed in the divertor region. On bottom, simple MANTIS diagnostic plots. The
appearance and the spatial location of the tentacles seem to correlate with the manifolds*
