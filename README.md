# Error-field-induced magnetic tangles in TCV

This repository contains code and results from my [master’s thesis (PDF)](V_Despland_master_project_compressed.pdf), carried out at the Swiss Plasma Center, EPFL, under the supervision of Dr. Christopher Smiet and Dr. Joaquim Loizu.

The project investigated homoclinic tangles in the divertor region of the TCV tokamak and their spatial correspondence with anomalous light-emission patterns observed by the MANTIS diagnostic.

- [Project overview](RESUME.md): methods, results and my contributions.
- [Numerical documentation and runnable example](https://gitlab.epfl.ch/spc/tcv/analysis/tokatangle): documentation of the numerical workflow. (EPFL access required)

## Method

The axisymmetric magnetic equilibrium of a TCV discharge at a selected time was combined with the $n = 1$ component of the non-axisymmetric magnetic perturbation induced by poloidal field (PF) coil displacements. These displacements were estimated from magnetic measurements using an adapted implementation of an [existing TCV method](https://www.sciencedirect.com/science/article/abs/pii/S0920379610001900).

Poincaré plots and invariant manifolds were then computed using tools from [pyoculus](https://pypi.org/project/pyoculus/), adapted to TCV magnetic configurations. The computed manifolds formed homoclinic tangles near the divertor X-point.

## Main result

The computed magnetic manifolds were overlaid on tomographic reconstructions of the MANTIS light emission. The comparison revealed a qualitative spatial correspondence between the manifolds and the observed tentacle-like emission structures, supporting further investigation of the proposed turnstile transport mechanism.

![Comparison of computed magnetic manifolds and MANTIS emission](figures/MANTIS_MF_He.png)

*Figure 1. Comparison of computed magnetic manifolds and MANTIS emission for four “jellyfish” discharges and one standard divertor discharge. Top row: Poincaré plots showing invariant manifolds and fixed points. Middle row: manifolds overlaid on MANTIS tomographic reconstructions, with a close-up of the divertor region. Bottom row: corresponding MANTIS reconstructions without overlays. The tentacle-like emission structures show a qualitative spatial correspondence with the computed manifolds.*