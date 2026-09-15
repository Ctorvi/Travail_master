# This is the repository for my Master's thesis "[Error-field induced magnetic tangles in TCV (PDF)](V_Despland_master_project_compressed.pdf)".

This project consisted of calculating the homoclinic tangles in the TCV divertor, and correlate these calculations with MANTIS observations of anomalous light emmission in this region.
This work was produced by Victor Léonard Despland.
The analysis consists of extracting the equilibrium axisymmetric field from a TCV discharge at a specific time from [LIUQUE](Papers/Theory chapter/LIUQE or psi axi/LIUQE_Theory.pdf), extracting the :math:n=1n=1n=1 component of the non-axisymmetric perturbation due to coil displacements using a modification of Piras's method. Then compute the Poincaré plot and the manifolds using the method implemented in PyOculus. An homoclinic tangle is forming around the divertor X-point.
The results can be plotted on top of the tomographically reconstructed light emission from MANTIS, and the correlation between the tangle and the anomalous light emission can be studied.


![Poincaré](figures/MANTIS_MF_He.png)
