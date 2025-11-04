from pyoculus.fields import AxisymmetricCylindricalGridField
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
import numpy as np
from matplotlib import pyplot as plt
import logging
from scipy.io import loadmat

data = loadmat('./Script/Snowflake/Liuqe_SF_09.mat', squeeze_me=True, struct_as_record=False)


flux=data['flux']
r = flux.r
z = flux.z
psi = flux.psi.T
psi_r = flux.psi_r
psi_z = flux.psi_z.T
psi_axis = flux.psi_axis



fig, ax = plt.subplots(1, 1, figsize=(6, 6.6))



levels1 = np.linspace(-flux.psi_axis, 0, 15)
levels2 = np.linspace(flux.psi_axis, 0, 15)

# Tri pour éviter l'erreur
levels1.sort()
levels2.sort()

# Premier contour
plt.contour(
    psi_r, psi_z, psi,
    levels=levels1,
    colors='k', linestyles='--',
    linewidths=1

)

# Deuxième contour
plt.contour(
    psi_r, psi_z, psi,
    levels=levels2,
    colors='k',
    linewidths=1
)

# Troisième contour : niveau 0, trait plus épais
plt.contour(
    psi_r, psi_z, psi,
    levels=[0],
    colors='k',
    linewidths=1
)

plt.xlabel("R (m)")
plt.ylabel("Z (m)")
plt.axis('equal')
plt.show()


