
from pyoculus.fields import AnalyticCylindricalBfield
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
from pyoculus.maps import TokaMap
from pyoculus.solvers.resonnance_zone import ResonnanceZone
from matplotlib import pyplot as plt
import numpy as np
import logging


toncaca= TokaMap(K=2.6,w=1.5)

Npoint=200
Niter=4000

psi= np.linspace(0, 1, Npoint)
theta= np.ones(psi.shape)*0.55

y=np.column_stack([theta, psi])
coords = np.zeros((Npoint, Niter+1, 2))  # (nombre de points, nombre d'itérations+1, 2)
coords[:, 0, :] = y  # point initial

for i in range(Npoint):
    for j in range(Niter):
        coords[i, j+1, :] = toncaca._f(coords[i, j, :])


opointL = FixedPoint(toncaca)
opointL.find(t=1, guess=[0.55, 0.73], method='scipy.root')
opointL.plot(ax=None, marker='o', color="x")

xpointL= FixedPoint(toncaca)
xpointL.find(t=1, guess=[0.28, 0.85], method='scipy.root')
x_pointL_coord = xpointL.coords[0]





fig, ax4 = plt.subplots(figsize=(6, 6))
ax4.scatter(coords[:, :, 0], coords[:, :, 1], s=0.1, color='black')
ax4.set_xlim(0.25, 0.85)
ax4.set_ylim(0.4, 1.0)
ax4.set_xlabel(r'$\psi$')
ax4.set_ylabel(r'$\theta$')
ax4.set_aspect('equal', 'box')
ax4.set_title('TokaMap: Many island chains')
plt.show()
# Résonances