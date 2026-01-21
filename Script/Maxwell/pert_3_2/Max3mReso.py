
from pyoculus.fields import AnalyticCylindricalBfield
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
from pyoculus.solvers.resonnance_zone import ResonnanceZone
from matplotlib import pyplot as plt
import numpy as np
import logging


plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.size": 10,
    }
)
logging.basicConfig(level=logging.DEBUG)

POINCARE_ITS = 500
SHEAR =  0.4
SF = 0.8875

maxwellboltzmann = {
    "m": 3,
    "n": 2,
    "d": 1.75/np.sqrt(2),
    "type": "maxwell-boltzmann",
    "amplitude": 0.1,
}

Hits=np.load('Script/ResonnanceZone/Reso_Hits.npy')

separatrix = {"type": "circular-current-loop", "amplitude": -10, "R": 6, "Z": -5.5}

fig, ax = plt.subplots(1, 1, figsize=(6, 6.6))



perturbedfield = AnalyticCylindricalBfield.with_new_axis(
    R=6, Z=0, sf=SF, shear=SHEAR, perturbations_args=[maxwellboltzmann]
)
perturbedmap = CylindricalBfieldSection.without_axis(perturbedfield, guess=[6., 0])
axiscoords = np.array((perturbedmap.R0, perturbedmap.Z0))


opointT = FixedPoint(perturbedmap)
opointT.find(t=3, guess=[7.6, 0.], method='scipy.root')
opointT.plot(ax=ax, marker='o', color="xkcd:green")
opointTcoords = opointT.coords[0]

xpointL = FixedPoint(perturbedmap)
xpointL.find(t=3, guess=[6.4, -1.6], method='scipy.root')
xpointL.plot(ax=ax, marker='x', color="xkcd:green")
xpointLcoords = xpointL.coords[0]

xpointT = FixedPoint(perturbedmap)
xpointT.find(t=3, guess=[6.4, 1.6], method='scipy.root')
xpointT.plot(ax=ax, marker='x', color="xkcd:green")
xpointTcoords = xpointT.coords[0]




ResoManifoldDown = Manifold(perturbedmap, xpointL, xpointT, opointTcoords-xpointLcoords, opointTcoords-xpointTcoords)
ResoManifoldUp = Manifold(perturbedmap, xpointT, xpointL, opointTcoords-xpointTcoords, opointTcoords-xpointLcoords)


REZO3=ResonnanceZone(ResoManifoldDown,ResoManifoldUp)
contour_points=REZO3.contour()

FluxMeiss=REZO3.area()
FluxMP=ResoManifoldUp._AdL_integral_points(contour_points,dl=None, is_closed=True)
FluxSL=REZO3.flux_approx_SL()


REZO3.plot(with_clinics=False, ax=ax, markersize=0, lw=1.5)

ax.scatter(Hits[:,:, 0], Hits[:,:, 1], color="xkcd:dark grey", s=1.4, linewidths=0)
texte = (
    r"$\phi_{meiss}=$" + f"{FluxMeiss:.5f}\n"
    r"$\phi_{MP}=$" + f"{FluxMP:.5f}\n"
    r"$\phi_{SL}=$" + f"{FluxSL:.5f}"
)
ax.text(6.2, -0.1, texte, fontsize=12, color='blue')  # moved to south-east

ax.set_xlabel(r"$R[m]$")
ax.set_ylabel(r"$Z[m]$")
ax.set_xlim(6, 8)
ax.set_ylim(-1.9, 1.9)
ax.set_aspect('equal') 

ax.plot(contour_points[:,0], contour_points[:,1], color='black', linewidth=0.8)

#plt.savefig('Doublet_Top_island.png', bbox_inches='tight', dpi=720)
plt.show()


