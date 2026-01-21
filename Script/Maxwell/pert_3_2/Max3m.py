from pyoculus.fields import AnalyticCylindricalBfield
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
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

repository_path = './script/maxwell/pert_3_2/'

POINCARE_ITS = 400
SHEAR =  0.4
SF = 0.8875

maxwellboltzmann = {
    "m": 3,
    "n": 2,
    "d": 1.75/np.sqrt(2),
    "type": "maxwell-boltzmann",
    "amplitude": 0.1,
}

separatrix = {"type": "circular-current-loop", "amplitude": -10, "R": 6, "Z": -5.5}

fig, ax = plt.subplots(1, 1, figsize=(6, 6.6))



perturbedfield = AnalyticCylindricalBfield.with_new_axis(
    R=6, Z=0, sf=SF, shear=SHEAR, perturbations_args=[maxwellboltzmann]
)
perturbedmap = CylindricalBfieldSection.without_axis(perturbedfield, guess=[6., 0])
axiscoords = np.array((perturbedmap.R0, perturbedmap.Z0))

pplot = PoincarePlot.with_linspace(
    perturbedmap,
    [6.0, 0],
    [8.0, 0],
    60
   )




# pplot.compute(POINCARE_ITS)
# # pplot.plot(ax=ax, color="xkcd:dark grey", s=1.4, linewidths=0)

# np.save(f"{repository_path}hits_max_3_2_P_v2.npy", pplot._hits)


Hits=np.load(f"{repository_path}hits_max_3_2_P_v2.npy")
ax.scatter(Hits[:,:, 0], Hits[:,:, 1], color="xkcd:dark grey", s=1, linewidths=0,zorder=10)




opoint1 = FixedPoint(perturbedmap)
opoint1.find(t=1, guess=[6.0, 0], method='scipy.root')
opointcoords1 = opoint1.coords[0]


opointT = FixedPoint(perturbedmap)
opointT.find(t=3, guess=[4.3, 1.3], method='scipy.root')
opointTcoords = opointT.coords[0]

opointR= FixedPoint(perturbedmap)
opointR.find(t=3, guess=[7.7, 0.], method='scipy.root')
opointRcoords = opointR.coords[0]

opointL= FixedPoint(perturbedmap)
opointL.find(t=3, guess=[4.5, -1.3], method='scipy.root')
opointLcoords = opointL.coords[0]

xpointB = FixedPoint(perturbedmap)
xpointB.find(t=3, guess=[6.4, -1.6], method='scipy.root')
xpointBcoords = xpointB.coords[0]

xpointT = FixedPoint(perturbedmap)
xpointT.find(t=3, guess=[6.4, 1.6], method='scipy.root')
xpointTcoords = xpointT.coords[0]

xpointL = FixedPoint(perturbedmap)
xpointL.find(t=3, guess=[4.5, 0.], method='scipy.root')
xpointLcoords = xpointL.coords[0]



opoint1.plot(ax=ax, s=50,marker='o', color="xkcd:dark blue", zorder=20, label='O-point')
opointT.plot(ax=ax, s=50, marker='o', color="xkcd:dark blue", zorder=20, label=None)
#opointR.plot(ax=ax, s=50, marker='x', color="xkcd:dark blue", zorder=20)
xpointB.plot(ax=ax, s=50, marker='x', color="xkcd:dark blue", zorder=20, label='X-points')
#opointL.plot(ax=ax, s=50, marker='x', color="xkcd:dark blue", zorder=20)
#xpointT.plot(ax=ax, s=50, marker='x', color="xkcd:dark blue", zorder=20)
#xpointL.plot(ax=ax, s=50, marker='x', color="xkcd:dark blue", zorder=20)


ax.legend(loc='upper right', bbox_to_anchor=(1., 0.99), ncol=1, fontsize=9)
ax.set_xlabel(r"$R[m]$")
ax.set_ylabel(r"$Z[m]$")
plt.savefig('./script/maxwell/pert_3_2/ppMax3m.png', bbox_inches='tight', dpi=720)
plt.show()