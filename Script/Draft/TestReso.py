from pyoculus.fields import AnalyticCylindricalBfield
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
from matplotlib import pyplot as plt
import numpy as np
import logging
from matplotlib.path import Path


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

separatrix = {"type": "circular-current-loop", "amplitude": -10, "R": 6, "Z": -5.5}

fig, ax = plt.subplots(1, 1, figsize=(6, 6.6))



perturbedfield = AnalyticCylindricalBfield.with_new_axis(
    R=6, Z=0, sf=SF, shear=SHEAR, perturbations_args=[maxwellboltzmann]
)
perturbedmap = CylindricalBfieldSection.without_axis(perturbedfield, guess=[6., 0])
axiscoords = np.array((perturbedmap.R0, perturbedmap.Z0))


#opoint1 = FixedPoint(perturbedmap)
#opoint1.find(t=1, guess=[6.0, 0.0], method='scipy.root')
#opointcoords1 = opoint1.coords[0]



opointT = FixedPoint(perturbedmap)
#opointT.find(t=3, guess=[7.6, 0.], method='scipy.root')
opointT.find(t=3, guess=[4.7, 1.3], method='scipy.root')
#opointT.find(t=3, guess=[4.7, -1.2], method='scipy.root')
opointTcoords = opointT.coords[0]


#opointR = FixedPoint(perturbedmap)
#opointR.find(t=3, guess=[7.6, 0.], method='scipy.root')
#opointRcoords = opointR.coords[0]

#opointB = FixedPoint(perturbedmap)
#opointB.find(t=3, guess=[4.7, -1.2], method='scipy.root')
#opointBcoords = opointB.coords[0]


xpointL = FixedPoint(perturbedmap)
xpointL.find(t=3, guess=[4.5, 0.0], method='scipy.root')
#xpointL.find(t=3, guess=[6.4, -1.6], method='scipy.root')
xpointLcoords = xpointL.coords[0]

xpointT = FixedPoint(perturbedmap)
xpointT.find(t=3, guess=[6.4, 1.6], method='scipy.root')
#xpointT.find(t=3, guess=[6.4, -1.6], method='scipy.root')
xpointTcoords = xpointT.coords[0]

#xpointB = FixedPoint(perturbedmap)
#xpointB.find(t=3, guess=[6.4, -1.6], method='scipy.root')
#xpointBcoords = xpointB.coords[0]


opointT.plot(ax=ax, marker='o', color="xkcd:green")
#opointB.plot(ax=ax, marker='o', color="xkcd:green")
#opointR.plot(ax=ax, marker='o', color="xkcd:green")
xpointL.plot(ax=ax, marker='x', color="xkcd:green")
xpointT.plot(ax=ax, marker='x', color="xkcd:green")
#xpointB.plot(ax=ax, marker='x', color="xkcd:green")









ResoManifoldDown = Manifold(perturbedmap, xpointL, xpointT)
ResoManifoldDown.compute(
     eps_s=9e-6, eps_u=8e-6, nint_s=8, nint_u=8, neps_s=80, neps_u=80
   )
ResoManifoldDown.plot(ax=ax, markersize=0, lw=1.5)

ResoManifoldDown.find_clinics(first_guess_eps_s=9e-6, first_guess_eps_u=8e-6)
#ResoManifoldDown.plot_clinics(ax=ax)
traj1=ResoManifoldDown.clinics[0].trajectory




ResoManifoldUp = Manifold(perturbedmap, xpointT, xpointL)
ResoManifoldUp.compute(
    eps_s=9e-6, eps_u=8e-6, nint_s=8, nint_u=8, neps_s=80, neps_u=80
)
ResoManifoldUp.plot(ax=ax, markersize=0, lw=1.5)#

ResoManifoldUp.find_clinics(first_guess_eps_s=9e-6, first_guess_eps_u=8e-6)
#ResoManifoldUp.plot_clinics(ax=ax)
traj2=ResoManifoldUp.clinics[0].trajectory

n=100

xL_m1 = np.linspace(xpointLcoords, traj1[-1], n, endpoint=True)
int1 = ResoManifoldDown._AdL_integral_points(xL_m1, is_closed=False)

A=np.linspace(traj1[-1], traj1[-2], n, endpoint=True)
intMid1 = (len(traj1)-1) * ResoManifoldDown._AdL_integral_points(A, is_closed=False)

mn_xR = np.linspace(traj1[0], xpointTcoords, n, endpoint=True)
int2 = ResoManifoldDown._AdL_integral_points(mn_xR, is_closed=False)

xT_m1 = np.linspace(xpointTcoords, traj2[-1], n, endpoint=True)
int3 = ResoManifoldUp._AdL_integral_points(xT_m1, is_closed=False)

B=np.linspace(traj2[-1], traj2[-2], n, endpoint=True)
intMid2 = (len(traj2)-1) * ResoManifoldUp._AdL_integral_points(B, is_closed=False)

mn_xB = np.linspace(traj2[0], xpointLcoords, n, endpoint=True)
int4 = ResoManifoldUp._AdL_integral_points(mn_xB, is_closed=False)

ResoArea= (int1+intMid1+int2+intMid2+int3+int4)


ax.scatter([traj1[0][0], traj1[-1][0],traj1[1][0], traj2[0][0], traj2[-1][0],traj2[1][0]] , [traj1[0][1], traj1[-1][1], traj1[1][1], traj2[0][1], traj2[-1][1], traj2[1][1]], color="red",marker='X', s=80, linewidths=1)

print(ResoArea)
print(int1,intMid1,int2,int3,intMid2,int4)
print(len(traj1), len(traj2))
plt.show()







