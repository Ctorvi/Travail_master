from pyoculus.fields import AnalyticCylindricalBfield
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
from matplotlib import pyplot as plt
import numpy as np
import logging
import logging
from scipy.io import loadmat

plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.size": 10,
    }
)
logging.basicConfig(level=logging.DEBUG)

repository_path = './script/maxwell/per_6_1/'

POINCARE_TRAJ = 15 
POINCARE_ITS = 500

SHEAR = 1.2
SF = 1.16

#fig, (ax, ax2) = plt.subplots(1, 2, figsize=(8, 7), width_ratios=[1,1], constrained_layout=False,sharex=True,sharey=True)

fig, ax = plt.subplots(1, 1, figsize=(6, 6.5))

fig2, ax2 = plt.subplots(1, 1, figsize=(6, 6.5))

separatrix = {"type": "circular-current-loop", "amplitude": -10, "R": 6, "Z": -5.5}

perturbedfield = AnalyticCylindricalBfield.with_new_axis(
    R=6, Z=0, sf=SF, shear=SHEAR, perturbations_args=[separatrix]
)
perturbedmap = CylindricalBfieldSection.without_axis(perturbedfield, guess=[6.0, 0])
axiscoords = np.array((perturbedmap.R0, perturbedmap.Z0))

xpointguess = [6.205, -4.50]
xpoint = FixedPoint(perturbedmap)
xpoint.find(t=1, guess=xpointguess)
xpoint.m = 1
xpointcoords = xpoint.coords[0]
xpoint.plot(ax=ax, s=100,marker='x', color="xkcd:darkgreen",zorder=20,label='X-point')


opoint = FixedPoint(perturbedmap)
opoint.find(t=1, guess=[6., 0.])
opointcoords = opoint.coords[0]
opoint.plot(ax=ax, s=70,marker='o', color="xkcd:darkgreen",zorder=20,label='O-point')

separatrixcoords = np.array([6.0, -5.4])

# pplot = PoincarePlot.with_segments(
#    perturbedmap,
#    [opointcoords, xpointcoords, xpointcoords, separatrixcoords, xpointcoords, xpointcoords + np.array([0.1, 0])],
#    neps=[POINCARE_TRAJ, 5, 5],
#    connected=False)

# pplot.compute(POINCARE_ITS)
# pplot.plot(ax=ax, color="xkcd:dark grey", s=1.4, linewidths=0)
# np.save(f"{repository_path}hits_max_6_1_NP.npy", pplot._hits)


Hits_NP=np.load(f"{repository_path}hits_max_6_1_NP.npy")
ax.scatter(Hits_NP[:,:, 0], Hits_NP[:,:, 1], color="xkcd:dark grey", s=0.5, linewidths=0,zorder=10)


Hits=np.load(f"{repository_path}hits_max_6_1.npy")
ax2.scatter(Hits[:,:, 0], Hits[:,:, 1], color="xkcd:dark grey", s=0.5, linewidths=0,zorder=10)


ax.set_xlabel(r"$R[m]$")
#ax2.set_yticks([])
ax.set_aspect('equal')
ax2.set_xlabel(r"$R[m]$")
ax2.set_aspect('equal')
ax.set_ylabel(r"$Z[m]$")
ax2.set_ylabel(r"$Z[m]$")

# perturbedmanifold = Manifold(perturbedmap, xpoint, xpoint)
# perturbedmanifold.compute(
#     eps_s=9e-6, eps_u=8e-6, nint_s=8, nint_u=8, neps_s=80, neps_u=80
# )
# perturbedmanifold.save(f"{repository_path}mf_NP.pkl")


maxwellboltzmann = {
    "m": 6,
    "n": -1,
    "d": 2,
    "type": "maxwell-boltzmann",
    "amplitude": 0.32,
}


perturbedfield = AnalyticCylindricalBfield.with_new_axis(
    R=6, Z=0, sf=SF, shear=SHEAR, perturbations_args=[separatrix, maxwellboltzmann]
)
perturbedmap = CylindricalBfieldSection.without_axis(perturbedfield, guess=[6.0, 0])
axiscoords = np.array((perturbedmap.R0, perturbedmap.Z0))

xpointguess = [6.205, -4.50]
xpoint = FixedPoint(perturbedmap)
xpoint.find(t=1, guess=xpointguess)
xpoint.m = 1
xpointcoords = xpoint.coords[0]
xpoint.plot(ax=ax2, s=100,marker='x', color="xkcd:darkgreen",zorder=20,label='X-point')


opoint = FixedPoint(perturbedmap)
opoint.find(t=1, guess=[6., 0.])
opointcoords = opoint.coords[0]
opoint.plot(ax=ax2, s=70,marker='o', color="xkcd:darkgreen",zorder=20,label='O-point')



NP_manifold  = Manifold.load(f"{repository_path}mf_NP.pkl")
#NP_manifold.plot(ax=ax, markersize=0, lw=1.2)

P_manifold  = Manifold.load(f"{repository_path}mf.pkl")
P_manifold.plot_clinics(ax=ax2, s=20,label='Homoclinics')
P_manifold.plot(ax=ax2, markersize=0, lw=1.2)

# ax.text( 0.05, 0.93, 
#    '(a)', transform=ax.transAxes,fontsize=12,fontweight='bold',ha='left',va='top',color='black', zorder=200)

# ax2.text( 0.05, 0.93, 
#    '(b)', transform=ax2.transAxes,fontsize=12,fontweight='bold',ha='left',va='top',color='black', zorder=200)

# ax2.legend(loc='lower left', bbox_to_anchor=(-0.02, -0.01), ncol=1, fontsize=7)

#texte=f"Turnstile flux: $\phi = ${perturbedmanifold.turnstile_areas[0]:.3e}"

#ax.text(2, 3, texte, fontsize=12, color='blue')

#plt.savefig(f"{repository_path}figure_chaos.png", dpi=720,bbox_inches="tight", pad_inches=0, facecolor=fig.get_facecolor())
#print(perturbedmanifold.turnstile_areas)

ax.set_xlim(3.2, 9.4)
ax.set_ylim(-6.2, 2.8)
ax2.set_xlim(3.2, 9.4)
ax2.set_ylim(-6.2, 2.8)


#fig.savefig(f"{repository_path}oral_max6_1_axi.png", bbox_inches='tight', pad_inches=0.0, dpi=720)
fig2.savefig(f"{repository_path}oral_max6_1_nonaxi_w_MF.png", bbox_inches='tight', pad_inches=0.0, dpi=720)


plt.show()

