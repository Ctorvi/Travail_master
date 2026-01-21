from pyoculus.fields import AnalyticCylindricalBfield
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
from matplotlib import pyplot as plt
import numpy as np
import logging
import logging
from scipy.io import loadmat


def _place_coord_label(ax, x, y, fmt="({:.3f}, {:.3f})", dx_frac=0.02, dy_frac=0.0, color="black", fontsize=10):
    """
    Place a label near (x,y) in data coordinates.
    - dx_frac : horizontal offset as fraction of x-axis width (default 0.02)
    - dy_frac : vertical offset as fraction of y-axis height (default 0.0)
    Uses matplotlib mathtext (r'$\ldots$') if fmt contains LaTeX-like strings.
    """
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    dx = dx_frac * (xlim[1] - xlim[0])
    dy = dy_frac * (ylim[1] - ylim[0])
    ax.text(x + dx, y + dy, fmt,
            color=color, fontsize=fontsize, fontweight="normal",
            ha="left", va="center",
            zorder=120)

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

fig, ax = plt.subplots(1, 1, figsize=(6, 6.5))

separatrix = {"type": "circular-current-loop", "amplitude": -10, "R": 6, "Z": -5.5}

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
xpoint.plot(ax=ax, s=100,marker='x', color="xkcd:crimson",zorder=20,label='X-point')


opoint = FixedPoint(perturbedmap)
opoint.find(t=1, guess=[6., 0.])
opointcoords = opoint.coords[0]
opoint.plot(ax=ax, s=70,marker='o', color="xkcd:crimson",zorder=20,label='O-point')

separatrixcoords = np.array([6.0, -5.4])


pplot = PoincarePlot.with_segments(
   perturbedmap,
   [opointcoords, xpointcoords, xpointcoords, separatrixcoords, xpointcoords, xpointcoords + np.array([0.1, 0])],
   neps=[POINCARE_TRAJ, 5, 5],
   connected=False)

# pplot.compute(POINCARE_ITS)
# pplot.plot(ax=ax, color="xkcd:dark grey", s=1.4, linewidths=0)
# np.save(f"{repository_path}hits_max_6_1.npy", pplot._hits)





Hits=np.load(f"{repository_path}hits_max_6_1.npy")
ax.scatter(Hits[:,:, 0], Hits[:,:, 1], color="xkcd:dark grey", s=0.5, linewidths=0,zorder=10)




ax.set_xlabel(r"$R[m]$")
ax.set_ylabel(r"$Z[m]$")
ax.set_aspect('equal')




# perturbedmanifold = Manifold(perturbedmap, xpoint, xpoint)
# perturbedmanifold.compute(
#     eps_s=9e-6, eps_u=8e-6, nint_s=8, nint_u=8, neps_s=80, neps_u=80
# )
# perturbedmanifold.find_clinics(first_guess_eps_s=9e-6, first_guess_eps_u=8e-6)
# perturbedmanifold.save(f"{repository_path}mf.pkl")


perturbedmanifold  = Manifold.load(f"{repository_path}mf.pkl")
perturbedmanifold.plot(ax=ax, markersize=0, lw=1.2)

fig, ax=perturbedmanifold.plot_clinics(ax=ax,label='homoclinics')


#color1=['darkgreen','green','mediumseagreen','limegreen','lightgreen','palegreen']
color1=['darkblue','blue','deepskyblue','lightskyblue','lightblue','paleturquoise']
color2=['darkred','red','orangered','salmon','lightsalmon','peachpuff']

color2.reverse()

# perturbedmanifold.plot_filled_lobe(ax=ax, lobe_number=3,which_section=2,alpha=1,color='darkgreen')


text1=[r'$f^{-2}(E)$', r'$f^{-1}(E)$',r'$E$', r'$f^1(E)$', r'$f^2(E)$',r'$f^3(E)$']
text2=[r'$f^{-3}(I)$', r'$f^{-2}(I)$', r'$f^{-1}(I)$',r'$I$', r'$f(I)$',r'$f^2(I)$']


dx=[0.04,0.15,-0.05,-0.2,-0.08,-0.11]
dy=[-0.02,0.06,0.15,0,-0.1,-0.06]

dx2=[0.09,0.07,-0.1,-0.07,-0.04,-0.08]
dy2=[0,0.06,0.04,-0.04,-0.08,-0.05]

a=perturbedmanifold.clinics[0].trajectory
b=perturbedmanifold.clinics[1].trajectory

for i in range(len(text1)):
    perturbedmanifold.plot_filled_lobe(ax=ax, lobe_number=i+3,alpha=0.4,color=color1[i])
    perturbedmanifold.plot_filled_lobe(ax=ax, lobe_number=i+3,which_section=2,alpha=0.4,color=color2[i])
    
    _place_coord_label(ax, a[i+3,0], a[i+3,1], fmt=text1[i], dx_frac=dx[i], dy_frac=dy[i],color='darkblue')
    _place_coord_label(ax, b[i+3,0], b[i+3,1], fmt=text2[i], dx_frac=dx2[i], dy_frac=dy2[i],color='darkred')

_place_coord_label(ax, 5.5, 1.9, fmt=r'm', dx_frac=-0.05, dy_frac=-0.03, color="xkcd:black",fontsize=17)    
#_place_coord_label(ax, xpointcoords[0], xpointcoords[1], fmt=r'h', dx_frac=-0.03, dy_frac=-0.06, color="xkcd:black",fontsize=17)
#perturbedmanifold.compute_turnstile_areas()


ax.legend(loc='lower left', fontsize=8)

#texte=f"Turnstile flux: $\phi = ${perturbedmanifold.turnstile_areas[0]:.3e}"

#ax.text(2, 3, texte, fontsize=12, color='blue')

#plt.savefig(f"{repository_path}figure_turnstile_v2.png", bbox_inches='tight', dpi=720)
#print(perturbedmanifold.turnstile_areas)
plt.show()

