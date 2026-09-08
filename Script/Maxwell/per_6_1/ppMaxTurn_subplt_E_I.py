from pyoculus.fields import AnalyticCylindricalBfield
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
from matplotlib import pyplot as plt
import numpy as np
import logging
import logging
from scipy.io import loadmat
from matplotlib.patches import FancyArrowPatch

def draw_curved_arrow(ax, p0, p1, rad=0.2, color='black', lw=1.5, arrowstyle='->', mutation_scale=15, zorder=300):
    """
    Draw a curved arrow from p0 to p1 (p0/p1 are (x,y) in data coords).
    rad: curvature (positive -> curve one side, negative -> other side).
    """
    arrow = FancyArrowPatch(posA=(p0[0], p0[1]), posB=(p1[0], p1[1]),
                            connectionstyle=f"arc3,rad={rad}",
                            arrowstyle=arrowstyle,
                            mutation_scale=mutation_scale,
                            color=color, linewidth=lw, zorder=zorder)
    ax.add_patch(arrow)
    return arrow

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

#fig, (ax2, ax) = plt.subplots(1, 2, figsize=(8, 7), width_ratios=[1,1], constrained_layout=False,sharex=True,sharey=True)


fig, ax = plt.subplots(1, 1, figsize=(6, 6.5))

fig2, ax2 = plt.subplots(1, 1, figsize=(6, 6.5))

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
# xpoint.plot(ax=ax, s=100,marker='x', color="xkcd:crimson",zorder=20,label='X-point')
# xpoint.plot(ax=ax2, s=100,marker='x', color="xkcd:crimson",zorder=20,label='X-point')


opoint = FixedPoint(perturbedmap)
opoint.find(t=1, guess=[6., 0.])
opointcoords = opoint.coords[0]
# opoint.plot(ax=ax, s=70,marker='o', color="xkcd:crimson",zorder=20,label='O-point')
# opoint.plot(ax=ax2, s=70,marker='o', color="xkcd:crimson",zorder=20,label='O-point')


xpoint.plot(ax=ax, s=100,marker='x', color="xkcd:darkgreen",zorder=20,label='X-point')
xpoint.plot(ax=ax2, s=100,marker='x', color="xkcd:darkgreen",zorder=20,label='X-point')

opoint.plot(ax=ax, s=70,marker='o', color="xkcd:darkgreen",zorder=20,label='O-point')
opoint.plot(ax=ax2, s=70,marker='o', color="xkcd:darkgreen",zorder=20,label='O-point')


separatrixcoords = np.array([6.0, -5.4])


# pplot = PoincarePlot.with_segments(
#    perturbedmap,
#    [opointcoords, xpointcoords, xpointcoords, separatrixcoords, xpointcoords, xpointcoords + np.array([0.1, 0])],
#    neps=[POINCARE_TRAJ, 5, 5],
#    connected=False)

# pplot.compute(POINCARE_ITS)
# pplot.plot(ax=ax, color="xkcd:dark grey", s=1.4, linewidths=0)
# np.save(f"{repository_path}hits_max_6_1.npy", pplot._hits)





Hits=np.load(f"{repository_path}hits_max_6_1.npy")
ax.scatter(Hits[:,:, 0], Hits[:,:, 1], color="xkcd:dark grey", s=0.5, linewidths=0,zorder=10)
ax2.scatter(Hits[:,:, 0], Hits[:,:, 1], color="xkcd:dark grey", s=0.5, linewidths=0,zorder=10)


ax.set_xlabel(r"$R[m]$")
ax2.set_ylabel(r"$Z[m]$")
ax.set_aspect('equal')

ax2.set_xlabel(r"$R[m]$")
#ax2.set_ylabel(r"$Z[m]$")
ax2.set_aspect('equal')

# perturbedmanifold = Manifold(perturbedmap, xpoint, xpoint)
# perturbedmanifold.compute(
#     eps_s=9e-6, eps_u=8e-6, nint_s=8, nint_u=9, neps_s=80, neps_u=2000
# )
# perturbedmanifold.find_clinics(first_guess_eps_s=9e-6, first_guess_eps_u=8e-6)
# perturbedmanifold.save(f"{repository_path}mf_nu_9.pkl")



mf_E  = Manifold.load(f"{repository_path}mf_nu_9.pkl")
mf_I  = Manifold.load(f"{repository_path}mf_ns_9.pkl")


fig, ax=mf_E.plot_clinics(ax=ax,s=20,label='homoclinics')
fig2, ax2=mf_I.plot_clinics(ax=ax2,s=20,label='homoclinics')

mf_E.plot(ax=ax, markersize=0, lw=1.)
mf_I.plot(ax=ax2, markersize=0, lw=1.)

cmap = plt.get_cmap('plasma') 
cmap2 = plt.get_cmap('gist_earth')
n = 7

vals = np.linspace(0.0, 0.6, n) 
 
            # linear sampling from 0..1
color1 = cmap(vals)  


vals2 = np.linspace(0.5, 0.0, n)               # linear sampling from 0..1
vals2 = np.linspace(0.0, 0.5, n)   


color2= cmap2(vals2)


#color1=['darkgreen','green','mediumseagreen','limegreen','lightgreen','palegreen']
#color2=['darkblue','blue','royalblue','darkcyan','mediumseagreen','green','darkgreen']
#color1=['darkred','firebrick','red','mediumvioletred','deeppink','hotpink','magenta']




# perturbedmanifold.plot_filled_lobe(ax=ax, lobe_number=3,which_section=2,alpha=1,color='darkgreen')

text1=[r'$f^{-2}(E)$', r'$f^{-1}(E)$',r'$E$', r'$f(E)$', r'$f^2(E)$',r'$f^3(E)$',r'$f^4(E)$']

text2=[r'$f^{-4}(I)$',r'$f^{-3}(I)$', r'$f^{-2}(I)$', r'$f^{-1}(I)$',r'$I$', r'$f(I)$',r'$f^2(I)$']

dx=[0.06,0.15,-0.05,-0.2,-0.08,-0.11,-0.2]
dy=[+0.01,0.06,0.15,0,-0.1,-0.06,-0.15]

# dx2=[0.08,0.09,0.07,-0.1,-0.07,-0.04,-0.1]
# dy2=[-0.12,0,0.06,0.04,-0.04,-0.08,-0.04]

dx2=[0.08,0.09,-0.1,-0.1,-0.07,-0.04,-0.1] #3eme 0.07
dy2=[-0.12,0,0.06,0.04,-0.04,-0.08,-0.04]


a=mf_E.clinics[0].trajectory

b=mf_I.clinics[1].trajectory

for i in range(len(text1)):

    mf_E.plot_filled_lobe(ax=ax, lobe_number=i+3,neps=10 if i!=6 else 500,alpha=0.7,color=color1[i])
    mf_I.plot_filled_lobe(ax=ax2, lobe_number=i+2,neps=10 if i!=0 else 500,which_section=2,alpha=0.7,color=color2[i])

    _place_coord_label(ax, a[i+3,0], a[i+3,1], fmt=text1[i], dx_frac=dx[i], dy_frac=dy[i],color=color1[i])#color="darkred"
    _place_coord_label(ax2, b[i+2,0], b[i+2,1], fmt=text2[i], dx_frac=dx2[i], dy_frac=dy2[i],color=color2[i])#color="darkblue"

_place_coord_label(ax, 5.5, 1.9, fmt=r'm', dx_frac=-0.05, dy_frac=-0.03, color="xkcd:black",fontsize=17)    
_place_coord_label(ax2, 5.5, 1.9, fmt=r'm', dx_frac=-0.05, dy_frac=-0.03, color="xkcd:black",fontsize=17)

#_place_coord_label(ax, xpointcoords[0], xpointcoords[1], fmt=r'h', dx_frac=-0.03, dy_frac=-0.06, color="xkcd:black",fontsize=17)


draw_curved_arrow(ax, (7.9,1.25), (4.75,1.3), rad=-0.25, color='black', lw=1.6)

ax.text(6.3, 1.13, r'$f$',
            color='black', fontsize=17, fontweight="normal", 
            ha="left", va="center",
            zorder=120)


draw_curved_arrow(ax2, (8.6,-0.85), (6.3,1.95), rad=-0.25, color='black', lw=1.6)

ax2.text(7.3, 0.3, r'$f$',
            color='black', fontsize=17, fontweight="normal", 
            ha="left", va="center",
            zorder=120)

# ax.text( 0.05, 0.96, 
#    '(b)', transform=ax.transAxes,fontsize=12,fontweight='bold',ha='left',va='top',color='black', zorder=200)

# ax2.text( 0.05, 0.96, 
#    '(a)', transform=ax2.transAxes,fontsize=12,fontweight='bold',ha='left',va='top',color='black', zorder=200)

#perturbedmanifold.compute_turnstile_areas()

#ax2.legend(loc='lower left', bbox_to_anchor=(0.0, 0.0),fontsize=7)
#ax2.legend(loc='lower left', bbox_to_anchor=(-0.02, -0.01), fontsize=7, frameon=False)
#ax2.legend(loc='lower right', fontsize=8)

#texte=f"Turnstile flux: $\phi = ${perturbedmanifold.turnstile_areas[0]:.3e}"

#ax.text(2, 3, texte, fontsize=12, color='blue')

ax.set_xlim(3.2, 9.4)
ax.set_ylim(-6.2, 2.8) 

ax2.set_xlim(3.2, 9.4)
ax2.set_ylim(-6.2, 2.8)

#plt.savefig(f"{repository_path}figure_turnstile_both.png",  bbox_inches='tight', pad_inches=0.0, dpi=720)
fig.savefig(f"{repository_path}ORAL_turnstile_exit_set_v3.png", bbox_inches='tight', pad_inches=0.0, dpi=720)
fig2.savefig(f"{repository_path}ORAL_turnstile_entry_set_v3.png", bbox_inches='tight', pad_inches=0.0, dpi=720)

#print(perturbedmanifold.turnstile_areas)

plt.show()

