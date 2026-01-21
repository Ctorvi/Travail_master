import os, sys
proj_root = os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..','..'))  # remonte jusqu'à "Travail master"
if proj_root not in sys.path:
    sys.path.insert(0, proj_root)
from script.function.field_utils import plot_tomographic as tomo
from script.function.field_utils import plot_LIUQE as LIUQE
from script.function.field_utils import plot_vessel as vessel

from pyoculus.fields import AxisymmetricCylindricalGridField
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
import numpy as np
from matplotlib import pyplot as plt
import logging
from scipy.io import loadmat
from matplotlib.tri import Triangulation

#fig, ax = plt.subplots(1, 1, figsize=(5, 8))
fig, ax = plt.subplots(1, 1, figsize=(6, 6.6))

plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.size": 9,
    }
)
logging.basicConfig(level=logging.DEBUG)

repository_path = './script/snowflake/backoff/'

pert_mat_file = 'SF_70620_09_BO78_OS.mat'


JFField = AxisymmetricCylindricalGridField.from_matlab_file(f"{repository_path}{pert_mat_file}", with_perturbation=True)

section = CylindricalBfieldSection(JFField,phi0=1.9,R0=0.88, Z0=0)

top_o = FixedPoint(section)
top_o.find(1, [0.9,0.0], method='scipy.root')
top_o_coord = top_o.coords[0]
top_o.plot(ax=ax, marker='o', color="xkcd:crimson")

x_point1= FixedPoint(section)
x_point1.find(1, [0.75,-0.4], method='scipy.root')
x_point1_coord = x_point1.coords[0]
x_point1.plot(ax=ax, marker='x', color="xkcd:crimson")

x_point2= FixedPoint(section)
x_point2.find(1, [0.8,-0.6], method='scipy.root')
x_point2_coord = x_point2.coords[0]
x_point2.plot(ax=ax, marker='x', color="xkcd:crimson")


########### poncaré plot #########

# pplot = PoincarePlot.with_segments(section, [top_o_coord, x_point1_coord, x_point1_coord, x_point2_coord],[30, 10],connected=False)
# pplot.compute(400)
# np.save(f"{repository_path}Snowflake_Pert.npy", pplot._hits)

###########################.tomo #########################


# ####top fp top mf#########


# manifold_1T = Manifold(section, x_point1, x_point1,-x_point1_coord+top_o_coord, -x_point1_coord+top_o_coord)
# manifold_1T.compute(
#    eps_s=9e-6, eps_u=8e-6, nint_s=10, nint_u=8, neps_s=80, neps_u=80) #8 14 240
# manifold_1T.save(f"{repository_path}/manifolds_P/mf_1T.pkl")

# #####top fp bottom mf#########


# manifold_1B = Manifold(section, x_point1, x_point1,x_point1_coord-top_o_coord, x_point1_coord-top_o_coord)
# manifold_1B.compute(
#       eps_s=9e-6, eps_u=8e-6, nint_s=10, nint_u=12, neps_s=80, neps_u=80)
# manifold_1B.save(f"{repository_path}/manifolds_P/mf_1B.pkl")


# # ##### bottom fp top mf

# manifold_2T = Manifold(section, x_point2, x_point2,-x_point2_coord+x_point1_coord, -x_point2_coord+x_point1_coord)
# manifold_2T.compute(
#      eps_s=9e-6, eps_u=8e-6, nint_s=8, nint_u=15, neps_s=80, neps_u=80)
# manifold_2T.save(f"{repository_path}/manifolds_P/mf_2T.pkl")

# # ##### bottom fp bottom mf

# manifold_2B = Manifold(section, x_point2, x_point2,x_point2_coord-x_point1_coord, x_point2_coord-x_point1_coord)
# manifold_2B.compute(
#      eps_s=9e-6, eps_u=8e-6, nint_s=40, nint_u=40, neps_s=80, neps_u=80)
# manifold_2B.save(f"{repository_path}/manifolds_P/mf_2B.pkl")

# ### loading and plotting  ###

manifold_1T  = Manifold.load(f"{repository_path}/manifolds_P/mf_1T.pkl")
manifold_1B  = Manifold.load(f"{repository_path}/manifolds_P/mf_1B.pkl")
manifold_2T  = Manifold.load(f"{repository_path}/manifolds_P/mf_2T.pkl")
manifold_2B  = Manifold.load(f"{repository_path}/manifolds_P/mf_2B.pkl")

manifold_1T.plot(stepsize_limit=0.1, ax=ax, markersize=0, lw=0.7,colors=["xkcd:royal blue", "xkcd:red"],labels=["Stable MF 1T","Unstable MF 1T"])
manifold_1B.plot(stepsize_limit=0.3, ax=ax, markersize=0, lw=0.7,colors=["xkcd:royal blue", "xkcd:red"],labels=[None,None])
manifold_2T.plot(stepsize_limit=0.1, ax=ax, markersize=0, lw=0.7,colors=["xkcd:royal blue", "xkcd:red"],labels=[None,None])
manifold_2B.plot(stepsize_limit=0.3, ax=ax, markersize=0, lw=0.7,colors=["xkcd:royal blue", "xkcd:red"],labels=[None,None])


top_o.plot(ax=ax, marker='o',s=40, color="xkcd:dark blue",label='O-points', zorder=10)
x_point1.plot(ax=ax, marker='x', s=50, color="xkcd:dark blue",label='X-points', zorder=10)
x_point2.plot(ax=ax, marker='x', s=50, color="xkcd:dark blue",label=None, zorder=10)

Hits=np.load(f"{repository_path}Snowflake_Pert.npy")
ax.scatter(Hits[:,:, 0], Hits[:,:, 1], color="xkcd:dark grey", s=0.5, linewidths=0)
 
#vessel(ax)

###figure t

#ax.legend(handles, labels,loc='lower right', bbox_to_anchor=(1.05, 0.0), ncol=1, fontsize=9)
ax.legend(loc='lower right', bbox_to_anchor=(1., 0.76), ncol=1, fontsize=9)

ax.set_xlim(0.62, 1.15)
ax.set_ylim(-0.75, 0.75)
ax.set_xlabel(r"$R[m]$")
ax.set_ylabel(r"$Z[m]$")
ax.set_aspect('equal') 
ax.set_title('Snowflake 70620 at 0.9 s')
plt.savefig(f"{repository_path}/figures/Snowflake_Pert.png", bbox_inches="tight", pad_inches=0, facecolor=fig.get_facecolor())
plt.show()

### zoomed

# ax.set_xlim(0.7, 1.0)
# ax.set_ylim(-0.6, -0.1)
# ax.set_xlabel(r"$R[m]$")
# ax.set_ylabel(r"$Z[m]$")
# ax.set_aspect('equal') 
# ax.set_title('70620 at 0.9 s, angle=1.9 rad')
# plt.savefig(f"{repository_path}/figures/Snowflake_Pert_mf_patch_Zoom.png", bbox_inches='tight', dpi=720)
# plt.show()



