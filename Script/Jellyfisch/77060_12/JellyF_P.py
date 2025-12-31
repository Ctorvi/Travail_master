import os, sys
proj_root = os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..','..'))  # remonte jusqu'à "Travail master"
if proj_root not in sys.path:
    sys.path.insert(0, proj_root)
from script.function.field_utils import plot_tomographic as tomo
from script.function.field_utils import plot_Br_Bz
from script.function.field_utils import plot_LIUQE as LIUQE

from pyoculus.fields import AxisymmetricCylindricalGridField
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
import numpy as np
from matplotlib import pyplot as plt
import logging
from scipy.io import loadmat
from matplotlib.tri import Triangulation



patch_mat_file = './script/jellyfisch/77060_12/JF_patch_77060_120_C3.mat'

repository_path = './script/jellyfisch/77060_12/'

pert_mat_file = 'JF_77060_120_BO78_OS.mat'

mf_list = ['mf_1T', 'mf_1B', 'mf_3L', 'mf_3R']

file_manifolds = 'manifolds_P/OS_BO78/'

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(5, 8))
#fig, ax = plt.subplots(1, 1, figsize=(6, 6.6))

plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.size": 9,
    }
)

logging.basicConfig(level=logging.DEBUG)

JFField = AxisymmetricCylindricalGridField.from_matlab_file(f"{repository_path}{pert_mat_file}", with_perturbation=True)

section = CylindricalBfieldSection(JFField,phi0=1.9,R0=0.88, Z0=0)

# top_o = FixedPoint(section)
# top_o.find(1, [0.9,0.18], method='scipy.root')
# top_o_coord = top_o.coords[0]

# x_point1= FixedPoint(section)
# x_point1.find(1, [0.8,-0.27], method='scipy.root')
# x_point1_coord = x_point1.coords[0]

# x_point2= FixedPoint(section)
# x_point2.find(1, [1.05,-0.58], method='scipy.root')
# x_point2_coord = x_point2.coords[0]

###########################  tomographic reconstruction  #########################



###################.  poincaré plot #####################

# pplot = PoincarePlot.with_linspace(section, top_o_coord, x_point1_coord,40)
# pplot.compute(400)
# np.save(f"{repository_path}JF_77060_pp_hits.npy", pplot._hits)


# Hits=np.load(f"{repository_path}JF_77060_pp_hits.npy")
# ax.scatter(Hits[:,:, 0], Hits[:,:, 1], color="xkcd:dark grey", s=1., linewidths=0)


# tpc= tomo(f"{patch_mat_file}",ax,emi_vmin=0, emi_vmax=6e20)

###########################. fixed point  #########################

#####top fp top mf#########


# manifold_1T = Manifold(section, x_point1, x_point1,-x_point1_coord+top_o_coord, -x_point1_coord+top_o_coord)
# manifold_1T.compute( eps_s=9e-6, eps_u=8e-6, nint_s=10, nint_u=14, neps_s=80, neps_u=240) #8 14 240
# manifold_1T.save(f"{repository_path}{file_manifolds}/mf_1T.pkl")


# #####top fp bottom mf#########

# manifold_1B = Manifold(section, x_point1, x_point1)
# manifold_1B.compute(eps_s=9e-6, eps_u=8e-6, nint_s=10, nint_u=12, neps_s=80, neps_u=80)
# manifold_1B.save(f"{repository_path}{file_manifolds}/mf_1B.pkl")

# ##### right fp left mf

# manifold_3L = Manifold(section, x_point2, x_point2,-x_point2_coord+x_point1_coord, -x_point2_coord+x_point1_coord)
# manifold_3L.compute(eps_s=9e-6, eps_u=8e-6, nint_s=40, nint_u=40, neps_s=80, neps_u=80)
# manifold_3L.save(f"{repository_path}{file_manifolds}/mf_3L.pkl")

# # # ##### right fp right mf

# manifold_3R = Manifold(section, x_point2, x_point2,x_point2_coord-top_o_coord, x_point2_coord-x_point1_coord)
# manifold_3R.compute(eps_s=9e-6, eps_u=8e-6, nint_s=40, nint_u=40, neps_s=80, neps_u=80)
# manifold_3R.save(f"{repository_path}{file_manifolds}/mf_3R.pkl")

##### loading and plotting  ###



# manifs = {name: {name} for name in mf_list}

# for i in mf_list:
#      manifs[i] = Manifold.load(f"{repository_path}{file_manifolds}{i}.pkl")
#      manifs[i].plot(stepsize_limit=0.05, ax=ax, markersize=0, lw=0.5,colors=["rosybrown", "xkcd:red"],labels=['stable MF','unstable MF'] if i=='mf_1T' else [None,None])


# top_o.plot(ax=ax, marker='o', color="xkcd:white")
# x_point1.plot(ax=ax, marker='x', color="xkcd:white")
# x_point2.plot(ax=ax, marker='x', color="xkcd:white")

#####figure

# ax.set_xlim(0.62, 1.15)
# ax.set_ylim(-0.75, 0.75)
# ax.set_xlabel(r"$R[m]$")
# ax.set_ylabel(r"$Z[m]$")
# ax.set_aspect('equal') 
# ax.set_title('(77060 / 1.2s / 1.9 rad)')
# #plt.savefig('./Script/Jellyfisch/77060_12/figures/JF_77062_mf_C3.png', bbox_inches='tight', dpi=720)
# plt.show()

#####zoomed

# ax.set_xlim(0.7, 1.0)
# ax.set_ylim(-0.6, -0.1)
# ax.set_xlabel(r"$R[m]$")
# ax.set_ylabel(r"$Z[m]$")
# ax.set_aspect('equal') 
# ax.set_title('(77060 / 1.2s / 1.9 rad)')
# plt.savefig('./Script/Jellyfisch/77060_12/figures/JF_77062_mf_C3_zoom.png', bbox_inches='tight', dpi=720)
# plt.show()




