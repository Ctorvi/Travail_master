import os, sys
proj_root = os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..'))  # remonte jusqu'à "Travail master"
if proj_root not in sys.path:
    sys.path.insert(0, proj_root)
from script.function.field_utils import plot_tomographic as tomo

from pyoculus.fields import AxisymmetricCylindricalGridField
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
import numpy as np
from matplotlib import pyplot as plt
import logging
from scipy.io import loadmat
from matplotlib.tri import Triangulation

fig, ax = plt.subplots(1, 1, figsize=(6, 6.6))
#fig,ax = plt.subplots(1,1, figsize=(5,8))
plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.size": 9,
    }
)
logging.basicConfig(level=logging.DEBUG)

repository_path = './script/standard/'

pert_mat_file = 'ST_80064_11_BO78_OS.mat'

patch_mat_file = 'ST_80064_patch_C3.mat'

file_manifolds = 'manifolds_P/OS_BO78/'

JFField = AxisymmetricCylindricalGridField.from_matlab_file(f"{repository_path}{pert_mat_file}", with_perturbation=True)

section = CylindricalBfieldSection(JFField,phi0=1.9,R0=0.88, Z0=0)

top_o = FixedPoint(section)
top_o.find(1, [0.89,0.3], method='scipy.root')
top_o_coord = top_o.coords[0]
top_o.plot(ax=ax, marker='o', color="xkcd:crimson")


x_point1= FixedPoint(section)
x_point1.find(1, [0.7,-0.16], method='scipy.root')
x_point1_coord = x_point1.coords[0]
x_point1.plot(ax=ax, marker='x', color="xkcd:crimson")


######################################################################################### PERTURBED SCRIPT ######################################################################################



###################.  poincaré plot #####################

pplot = PoincarePlot.with_linspace(section, top_o_coord, x_point1_coord,40)
pplot.compute(400)
np.save(f"{repository_path}STD_80064_pp_hits.npy.npy", pplot._hits)


Hits=np.load(f"{repository_path}STD_80064_pp_hits.npy")
ax.scatter(Hits[:,:, 0], Hits[:,:, 1], color="xkcd:dark grey", s=1., linewidths=0)



###########################  tomographic reconstruction  #########################


tpc= tomo(f"{repository_path}{patch_mat_file}",ax,emi_vmin=0, emi_vmax=8e20)


###########################. computation #########################

#####top fp top mf#########

# manifold_1T = Manifold(section, x_point1, x_point1,-x_point1_coord+top_o_coord, -x_point1_coord+top_o_coord)
# manifold_1T.compute(eps_s=9e-6, eps_u=8e-6, nint_s=8, nint_u=12, neps_s=80, neps_u=240) #8 14 240
# manifold_1T.save(f"{repository_path}{file_manifolds}mf_1T.pkl")

# # #####top fp bottom mf#########

# manifold_1B = Manifold(section, x_point1, x_point1,x_point1_coord-top_o_coord, x_point1_coord-top_o_coord)
# manifold_1B.compute(eps_s=9e-6, eps_u=8e-6, nint_s=10, nint_u=12, neps_s=80, neps_u=80)
# manifold_1B.save(f"{repository_path}{file_manifolds}mf_1B.pkl")

####################### loading and plotting  #######################

manifold_1T  = Manifold.load(f"{repository_path}{file_manifolds}mf_1T.pkl")
manifold_1B  = Manifold.load(f"{repository_path}{file_manifolds}mf_1B.pkl")

manifold_1T.plot(stepsize_limit=0.3, ax=ax, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"])
manifold_1B.plot(stepsize_limit=0.3, ax=ax, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"])

top_o.plot(ax=ax, marker='o', color="xkcd:white")
x_point1.plot(ax=ax, marker='x', color="xkcd:white")

#####figure

# ax.set_xlim(0.62, 1.15)
# ax.set_ylim(-0.75, 0.75)
# ax.set_xlabel(r"$R[m]$")
# ax.set_ylabel(r"$Z[m]$")
# ax.set_aspect('equal') 
# ax.set_title('(80064 / 1.1s. / 1.9 rad)')
# plt.savefig(f"{repository_path}/figures/ST_Pert_80064_C3.png", bbox_inches='tight', dpi=720)
# plt.show()


#####zoomed



# ax.set_xlim(0.7, 1.0)
# ax.set_ylim(-0.7, -0.0)
# ax.set_xlabel(r"$R[m]$")
# ax.set_ylabel(r"$Z[m]$")
# ax.set_aspect('equal') 
# ax.set_title('(80064 / 1.1s. / 1.9 rad)')
# plt.savefig(f"{repository_path}/figures/ST_Pert_80064_C3_zoom.png", bbox_inches='tight', dpi=720)
# plt.show()





