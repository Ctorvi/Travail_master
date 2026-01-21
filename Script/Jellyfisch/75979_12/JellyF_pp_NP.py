import os, sys
proj_root = os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..','..'))  # remonte jusqu'à "Travail master"
if proj_root not in sys.path:
    sys.path.insert(0, proj_root)
from script.function.field_utils import plot_tomographic as tomo
from script.function.field_utils import plot_LIUQE as LIUQE


from pyoculus.fields import AxisymmetricCylindricalGridField
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
import numpy as np
from matplotlib import pyplot as plt
import logging
from scipy.io import loadmat
from matplotlib.lines import Line2D
from matplotlib.tri import Triangulation


fig, ax2 = plt.subplots(1, 1, figsize=(6, 6.6))


plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.size": 9,
    }
)


file_manifolds = 'manifolds_NP/'

logging.basicConfig(level=logging.DEBUG)

repository_path = './script/jellyfisch/75979_12/'

pert_mat_file = 'mat_files/JF_75979_120.mat'

patch_mat_file = './script/jellyfisch/75979_12/JF_75979_patch_120_He2.mat'


mf_list = ['mf_1T', 'mf_1B', 'mf_2T', 'mf_2B', 'mf_3L', 'mf_3R', 'mf_4T', 'mf_4B']

LIUQE_mat_file = 'mat_files/Liuqe_JF_75979_120.mat'

JFField = AxisymmetricCylindricalGridField.from_matlab_file(f"{repository_path}{pert_mat_file}", with_perturbation=False)

section = CylindricalBfieldSection(JFField,phi0=1.9,R0=0.88, Z0=0)

top_o = FixedPoint(section)
top_o.find(1, [0.9,0.18], method='scipy.root')
top_o_coord = top_o.coords[0]

x_point1= FixedPoint(section)
x_point1.find(1, [0.8,-0.27], method='scipy.root')
x_point1_coord = x_point1.coords[0]


x_point2= FixedPoint(section)
x_point2.find(1, [0.8,-0.57], method='scipy.root')
x_point2_coord = x_point2.coords[0]

x_point3= FixedPoint(section)
x_point3.find(1, [1.05,-0.58], method='scipy.root')
x_point3_coord = x_point3.coords[0]

x_point4= FixedPoint(section)
x_point4.find(1, [0.75,0.65], method='scipy.root')
x_point4_coord = x_point4.coords[0]


###################.  poincaré plot #####################

# pplot = PoincarePlot.with_linspace(section, top_o_coord, x_point1_coord,40)
# pplot.compute(400)
# np.save('./script/jellyfisch/75979_12/JF_75979_pp_hits_P.npy', pplot._hits)

Hits=np.load('./script/jellyfisch/75979_12/JF_75979_pp_hits_NP.npy')
ax2.scatter(Hits[:,:, 0], Hits[:,:, 1], color="xkcd:dark grey", s=1., linewidths=0)


###########################  tomographic reconstruction  #########################


#tpc= tomo(f"{patch_mat_file}",ax,emi_vmin=0, emi_vmax=2.25e19)


################# plot LIUQE ###################

data_LIUQE, lgd1,lgd2= LIUQE(f"{repository_path}{LIUQE_mat_file}",n_levels=15,ax=ax2,colors=["springgreen", "springgreen", "black"])

#data_LIUQE, lgd1b,lgd2b= LIUQE(f"{repository_path}{LIUQE_mat_file}",n_levels=15,ax=ax1, colors=["black", "black", "black"],linestyles=['--', '--', '-'], lw=[0.5,0.5,2])



###top fp top mf#########

# manifold_1T = Manifold(section, x_point1, x_point1,-x_point1_coord+top_o_coord, -x_point1_coord+top_o_coord)
# manifold_1T.compute(
#    eps_s=9e-6, eps_u=8e-6, nint_s=30, nint_u=21, neps_s=80, neps_u=240) #8 14 240

# manifold_1T.save(f"{repository_path}{file_manifolds}mf_1T_new.pkl")

# #####top fp bottom mf#########

# manifold_1B = Manifold(section, x_point1, x_point1)
# manifold_1B.compute(
#        eps_s=9e-6, eps_u=8e-6, nint_s=17, nint_u=18, neps_s=80, neps_u=80)
# manifold_1B.save(f"{repository_path}{file_manifolds}mf_1B.pkl")


# #### bottom fp top mf


# manifold_2T = Manifold(section, x_point2, x_point2,-x_point2_coord+top_o_coord, -x_point2_coord+top_o_coord)
# manifold_2T.compute(
#       eps_s=9e-6, eps_u=8e-6, nint_s=32, nint_u=30, neps_s=80, neps_u=80)
# manifold_2T.save(f"{repository_path}{file_manifolds}mf_2T.pkl")

# #### bottom fp bottom mf


# manifold_2B = Manifold(section, x_point2, x_point2,x_point2_coord-top_o_coord, x_point2_coord-top_o_coord)
# manifold_2B.compute(
#       eps_s=9e-6, eps_u=8e-6, nint_s=20, nint_u=20, neps_s=80, neps_u=80)
# manifold_2B.save(f"{repository_path}{file_manifolds}mf_2B.pkl")

# #### right fp left mf

# manifold_3L = Manifold(section, x_point3, x_point3,-x_point3_coord+top_o_coord, -x_point3_coord+top_o_coord)
# manifold_3L.compute(
#       eps_s=9e-6, eps_u=8e-6, nint_s=33, nint_u=30, neps_s=80, neps_u=80)
# manifold_3L.save(f"{repository_path}{file_manifolds}mf_3L.pkl")

# # ##### right fp right mf

# manifold_3R = Manifold(section, x_point3, x_point3,x_point3_coord-top_o_coord, x_point3_coord-top_o_coord)
# manifold_3R.compute(
#       eps_s=9e-6, eps_u=8e-6, nint_s=20, nint_u=30, neps_s=80, neps_u=80)
# manifold_3R.save(f"{repository_path}{file_manifolds}mf_3R.pkl")

# ## top fp top mf 

# manifold_4T = Manifold(section, x_point4, x_point4,x_point4_coord-top_o_coord, x_point4_coord-top_o_coord)
# manifold_4T.compute(
#      eps_s=9e-6, eps_u=8e-6, nint_s=8, nint_u=8, neps_s=80, neps_u=80)
# manifold_4T.save(f"{repository_path}{file_manifolds}mf_4T.pkl")

# # # # ##### over the top fp bottom mf


# manifold_4B = Manifold(section, x_point4, x_point4,-x_point4_coord+top_o_coord, -x_point4_coord+top_o_coord)
# manifold_4B.compute(
#      eps_s=9e-6, eps_u=8e-6, nint_s=8, nint_u=30, neps_s=80, neps_u=80)
# manifold_4B.save(f"{repository_path}{file_manifolds}mf_4B.pkl")

##### loading and plotting  ###

manifs = {name: {name} for name in mf_list}

for i in mf_list:
    manifs[i] = Manifold.load(f"{repository_path}{file_manifolds}{i}.pkl")
   # manifs[i].plot(ax=ax2, markersize=0, lw=0.7,labels=[None,None] if i!=0 else ['stable MF','unstable MF'])

    if i=='mf_1T':
     trajectory = manifs[i].stable
     Aire=manifs[i]._AdL_integral_points(trajectory)
    manifs[i].plot(stepsize_limit=0.05, ax=ax2, markersize=0, lw=0.5,colors=["xkcd:royal blue", "xkcd:red"],labels=['stable MF','unstable MF'] if i=='mf_1T' else [None,None])



top_o.plot(ax=ax2, marker='o',s=40, color="xkcd:dark blue",label=None, zorder=10)
x_point1.plot(ax=ax2, marker='x', s=50, color="xkcd:dark blue",label='X-points', zorder=10)
x_point2.plot(ax=ax2, marker='x', s=50, color="xkcd:dark blue",label=None, zorder=10)
x_point3.plot(ax=ax2, marker='x', s=50, color="xkcd:dark blue",label=None, zorder=10)
x_point4.plot(ax=ax2, marker='x', s=50, color="xkcd:dark blue",label=None, zorder=10)
ratio=2.8158

# top_o.plot(ax=ax2, marker='o', s=50, color="xkcd:green",label='O-point', zorder=10)
# x_point1.plot(ax=ax2, marker='x', s=70, color="xkcd:green",label='X-points', zorder=10)
# x_point2.plot(ax=ax2, marker='x', s=50, color="xkcd:green",label=None, zorder=10)
# x_point3.plot(ax=ax2, marker='x', s=50, color="xkcd:green",label=None, zorder=10)
# x_point4.plot(ax=ax2, marker='x', s=50, color="xkcd:green",label=None, zorder=10)
# ratio=2.8158



handles, labels = ax2.get_legend_handles_labels()
handles.append(lgd1)
handles.append(lgd2)
labels.append("LCFS")
labels.append("LIUQE")
ax2.legend(handles, labels,loc='lower right', bbox_to_anchor=(1.05, 0.0), ncol=1, fontsize=9)


ax2.set_xlim(0.62, 1.15)
ax2.set_ylim(-0.75, 0.75)
ax2.set_xlabel(r"$R[m]$")
ax2.set_ylabel(r"$Z[m]$")
ax2.set_aspect('equal') 
ax2.set_title('Poincaré and LIUQE')

#plt.savefig(f"{repository_path}figures/JF_75979_NP_pp_LIU.png", bbox_inches='tight', dpi=720)
plt.show()

print(Aire)
#####zoomed


# ax.set_xlim(0.7, 1.0)
# ax.set_ylim(-0.58, -0.08)
# ax.set_xlabel(r"$R[m]$")
# ax.set_ylabel(r"$Z[m]$")
# ax.set_aspect('equal') 
# ax.legend(loc='upper center', bbox_to_anchor=(0.5, 0.99), ncol=1, fontsize=9)
# ax.set_title('(75979 / 1.2s / 1.9 rad)')
# #plt.savefig(f"{repository_path}figures/Jellyfisch_75979_BO_Pert_mf_patch_zoom.png", dpi=150, bbox_inches=None, facecolor=fig.get_facecolor())
# plt.show()




