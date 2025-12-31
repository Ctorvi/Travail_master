from pyoculus.fields import AxisymmetricCylindricalGridField
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
import numpy as np
from matplotlib import pyplot as plt
import logging
from scipy.io import loadmat
from matplotlib.tri import Triangulation


fig, ax = plt.subplots(1, 1, figsize=(5, 8))
#fig, ax = plt.subplots(1, 1, figsize=(6, 6.6))


plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.size": 9,
    }
)


logging.basicConfig(level=logging.INFO)

mf_list = ['mf_1T', 'mf_1B', 'mf_2T', 'mf_2B', 'mf_3L', 'mf_3R', 'mf_4T', 'mf_4B']


repository_path = './Script/Jellyfisch/77062_12/KINETIC/'

pert_mat_file = '77062_12_KINETIC.mat'

patch_mat_file = './Script/Jellyfisch/77062_12/JF_77062_patch_12_C3.mat'

JFField = AxisymmetricCylindricalGridField.from_matlab_file(f"{repository_path}{pert_mat_file}", with_perturbation=True)


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


###########################  tomographic reconstruction  #########################

# tpc= tomo(f"{patch_mat_file}",ax,emi_vmin=0, emi_vmax=9e20)


###########################. fixed point  #########################


# ####top fp top mf#########

# manifold_1T = Manifold(section, x_point1, x_point1,-x_point1_coord+top_o_coord, -x_point1_coord+top_o_coord)
# manifold_1T.compute(
#    eps_s=9e-6, eps_u=8e-6, nint_s=10, nint_u=14, neps_s=80, neps_u=240) #8 14 240

# manifold_1T.save(f"{repository_path}manifolds_P/mf_1T_add.pkl")

# #####top fp bottom mf#########

# manifold_1B = Manifold(section, x_point1, x_point1)
# manifold_1B.compute(
#        eps_s=9e-6, eps_u=8e-6, nint_s=10, nint_u=12, neps_s=80, neps_u=80)
# manifold_1B.save(f"{repository_path}manifolds_P/mf_1B.pkl")


# #### bottom fp top mf


# manifold_2T = Manifold(section, x_point2, x_point2,-x_point2_coord+top_o_coord, -x_point2_coord+top_o_coord)
# manifold_2T.compute(
#       eps_s=9e-6, eps_u=8e-6, nint_s=32, nint_u=30, neps_s=80, neps_u=80)
# manifold_2T.save(f"{repository_path}manifolds_P/mf_2T.pkl")

# #### bottom fp bottom mf


# manifold_2B = Manifold(section, x_point2, x_point2,x_point2_coord-top_o_coord, x_point2_coord-top_o_coord)
# manifold_2B.compute(
#       eps_s=9e-6, eps_u=8e-6, nint_s=32, nint_u=32, neps_s=80, neps_u=80)
# manifold_2B.save(f"{repository_path}manifolds_P/mf_2B.pkl")

# #### right fp left mf

# manifold_3L = Manifold(section, x_point3, x_point3,-x_point3_coord+top_o_coord, -x_point3_coord+top_o_coord)
# manifold_3L.compute(
#       eps_s=9e-6, eps_u=8e-6, nint_s=33, nint_u=30, neps_s=80, neps_u=80)
# manifold_3L.save(f"{repository_path}manifolds_P/mf_3L.pkl")

# # ##### right fp right mf

# manifold_3R = Manifold(section, x_point3, x_point3,x_point3_coord-top_o_coord, x_point3_coord-top_o_coord)
# manifold_3R.compute(
#       eps_s=9e-6, eps_u=8e-6, nint_s=20, nint_u=20, neps_s=80, neps_u=80)
# manifold_3R.save(f"{repository_path}manifolds_P/mf_3R.pkl")

# ## top fp top mf 

# manifold_4T = Manifold(section, x_point4, x_point4,x_point4_coord-top_o_coord, x_point4_coord-top_o_coord)
# manifold_4T.compute(
#      eps_s=9e-6, eps_u=8e-6, nint_s=8, nint_u=8, neps_s=80, neps_u=80)
# manifold_4T.save(f"{repository_path}manifolds_P/mf_4T.pkl")

# # # # ##### over the top fp bottom mf


# manifold_4B = Manifold(section, x_point4, x_point4,-x_point4_coord+top_o_coord, -x_point4_coord+top_o_coord)
# manifold_4B.compute(
#      eps_s=9e-6, eps_u=8e-6, nint_s=8, nint_u=30, neps_s=80, neps_u=80)
# manifold_4B.save(f"{repository_path}manifolds_P/mf_4B.pkl")





##### loading and plotting  ###



#manifold_NK= Manifold.load('./Script/Jellyfisch/77062_12/manifolds_P/mf_1T.pkl')
manifs = {name: {name} for name in mf_list}

for i in mf_list:
     manifs[i] = Manifold.load(f"{repository_path}manifolds_P/{i}.pkl")
     manifs[i].plot(stepsize_limit=0.05, ax=ax, markersize=0, lw=0.6,colors=["rosybrown", "xkcd:red"],labels=['stable MF','unstable MF'] if i=='mf_1T' else [None,None])


#manifold_NK.plot(stepsize_limit=0.1, ax=ax, markersize=0, lw=0.7,colors=["grey", "xkcd:blue"])





top_o.plot(ax=ax, marker='o', color="xkcd:white",label=None)
x_point1.plot(ax=ax, marker='x', color="xkcd:white",label='X-points')
x_point2.plot(ax=ax, marker='x', color="xkcd:white",label=None)
x_point3.plot(ax=ax, marker='x', color="xkcd:white",label=None)
x_point4.plot(ax=ax, marker='x', color="xkcd:white",label=None)
ratio=2.8158

#Hits=np.load('./Script/Jellyfisch/77021_120/poincare_hits/jellyfisch_Pert_test.npy')
# ax.scatter(Hits[:,:, 0], Hits[:,:, 1], color="xkcd:dark grey", s=1.4, linewidths=0)

#####figure


# ax.set_xlim(0.62, 1.15)
# ax.set_ylim(-0.75, 0.75)
# ax.set_xlabel(r"$R[m]$")
# ax.set_ylabel(r"$Z[m]$")
# ax.set_aspect('equal') 
# ax.set_title('KINETIC 77062 at 1.2s, at 1.9 rad')
# plt.savefig(f"{repository_path}figures/Jellyfisch_77062_KINETIC_BO_Pert_mf_patch.png", bbox_inches='tight', dpi=720)
# plt.show()


#####zoomed



ax.set_xlim(0.7, 1.0)
ax.set_ylim(-0.58, -0.08)
ax.set_xlabel(r"$R[m]$")
ax.set_ylabel(r"$Z[m]$")
ax.set_aspect('equal') 
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 0.99), ncol=1, fontsize=9)
ax.set_title('KINETIC (77062 / 1.2s / 1.9 rad)')
#plt.savefig(f"{repository_path}figures/Jellyfisch_77062_KINETIC_BO_Pert_mf_patch_zoom.png", dpi=150, bbox_inches=None, facecolor=fig.get_facecolor())
plt.show()





