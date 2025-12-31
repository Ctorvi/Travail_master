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

fig, ax2 = plt.subplots(1, 1, figsize=(5, 8))


plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.size": 9,
    }
)


file_manifolds = 'manifolds_P/OS_BO78/'


logging.basicConfig(level=logging.DEBUG)

repository_path = './script/jellyfisch/75979_12/'

pert_mat_file = 'mat_files/JF_75979_120_BO78_OS.mat'

mf_list = ['mf_1T_nu22', 'mf_1B', 'mf_2T', 'mf_2B', 'mf_3L', 'mf_3R', 'mf_4T', 'mf_4B']

patch_mat_file = './script/jellyfisch/75979_12/mat_files/JF_75979_patch_120_He2.mat'

JFField = AxisymmetricCylindricalGridField.from_matlab_file(f"{repository_path}{pert_mat_file}", with_perturbation=True)

section = CylindricalBfieldSection(JFField,phi0=1.9,R0=0.88, Z0=0)



######## turnstile plot and turnstile flux calculation #########


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


##########################


Hits=np.load('./script/jellyfisch/75979_12/JF_75979_pp_hits_P.npy')
ax2.scatter(Hits[:,:, 0], Hits[:,:, 1], color="xkcd:dark grey", s=1., linewidths=0)



manifs = {name: {name} for name in mf_list}

for i in mf_list:
     manifs[i] = Manifold.load(f"{repository_path}{file_manifolds}{i}.pkl")
    # manifs[i].plot(ax=ax2, markersize=0, lw=0.7,labels=[None,None] if i!=0 else ['stable MF','unstable MF'])
     if i=='mf_1T_nu22':   
        manifs[i].find_clinics(first_guess_eps_s=9e-6, first_guess_eps_u=8e-6)
        A=manifs[i].compute_turnstile_areas()
        print(f"Turnstile area for {i} : {A} m^2")
        manifs[i].save(f"{repository_path}{file_manifolds}{i}.pkl")
        manifs[i].plot_clinics(ax=ax2)
     manifs[i].plot(stepsize_limit=0.05, ax=ax2, markersize=0, lw=0.5,colors=["xkcd:royal blue", "xkcd:red"],labels=['stable MF','unstable MF'] if i=='mf_1T' else [None,None])



top_o.plot(ax=ax2, marker='o', s=50, color="xkcd:green",label='O-point', zorder=10)
x_point1.plot(ax=ax2, marker='x', s=70, color="xkcd:green",label='X-points', zorder=10)
x_point2.plot(ax=ax2, marker='x', s=50, color="xkcd:green",label=None, zorder=10)
x_point3.plot(ax=ax2, marker='x', s=50, color="xkcd:green",label=None, zorder=10)
x_point4.plot(ax=ax2, marker='x', s=50, color="xkcd:green",label=None, zorder=10)
ratio=2.8158


ax2.set_xlim(0.62, 1.15)
ax2.set_ylim(-0.75, 0.75)
ax2.set_xlabel(r"$R[m]$")
ax2.set_ylabel(r"$Z[m]$")
ax2.set_aspect('equal') 
ax2.set_title('Poincaré perturbed')

ax2.legend(loc='lower right', bbox_to_anchor=(1., 0.86), ncol=1, fontsize=9)

#plt.savefig(f"{repository_path}figures/JF_75979_P_pp.png", bbox_inches='tight', dpi=720)
plt.show()




########### resonnance zone #############



# top_o = FixedPoint(section)
# top_o.find(1, [0.9,0.18], method='scipy.root')
# top_o_coord = top_o.coords[0]


# o_point_island= FixedPoint(section)
# o_point_island.find(1, [0.99,0.199], method='scipy.root')
# o_point_island_coord = o_point_island.coords[0]

# # x_point_island_top= FixedPoint(section)
# # x_point_island_top.find(1, [0.85,0.3], method='scipy.root')
# # x_point_island_top_coord = x_point_island_top.coords[0]

# # x_point_island_bottom= FixedPoint(section)
# # x_point_island_bottom.find(1, [0.85,0.06], method='scipy.root')
# # x_point_island_bottom_coord = x_point_island_bottom.coords[0]



# x_point_island_n= FixedPoint(section)
# x_point_island_n.find(1, [0.79,0.18], method='scipy.root')
# x_point_island_n_coord = x_point_island_n.coords[0]


# o_point_island.plot(ax=ax2, marker='o',s=40, color="xkcd:dark blue",label='O-points', zorder=10)
# # x_point_island_top.plot(ax=ax2, marker='x', s=50, color="xkcd:dark blue",label='X-points', zorder=10)
# x_point_island_n.plot(ax=ax2, marker='x', s=50, color="xkcd:dark blue",label=None, zorder=10)
# #x_point_island_bottom.plot(ax=ax2, marker='x', s=50, color="xkcd:dark blue",label=None, zorder=10)


# #####top fp bottom mf#########

# manifold_1T = Manifold(section, x_point_island_n, x_point_island_n, [0,1],[0,1] )
# manifold_1T.compute(
#        eps_s=9e-6, eps_u=8e-6, nint_s=80, nint_u=80, neps_s=40, neps_u=40)



# manifold_1T = Manifold.load(f"{repository_path}{file_manifolds}mf_1_island_n.pkl")


# manifold_1T.find_clinics(first_guess_eps_s=8e-6, first_guess_eps_u=9e-8)


# # manifold_1T.save(f"{repository_path}{file_manifolds}mf_1_island_n.pkl")



# manifold_1T.plot(ax=ax2, stepsize_limit=0.1, linewidth=0.7, markersize=0)

# # manifold_1B = Manifold(section, x_point_island_n, x_point_island_n, [0,-1],[0,-1] )
# # manifold_1B.compute(
# #        eps_s=9e-6, eps_u=8e-6, nint_s=80, nint_u=80, neps_s=40, neps_u=40)

# manifold_1B = Manifold.load(f"{repository_path}{file_manifolds}mf_1_island_n2.pkl")

# manifold_1B.find_clinics(first_guess_eps_s=9e-5, first_guess_eps_u=8e-5)




# # manifold_1B.save(f"{repository_path}{file_manifolds}mf_1_island_n2.pkl")

# manifold_1B.plot(ax=ax2, stepsize_limit=0.1, linewidth=0.7, markersize=0)

# manifold_1B = Manifold(section, x_point_island_bottom, x_point_island_bottom,-x_point_island_bottom.coords[0]+o_point_island.coords[0], -x_point_island_bottom.coords[0]+o_point_island.coords[0])
# manifold_1B.compute(
#        eps_s=9e-6, eps_u=8e-6, nint_s=17, nint_u=18, neps_s=80, neps_u=80)
# manifold_1B.save(f"{repository_path}{file_manifolds}mf_2_island.pkl")
# manifold_1B.plot(ax=ax2, stepsize_limit=0.1, linewidth=0.7, markersize=0, color="xkcd:red", label='manifold bottom FP')



# ax2.set_xlim(0.75, 1.0)
# ax2.set_ylim(0.0, 0.35)
# ax2.set_xlabel(r"$R[m]$")
# ax2.set_ylabel(r"$Z[m]$")
# ax2.set_aspect('equal') 
# ax2.set_title('Poincaré perturbed')

# ax2.legend(loc='lower right', bbox_to_anchor=(1., 0.86), ncol=1, fontsize=9)

# #plt.savefig(f"{repository_path}figures/JF_75979_P_pp.png", bbox_inches='tight', dpi=720)
# plt.show()


