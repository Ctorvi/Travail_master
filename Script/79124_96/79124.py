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

repository_path = './Script/79124_96/'
pert_mat_file = '79124_96.mat'
patch_mat_file = '79124_patch_96_C3.mat'

JFField = AxisymmetricCylindricalGridField.from_matlab_file(f"{repository_path}{pert_mat_file}", with_perturbation=True)

section = CylindricalBfieldSection(JFField,phi0=1.9,R0=0.88, Z0=0)

top_o = FixedPoint(section)
top_o.find(1, [0.9,0.], method='scipy.root')
top_o_coord = top_o.coords[0]



x_point1= FixedPoint(section)
x_point1.find(1, [0.8,-0.4], method='scipy.root')
x_point1_coord = x_point1.coords[0]



x_point2= FixedPoint(section)
x_point2.find(1, [0.9,-0.65], method='scipy.root')
x_point2_coord = x_point2.coords[0]



# pplot = PoincarePlot.with_linspace(section, top_o_coord, x_point1_coord,40)
# pplot.compute(400)
# np.save('./Script/Jellyfisch/77021_120/poincare_hits/jellyfisch_Pert_120.npy', pplot._hits)




###########################  tomographic reconstruction  #########################





data = loadmat(f"{repository_path}{patch_mat_file}", squeeze_me=True, struct_as_record=False)

inv_grid=data['inv_grid']

tri_x = inv_grid.tri_x
tri_y = inv_grid.tri_y
tri_nodes = inv_grid.tri_nodes
R_values = inv_grid.R_values
Z_values = inv_grid.Z_values
R_mean = inv_grid.R_mean
Z_mean = inv_grid.Z_mean
Area= inv_grid.Area

emi=data['emi']
t=data['t']
time=data['time']
time_vec = np.asarray(time)

# Find index of closest time
t_idx = np.argmin(np.abs(time_vec - t))
emi = emi[:, t_idx]
triangles = tri_nodes.astype(int)  # tri_nodes doit être zéro-indexé
tri = Triangulation(tri_x, tri_y, triangles)
tpc=ax.tripcolor(tri, emi, shading='flat', edgecolors='none', vmin=0, vmax=1.5e21)
#fig.colorbar(tpc, ax=ax, label='Émission')
   

###########################. fixed point  #########################


####top fp top mf#########

# manifold_1T = Manifold(section, x_point1, x_point1,-x_point1_coord+top_o_coord, -x_point1_coord+top_o_coord)
# manifold_1T.compute(
#    eps_s=9e-6, eps_u=8e-6, nint_s=10, nint_u=14, neps_s=80, neps_u=240) #8 14 240

# manifold_1T.save(f"{repository_path}manifolds_P/mf_1T.pkl")

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
#       eps_s=9e-6, eps_u=8e-6, nint_s=8, nint_u=8, neps_s=80, neps_u=80)
# manifold_2B.save(f"{repository_path}manifolds_P/mf_2B.pkl")










##### loading and plotting  ###


manifold_1T  = Manifold.load(f"{repository_path}manifolds_P/mf_1T.pkl")
manifold_1B  = Manifold.load(f"{repository_path}manifolds_P/mf_1B.pkl")
manifold_2T  = Manifold.load(f"{repository_path}manifolds_P/mf_2T.pkl")
manifold_2B  = Manifold.load(f"{repository_path}manifolds_P/mf_2B.pkl")




manifold_1T.plot(stepsize_limit=0.1, ax=ax, markersize=0, lw=0.7,colors=["magenta", "xkcd:red"])
manifold_1B.plot(stepsize_limit=0.3,ax=ax, markersize=0, lw=0.7,colors=["magenta", "xkcd:red"])
manifold_2T.plot(stepsize_limit=0.1, ax=ax, markersize=0, lw=0.7,colors=["magenta", "xkcd:red"])
manifold_2B.plot(stepsize_limit=0.3, ax=ax, markersize=0, lw=0.7,colors=["magenta", "xkcd:red"])





top_o.plot(ax=ax, marker='o', color="xkcd:white")
x_point1.plot(ax=ax, marker='x', color="xkcd:white")
x_point2.plot(ax=ax, marker='x', color="xkcd:white")


ratio=2.8158

#Hits=np.load('./Script/Jellyfisch/77021_120/poincare_hits/jellyfisch_Pert_test.npy')
# ax.scatter(Hits[:,:, 0], Hits[:,:, 1], color="xkcd:dark grey", s=1.4, linewidths=0)

#####figure


# ax.set_xlim(0.62, 1.15)
# ax.set_ylim(-0.75, 0.75)
# ax.set_xlabel(r"$R[m]$")
# ax.set_ylabel(r"$Z[m]$")
# ax.set_aspect('equal') 
# ax.set_title('79124 at 0.96s, at 1.9 rad')
# plt.savefig(f"{repository_path}figures/79124_BO_Pert_mf_patch.png", bbox_inches='tight', dpi=720)
# plt.show()


#####zoomed



ax.set_xlim(0.7, 1.15)
ax.set_ylim(-0.6, -0.1)
ax.set_xlabel(r"$R[m]$")
ax.set_ylabel(r"$Z[m]$")
ax.set_aspect('equal') 
ax.set_title('79124 at 0.96s, at 1.9 rad')
plt.savefig(f"{repository_path}figures/79124_BO_Pert_mf_patch_zoom.png", bbox_inches='tight', dpi=720)
plt.show()





