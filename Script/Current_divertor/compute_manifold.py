from pyoculus.fields import AxisymmetricCylindricalGridField
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
import numpy as np
from matplotlib import pyplot as plt
import logging
from scipy.io import loadmat
from matplotlib.tri import Triangulation
from scipy.special import ellipk, ellipe


DATA=np.load('./Script/Current_divertor/Current_div_B.npz')

B_R_tot=DATA['B_R_tot']
B_Z_tot=DATA['B_Z_tot']
psi=DATA['psi']



fig, ax = plt.subplots(1, 1, figsize=(6, 6.6))


plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.size": 9,
    }
)


logging.basicConfig(level=logging.DEBUG)

repository_path='./Script/Current_divertor/'


JFField = AxisymmetricCylindricalGridField.from_matlab_file_with_separatrix_current('./Script/Jellyfisch/77062_12/normal/JF_77062_12.mat', B_R_tot, B_Z_tot, psi, with_perturbation=True)
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


# pplot = PoincarePlot.with_linspace(section, top_o_coord, x_point1_coord,40)
# pplot.compute(400)
# np.save('./Script/Jellyfisch/77021_120/poincare_hits/jellyfisch_Pert_120.npy', pplot._hits)




###########################  tomographic reconstruction  #########################





# data = loadmat(f"{repository_path}{patch_mat_file}", squeeze_me=True, struct_as_record=False)

# inv_grid=data['inv_grid']

# tri_x = inv_grid.tri_x
# tri_y = inv_grid.tri_y
# tri_nodes = inv_grid.tri_nodes
# R_values = inv_grid.R_values
# Z_values = inv_grid.Z_values
# R_mean = inv_grid.R_mean
# Z_mean = inv_grid.Z_mean
# Area= inv_grid.Area

# emi=data['emi']
# t=data['t']
# time=data['time']
# time_vec = np.asarray(time)

# # Find index of closest time
# t_idx = np.argmin(np.abs(time_vec - t))
# emi = emi[:, t_idx]
# triangles = tri_nodes.astype(int)  # tri_nodes doit être zéro-indexé
# tri = Triangulation(tri_x, tri_y, triangles)
# tpc=ax.tripcolor(tri, emi, shading='flat', edgecolors='none', vmin=0, vmax=8e20)
# #fig.colorbar(tpc, ax=ax, label='Émission')
   

###########################. fixed point  #########################


#####top fp top mf#########

manifold_1T = Manifold(section, x_point1, x_point1,-x_point1_coord+top_o_coord, -x_point1_coord+top_o_coord)
manifold_1T.compute(
   eps_s=9e-6, eps_u=8e-6, nint_s=10, nint_u=14, neps_s=80, neps_u=240) #8 14 240

manifold_1T.save(f"{repository_path}manifolds_P/mf_1T.pkl")

# #####top fp bottom mf#########

# manifold_1B = Manifold(section, x_point1, x_point1)
# manifold_1B.compute(
#        eps_s=9e-6, eps_u=8e-6, nint_s=10, nint_u=12, neps_s=80, neps_u=80)
# manifold_1B.save(f"{repository_path}manifolds_P/mf_1B.pkl")


##### bottom fp top mf


# manifold_2T = Manifold(section, x_point2, x_point2,-x_point2_coord+top_o_coord, -x_point2_coord+top_o_coord)
# manifold_2T.compute(
#       eps_s=9e-6, eps_u=8e-6, nint_s=32, nint_u=30, neps_s=80, neps_u=80)
# manifold_2T.save(f"{repository_path}manifolds_P/mf_2T.pkl")

##### bottom fp bottom mf


# manifold_2B = Manifold(section, x_point2, x_point2,x_point2_coord-top_o_coord, x_point2_coord-top_o_coord)
# manifold_2B.compute(
#       eps_s=9e-6, eps_u=8e-6, nint_s=8, nint_u=8, neps_s=80, neps_u=80)
# manifold_2B.save(f"{repository_path}manifolds_P/mf_2B.pkl")

##### right fp left mf

# manifold_3L = Manifold(section, x_point3, x_point3,-x_point3_coord+top_o_coord, -x_point3_coord+top_o_coord)
# manifold_3L.compute(
#       eps_s=9e-6, eps_u=8e-6, nint_s=33, nint_u=30, neps_s=80, neps_u=80)
# manifold_3L.save(f"{repository_path}manifolds_P/mf_3L.pkl")

# # ##### right fp right mf

# manifold_3R = Manifold(section, x_point3, x_point3,x_point3_coord-top_o_coord, x_point3_coord-top_o_coord)
# manifold_3R.compute(
#       eps_s=9e-6, eps_u=8e-6, nint_s=20, nint_u=20, neps_s=80, neps_u=80)
# manifold_3R.save(f"{repository_path}manifolds_P/mf_3R.pkl")

### top fp top mf 

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



manifold_init = Manifold.load('./Script/Jellyfisch/77062_12/normal/manifolds_P/mf_1T.pkl')

#
#manifold_1T  = Manifold.load(f"{repository_path}manifolds_P/mf_1T.pkl")
# manifold_1B  = Manifold.load(f"{repository_path}manifolds_P/mf_1B.pkl")
# manifold_2T  = Manifold.load(f"{repository_path}manifolds_P/mf_2T.pkl")
# manifold_2B  = Manifold.load(f"{repository_path}manifolds_P/mf_2B.pkl")
# manifold_3L = Manifold.load(f"{repository_path}manifolds_P/mf_3L.pkl")
# manifold_3R  = Manifold.load(f"{repository_path}manifolds_P/mf_3R.pkl")
# manifold_4T = Manifold.load(f"{repository_path}manifolds_P/mf_4T.pkl")
# manifold_4B  = Manifold.load(f"{repository_path}manifolds_P/mf_4B.pkl")

manifold_init.plot(stepsize_limit=0.1, ax=ax, markersize=0, lw=0.7,colors=["lightgrey", "darkgrey"])

manifold_1T.plot(stepsize_limit=0.1, ax=ax, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"])
# manifold_1B.plot(stepsize_limit=0.3,ax=ax, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"])
# manifold_2T.plot(stepsize_limit=0.1, ax=ax, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"])
# manifold_2B.plot(stepsize_limit=0.3, ax=ax, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"])
# manifold_3L.plot(stepsize_limit=0.1,ax=ax, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"])
# manifold_3R.plot(ax=ax, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"])
# manifold_4T.plot(ax=ax, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"])
# manifold_4B.plot(ax=ax, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"])




top_o.plot(ax=ax, marker='o', color="xkcd:black")
x_point1.plot(ax=ax, marker='x', color="xkcd:black")
x_point2.plot(ax=ax, marker='x', color="xkcd:black")
x_point3.plot(ax=ax, marker='x', color="xkcd:black")
x_point4.plot(ax=ax, marker='x', color="xkcd:black")
ratio=2.8158

#Hits=np.load('./Script/Jellyfisch/77021_120/poincare_hits/jellyfisch_Pert_test.npy')
# ax.scatter(Hits[:,:, 0], Hits[:,:, 1], color="xkcd:dark grey", s=1.4, linewidths=0)

#####figure



ax.set_xlim(0.7, 1.0)
ax.set_ylim(-0.6, -0.1)
ax.set_xlabel(r"$R[m]$")
ax.set_ylabel(r"$Z[m]$")
ax.set_aspect('equal') 
ax.set_title('77062 at 1.2s, with current in separatrix-strike-point (red)')
#plt.savefig(f"{repository_path}figures/Jellyfisch_77062_BO_Pert_mf_patch.png", bbox_inches='tight', dpi=720)
plt.show()


#####zoomed





