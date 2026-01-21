import os, sys
proj_root = os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..'))  # remonte jusqu'à "Travail master"
if proj_root not in sys.path:
    sys.path.insert(0, proj_root)
from script.function.field_utils import plot_tomographic as tomo
from script.function.field_utils import plot_vessel as vessel

from pyoculus.fields import AxisymmetricCylindricalGridField
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
import numpy as np
from matplotlib import pyplot as plt
import logging
from scipy.io import loadmat
from matplotlib.tri import Triangulation


def IndexClosestPoint(traj, point):
    """
    Return index of the point in traj closest to `point`.
    Robust to list/ndarray input, extra columns (e.g. phi), and empty traj.
    """
    import numpy as _np

    traj_arr = _np.asarray(traj)
    pt = _np.asarray(point).ravel()

    

    if traj_arr.size == 0:
        raise ValueError("traj is empty")

    # ensure traj is 2D: shape (N, d)
    if traj_arr.ndim == 1:
        # if traj is a flat array of same length as point -> single point
        if traj_arr.size == pt.size:
            traj_arr = traj_arr.reshape(1, -1)
        else:
            raise ValueError(f"traj is 1D of length {traj_arr.size}, point length {pt.size}")

    # if traj has extra columns (e.g. R, phi, Z), try to use first len(pt) cols
    if traj_arr.shape[1] != pt.size:
        if traj_arr.shape[1] > pt.size:
            traj_proc = traj_arr[:, :pt.size]
        else:
            raise ValueError(f"traj points have dimension {traj_arr.shape[1]} but point has length {pt.size}")
    else:
        traj_proc = traj_arr

    # compute squared distances (faster and avoid needless sqrt)
    d2 = _np.sum((traj_proc - pt) ** 2, axis=1)
    idx = int(_np.argmin(d2))
    return idx

#fig, ax = plt.subplots(1, 1, figsize=(6, 6.6))
fig,ax = plt.subplots(1,1, figsize=(4,6))
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

#JFField = AxisymmetricCylindricalGridField.from_matlab_file(f"{repository_path}{pert_mat_file}", with_perturbation=True)

JFField = AxisymmetricCylindricalGridField.from_matlab_file(f"{repository_path}{pert_mat_file}", with_perturbation=False)

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

# pplot = PoincarePlot.with_linspace(section, top_o_coord, x_point1_coord,40)
# pplot.compute(400)
# np.save(f"{repository_path}STD_80064_pp_NP_hits.npy", pplot._hits)



###########################  tomographic reconstruction  #########################


# tpc= tomo(f"{repository_path}{patch_mat_file}",ax,emi_vmin=0, emi_vmax=8e20)


###########################. computation #########################

#####top fp top mf#########

# manifold_1T = Manifold(section, x_point1, x_point1,-x_point1_coord+top_o_coord, -x_point1_coord+top_o_coord)
# manifold_1T.compute(eps_s=9e-6, eps_u=8e-6, nint_s=10, nint_u=12, neps_s=80, neps_u=240) #8 14 240
# manifold_1T.save(f"{repository_path}{file_manifolds}mf_1T_ns10.pkl")

# # #####top fp bottom mf#########

# manifold_1B = Manifold(section, x_point1, x_point1,x_point1_coord-top_o_coord, x_point1_coord-top_o_coord)
# manifold_1B.compute(eps_s=9e-6, eps_u=8e-6, nint_s=10, nint_u=12, neps_s=80, neps_u=80)
# manifold_1B.save(f"{repository_path}{file_manifolds}mf_1B_NP.pkl")

####################### loading and plotting  #######################



###### perturbed ######
# manifold_1T  = Manifold.load(f"{repository_path}{file_manifolds}mf_1T_ns10.pkl")
# manifold_1B  = Manifold.load(f"{repository_path}{file_manifolds}mf_1B.pkl")
# # manifold_1T.find_clinics(first_guess_eps_s=9e-6, first_guess_eps_u=8e-6)
# A=manifold_1T.compute_turnstile_areas()

# # manifold_1T.save(f"{repository_path}{file_manifolds}mf_1T_ns10.pkl")

# # manifold_1T.plot_clinics(ax=ax)

# Hits=np.load(f"{repository_path}ST_80064_pp_hits.npy")
# ax.scatter(Hits[:,:, 0], Hits[:,:, 1], color="xkcd:dark grey", s=1., linewidths=0,zorder=10)


###### non perturbed ######

manifold_1T  = Manifold.load(f"{repository_path}{file_manifolds}mf_1T_NP.pkl")
manifold_1B  = Manifold.load(f"{repository_path}{file_manifolds}mf_1B_NP.pkl")

Hits=np.load(f"{repository_path}ST_80064_pp_NP_hits.npy")
ax.scatter(Hits[:,:, 0], Hits[:,:, 1], color="xkcd:dark grey", s=1., linewidths=0,zorder=10)

#################


manifold_1T.plot(stepsize_limit=0.3, ax=ax, markersize=0, lw=3)
manifold_1B.plot(stepsize_limit=0.3, ax=ax, markersize=0, lw=2,labels=[None, None])


top_o.plot(ax=ax, marker='o',s=70, color="xkcd:dark blue",label='O-points', zorder=10)
x_point1.plot(ax=ax, marker='x', s=70, color="xkcd:dark blue",label='X-points', zorder=10)

#####figure
vessel(ax)

ax.legend(loc='center right', bbox_to_anchor=(1.02, 0.33), ncol=1, fontsize=7)

ax.set_xlim(0.62, 1.15)
ax.set_ylim(-0.75, 0.75)

ax.set_xlabel(r"$R[m]$")
ax.set_ylabel(r"$Z[m]$")
ax.set_aspect('equal') 
ax.set_title('(80064 / 1.1s.)')
plt.savefig(f"{repository_path}/figures/ST_NoPert_80064.png",  dpi=720,bbox_inches="tight", pad_inches=0, facecolor=fig.get_facecolor())
plt.show()


#####zoomed


# ax.set_xlim(0.7, 1.0)
# ax.set_ylim(-0.7, -0.0)
# ax.set_xlabel(r"$R[m]$")
# ax.set_ylabel(r"$Z[m]$")
# ax.set_aspect('equal') 
# ax.set_title('(80064 / 1.1s. / 1.9 rad)')
# plt.savefig(f"{repository_path}/figures/ST_Pert_80064_C3_zoom.png", bbox_inches='tight', dpi=720)
# plt.show()


######### tangle ########

# ax.set_xlim(0.61, 0.9)
# ax.set_ylim(-0.4, 0.1)

# ax.set_xlabel(r"$R[m]$")
# ax.set_ylabel(r"$Z[m]$")
# ax.set_aspect('equal') 
# ax.set_title('(80064 / 1.1s.)')
# #plt.savefig(f"{repository_path}figures/ST_tangle.png", dpi=150,bbox_inches="tight", pad_inches=0, facecolor=fig.get_facecolor())
# plt.show()




######## calculation #######

# trajectory_s = manifold_1T._stable_trajectory
# trajectory_u = manifold_1T._unstable_trajectory

# mask_s = np.isfinite(trajectory_s[:,0]) | np.isfinite(trajectory_s[:,1])
# mask_u = np.isfinite(trajectory_u[:,0]) | np.isfinite(trajectory_u[:,1])

# M1_0_idx=len(manifold_1T.clinics[0].trajectory)//2

# M1_0 = manifold_1T.clinics[0].trajectory[M1_0_idx]

# trajectory_s=trajectory_s[mask_s]
# trajectory_u=trajectory_u[mask_u]
       
# idx_1s = IndexClosestPoint(trajectory_s, M1_0)
# idx_1u = IndexClosestPoint(trajectory_u , M1_0)

# trajectory_u= trajectory_u[:idx_1u]  
# trajectory_s= trajectory_s[:idx_1s]  

# idx_1s = IndexClosestPoint(trajectory_s, M1_0)
# idx_1u = IndexClosestPoint(trajectory_u, M1_0)

# trajectory_u= trajectory_u[:idx_1u]  
# trajectory_s= trajectory_s[:idx_1s] 

# trajectory_u= trajectory_u[::-1]      
# contour_points = np.vstack([
#          manifold_1T.rfp_s,      
#          trajectory_s,              
#          trajectory_u,      
#  ])       

# Aire=manifold_1T._AdL_integral_points(contour_points)
# print(Aire)



