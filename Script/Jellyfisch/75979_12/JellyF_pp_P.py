import os, sys
proj_root = os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..','..'))  # remonte jusqu'à "Travail master"
if proj_root not in sys.path:
    sys.path.insert(0, proj_root)
from script.function.field_utils import plot_tomographic as tomo
from script.function.field_utils import plot_LIUQE as LIUQE

from matplotlib.ticker import MaxNLocator
from pyoculus.fields import AxisymmetricCylindricalGridField
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
import numpy as np
from matplotlib import pyplot as plt
import logging
from scipy.io import loadmat
from matplotlib.lines import Line2D
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



fig, ax2 = plt.subplots(1, 1, figsize=(4, 2.6))
#fig, ax2 = plt.subplots(1, 1, figsize=(5, 5))

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

mf_list = ['mf_1T', 'mf_1B', 'mf_2T', 'mf_2B', 'mf_3L', 'mf_3R', 'mf_4T', 'mf_4B']

patch_mat_file = './script/jellyfisch/75979_12/mat_files/JF_75979_patch_120_He2.mat'

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

o_point_island= FixedPoint(section)
o_point_island.find(1, [0.99,0.199], method='scipy.root')
o_point_island_coord = o_point_island.coords[0]


x_point_island_n= FixedPoint(section)
x_point_island_n.find(1, [0.79,0.18], method='scipy.root')
x_point_island_n_coord = x_point_island_n.coords[0]


###top fp top mf#########

# manifold_1T = Manifold(section, x_point1, x_point1,-x_point1_coord+top_o_coord, -x_point1_coord+top_o_coord)
# manifold_1T.compute(
#    eps_s=9e-6, eps_u=8e-6, nint_s=30, nint_u=22, neps_s=80, neps_u=240) #8 14 240

# manifold_1T.save(f"{repository_path}{file_manifolds}mf_1T_nu22.pkl")


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

###################.  poincaré plot #####################

# pplot = PoincarePlot.with_linspace(section, top_o_coord, x_point1_coord,40)
# pplot.compute(400)
# np.save('./script/jellyfisch/75979_12/JF_75979_pp_hits_P.npy', pplot._hits)

################### tomographic plot #####################

#tpc= tomo(f"{patch_mat_file}",ax2,emi_vmin=0, emi_vmax=2.25e19)


##########################


Hits=np.load('./script/jellyfisch/75979_12/JF_75979_pp_hits_P.npy')
ax2.scatter(Hits[:,:, 0], Hits[:,:, 1], color="xkcd:dark grey", s=0.6, linewidths=0)
eps_s1 = 8.2e-07
eps_u1 = 7.4e-07


manifs = {name: {name} for name in mf_list}

for i in mf_list:
     manifs[i] = Manifold.load(f"{repository_path}{file_manifolds}{i}.pkl")
    # manifs[i].plot(ax=ax2, markersize=0, lw=0.7,labels=[None,None] if i!=0 else ['stable MF','unstable MF'])
    #  if i=='mf_1T':
       
    #    trajectory_s = manifs[i]._stable_trajectory
    #    trajectory_u = manifs[i]._unstable_trajectory

    #    mask_s = np.isfinite(trajectory_s[:,0]) | np.isfinite(trajectory_s[:,1])
    #    mask_u = np.isfinite(trajectory_u[:,0]) | np.isfinite(trajectory_u[:,1])

    #    M1_0_idx=len(manifs[i].clinics[0].trajectory)//2

    #    M1_0 = manifs[i].clinics[0].trajectory[M1_0_idx]

    #    trajectory_s=trajectory_s[mask_s]
    #    trajectory_u=trajectory_u[mask_u]
       
    #    idx_1s = IndexClosestPoint(trajectory_s, M1_0)
    #    idx_1u = IndexClosestPoint(trajectory_u , M1_0)

    #    trajectory_u= trajectory_u[:idx_1u]  
    #    trajectory_s= trajectory_s[:idx_1s]  

    #    idx_1s = IndexClosestPoint(trajectory_s, M1_0)
    #    idx_1u = IndexClosestPoint(trajectory_u, M1_0)

    #    trajectory_u= trajectory_u[:idx_1u]  
    #    trajectory_s= trajectory_s[:idx_1s] 

    #    trajectory_u= trajectory_u[::-1]      

    #    contour_points = np.vstack([
    #      manifs[i].rfp_s,      
    #      trajectory_s,              
    #      trajectory_u,      
    #     ])       

    #    Aire=manifs[i]._AdL_integral_points(contour_points)
       #ax2.plot(contour_points[:,0], contour_points[:,1], 'xkcd:dark orange', lw=1.5)
     manifs[i].plot(stepsize_limit=0.05, ax=ax2, markersize=0, lw=0.5,colors=["xkcd:royal blue", "xkcd:red"],labels=['stable MF','unstable MF'] if i=='mf_1T_e' else [None,None])


top_o.plot(ax=ax2, marker='o',s=40, color="xkcd:dark blue",label=None, zorder=10)
x_point1.plot(ax=ax2, marker='x', s=50, color="xkcd:dark blue",label=None, zorder=10)
x_point2.plot(ax=ax2, marker='x', s=50, color="xkcd:dark blue",label=None, zorder=10)
x_point3.plot(ax=ax2, marker='x', s=50, color="xkcd:dark blue",label=None, zorder=10)
x_point4.plot(ax=ax2, marker='x', s=50, color="xkcd:dark blue",label=None, zorder=10)
o_point_island.plot(ax=ax2, marker='o',s=40, color="xkcd:purple",label='New O-point', zorder=10)
x_point_island_n.plot(ax=ax2, marker='x', s=50, color="xkcd:purple",label='New X-point', zorder=10)


# top_o.plot(ax=ax2, marker='o', s=50, color="xkcd:green",label='O-point', zorder=10)
# x_point1.plot(ax=ax2, marker='x', s=70, color="xkcd:green",label='X-points', zorder=10)
# x_point2.plot(ax=ax2, marker='x', s=50, color="xkcd:green",label=None, zorder=10)
# x_point3.plot(ax=ax2, marker='x', s=50, color="xkcd:green",label=None, zorder=10)
# x_point4.plot(ax=ax2, marker='x', s=50, color="xkcd:green",label=None, zorder=10)
# ratio=2.8158



ax2.set_xlim(0.62, 1.15)
# ax2.set_ylim(-0.75, 0.75)
ax2.set_ylim(0, 0.4)

ax2.set_xlabel(r"$R[m]$")
ax2.set_ylabel(r"$Z[m]$")
ax2.set_aspect('equal') 



#ax2.set_title('Poincaré perturbed')

ax2.legend(loc='upper right', bbox_to_anchor=(1., 1), ncol=1, fontsize=9)

ax2.yaxis.set_major_locator(MaxNLocator(nbins=6, prune='both'))  # change nbins as needed


plt.savefig(f"{repository_path}figures/JF_75979_croissant_v2.png",  bbox_inches="tight",dpi=720, pad_inches=0, facecolor=fig.get_facecolor())
plt.show()
print(Aire)




