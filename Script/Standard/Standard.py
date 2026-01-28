import os, sys
proj_root = os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..'))  # remonte jusqu'à "Travail master"
if proj_root not in sys.path:
    sys.path.insert(0, proj_root)
from script.function.field_utils import plot_tomographic as tomo
from script.function.field_utils import plot_vessel as vessel
from matplotlib.patches import FancyArrowPatch
from pyoculus.fields import AxisymmetricCylindricalGridField
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
import numpy as np
from matplotlib import pyplot as plt
import logging
from scipy.io import loadmat
from matplotlib.tri import Triangulation
from matplotlib.patches import FancyArrowPatch




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
fig,ax = plt.subplots(1,1, figsize=(5,8))
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

####top fp top mf#########

# manifold_1T = Manifold(section, x_point1, x_point1,-x_point1_coord+top_o_coord, -x_point1_coord+top_o_coord)
# manifold_1T.compute(eps_s=9e-8, eps_u=8e-8, nint_s=6, nint_u=6, neps_s=80, neps_u=80) #8 14 240
# manifold_1T.save(f"{repository_path}{file_manifolds}mf_1T_NP_dir.pkl")

# # #####top fp bottom mf#########

# manifold_1B = Manifold(section, x_point1, x_point1,x_point1_coord-top_o_coord, x_point1_coord-top_o_coord)
# manifold_1B.compute(eps_s=9e-6, eps_u=8e-6, nint_s=4, nint_u=4, neps_s=80, neps_u=80)
# manifold_1B.save(f"{repository_path}{file_manifolds}mf_1B_NP_dir.pkl")

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

manifold_1T  = Manifold.load(f"{repository_path}{file_manifolds}mf_1T_NP_dir.pkl") 
manifold_1B  = Manifold.load(f"{repository_path}{file_manifolds}mf_1B_NP_dir.pkl")

Hits=np.load(f"{repository_path}ST_80064_pp_NP_hits.npy")
ax.scatter(Hits[:,:, 0], Hits[:,:, 1], color="xkcd:dark grey", s=1., linewidths=0,zorder=10)

#################

# manifold_1B.plot(stepsize_limit=0.3, ax=ax, markersize=0, lw=2,labels=[None, None])
# top_o.plot(ax=ax, marker='o',s=70, color="xkcd:dark blue",label='O-points', zorder=10)
# x_point1.plot(ax=ax, marker='x', s=70, color="xkcd:dark blue",label='X-points', zorder=10)






##### figure oral

# manifold_1T.plot(stepsize_limit=0.3,which='stable', ax=ax, markersize=0, lw=2,colors=['blue','red'], labels=[None, None])
# manifold_1B.plot(stepsize_limit=0.3, which='stable', ax=ax, markersize=0, lw=2,colors=['blue','red'], labels=[None, None])

# arrow = FancyArrowPatch(
#     posA=(0.795,-0.11),posB=(0.755,-0.127),
#     arrowstyle='-|>',               # simple head
#     mutation_scale=25,              # head size (adjust)
#     linewidth=0,                    # hide shaft line
#     color='blue',
#     zorder=50,
#     shrinkA=0, shrinkB=0
# )
# ax.add_patch(arrow)

# arrow2 = FancyArrowPatch(
#     posA=(0.6567,-0.1365),posB=(0.6667,-0.1358),
#     arrowstyle='-|>',               # simple head
#     mutation_scale=25,              # head size (adjust)
#     linewidth=0,                    # hide shaft line
#     color='blue',
#     zorder=50,
#     shrinkA=0, shrinkB=0
# )
# ax.add_patch(arrow2)

# ##### figure oral unstable

# manifold_1T.plot(stepsize_limit=0.3,which='unstable', ax=ax, markersize=0, lw=2,colors=['blue','red'], labels=[None, None])
# manifold_1B.plot(stepsize_limit=0.3, which='unstable', ax=ax, markersize=0, lw=2,colors=['blue','red'], labels=[None, None])

# arrow = FancyArrowPatch(
#     posA=(0.6991,-0.0655),posB=(0.6967,-0.0516),
#     arrowstyle='-|>',               # simple head
#     mutation_scale=25,              # head size (adjust)
#     linewidth=0,                    # hide shaft line
#     color='red',
#     zorder=50,
#     shrinkA=0, shrinkB=0
# )
# ax.add_patch(arrow)

# arrow2 = FancyArrowPatch(
#     posA=(0.7184,-0.206),posB=(0.7195,-0.2209),
#     arrowstyle='-|>',               # simple head
#     mutation_scale=25,              # head size (adjust)
#     linewidth=0,                    # hide shaft line
#     color='red',
#     zorder=50,
#     shrinkA=0, shrinkB=0
# )
# ax.add_patch(arrow2)




top_o.plot(ax=ax, marker='o',s=160,facecolor='cyan',edgecolor='black',linewidth=1.5,label='O-points', zorder=10)
# # #x_point1.plot(ax=ax, marker='+', s=220, facecolor='cyan', edgecolor='black', linewidth=1.5, label='X-points', zorder=10)


# ax.scatter([x_point1_coord[0]], [x_point1_coord[1]],
#            marker='x', s=180, color='black', linewidths=5, zorder=11)
# # 2) smaller colored cross on top (interior color)
# ax.scatter([x_point1_coord[0]], [x_point1_coord[1]],
#            marker='x', s=150, color='cyan', linewidths=3, zorder=12, label='X-points')



#####figure
vessel(ax)


# ...existing code...
# ellipse parameters (ajoute `alpha` pour l'orientation en radians)
cx, cy = top_o_coord[0]-0.0015, top_o_coord[1]  # centre
a, b = 0.078, 0.102    # semi-major (a) et semi-minor (b)
alpha = np.deg2rad(5.)  # orientation de l'ellipse en degrés -> radians (ex: 30°)

# parametric ellipse non-rotée
t = np.linspace(0, 2*np.pi, 400)
xx = a * np.cos(t)
yy = b * np.sin(t)

# rotation
ca, sa = np.cos(alpha), np.sin(alpha)
x = cx + ca * xx - sa * yy
y = cy + sa * xx + ca * yy

# choix du point tangent (angle paramétrique)
theta = -0.2
px0 = a * np.cos(theta)
py0 = b * np.sin(theta)
# appliquer la même rotation pour obtenir le point réel
px = cx + ca * px0 - sa * py0
py = cy + sa * px0 + ca * py0

# dérivée paramétrique avant rotation
dx_dt0 = -a * np.sin(theta)
dy_dt0 =  b * np.cos(theta)
# appliquer rotation à la dérivée
dx_dt = ca * dx_dt0 - sa * dy_dt0
dy_dt = sa * dx_dt0 + ca * dy_dt0

v = np.array([dx_dt, dy_dt], dtype=float)
v_norm = np.linalg.norm(v)
v_unit = v / v_norm if v_norm != 0 else np.array([1.0, 0.0])

# tracés
ax.plot(x, y, color='blue', lw=2)                       # ellipse orientée
# ax.scatter([cx], [cy],
#            s=220,
#            facecolor='cyan',
#            edgecolor='black',
#            linewidth=1.5,
#            zorder=30)

# flèche (seulement la pointe)
small_len = 0.03 * max(a, b)
posB = (px - v_unit[0] * small_len, py - v_unit[1] * small_len)
posA = (px, py)
arrow = FancyArrowPatch(posA=posA, posB=posB, arrowstyle='-|>', mutation_scale=25,
                        linewidth=0, color='blue', zorder=50, shrinkA=0, shrinkB=0)
ax.add_patch(arrow)




#ax.legend(loc='center right', bbox_to_anchor=(1.02, 0.33), ncol=1, fontsize=7)

ax.set_xlim(0.62, 1.15)
ax.set_ylim(-0.75, 0.75)

ax.set_xlabel(r"$R[m]$")
ax.set_ylabel(r"$Z[m]$")
ax.set_aspect('equal') 
# ax.set_title('(80064 / 1.1s.)')
plt.savefig(f"{repository_path}/figures/oral_presentation/oral_ST_dir_stable.png",  dpi=720,bbox_inches="tight", pad_inches=0, facecolor=fig.get_facecolor())
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



