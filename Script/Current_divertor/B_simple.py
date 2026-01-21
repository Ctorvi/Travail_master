from pyoculus.fields import AxisymmetricCylindricalGridField
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
from matplotlib.ticker import FormatStrFormatter
import os, sys
proj_root = os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..'))  # remonte jusqu'à "Travail master"
if proj_root not in sys.path:
    sys.path.insert(0, proj_root)
import numpy as np
from matplotlib import pyplot as plt
import logging
from scipy.io import loadmat
from matplotlib.tri import Triangulation
from scipy.special import ellipk, ellipe
from script.function.field_utils import B_loop_cyl, B_sum_of_loops, BR_BZ_PSI,plot_B_parallel_along_line, plot_B_parallel_along_line_from_loops
from matplotlib.gridspec import GridSpec

def resample_by_arclength(pts, n_points, closed=False):
    """
    Rééchantillonne une courbe pts (N,2) en n_points uniformément espacés en arc-length.
    Args:
      pts : array-like (N,2) ou (2,N) -> coordonnées (x,y) ou (R,Z)
      n_points : int, nombre de points de sortie
      closed : bool, si True la courbe est considérée fermée (dernier->premier inclus)
    Returns:
      ndarray (n_points, 2)
    """
    import numpy as np
    pts = np.asarray(pts)
    if pts.ndim == 2 and pts.shape[0] == 2 and pts.shape[1] > 2:
        pts = pts.T
    if pts.ndim != 2 or pts.shape[1] != 2:
        raise ValueError("pts must be shape (N,2) or (2,N)")

    if closed:
        # append first point to close curve for length calc, will be removed after interp
        pts_closed = np.vstack([pts, pts[0]])
        diffs = np.diff(pts_closed, axis=0)
    else:
        diffs = np.diff(pts, axis=0)

    seglen = np.hypot(diffs[:, 0], diffs[:, 1])
    s = np.concatenate(([0.0], np.cumsum(seglen)))
    total = s[-1]
    if total == 0 or n_points <= 1:
        return np.tile(pts[0], (max(1, n_points), 1))

    s_uniform = np.linspace(0.0, total, n_points)

    # separate interpolation of x and y
    x = pts[:, 0]
    y = pts[:, 1]
    x_u = np.interp(s_uniform, s, x)
    y_u = np.interp(s_uniform, s, y)

    out = np.column_stack((x_u, y_u))
    return out

fig, ax=plt.subplots(1, 1, figsize=(5, 8))

plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.size": 9,
    }
)

logging.basicConfig(level=logging.DEBUG)


repository_path = './script/jellyfisch/77021/1_20/backoff/'

pert_mat_file = 'JF_77021_120_BO78_OS.mat'


repository_path2 = './script/jellyfisch/77062_12/normal/'

pert_mat_file2 = 'JF_77062_120_BO78_OS.mat'

JFField2 = AxisymmetricCylindricalGridField.from_matlab_file(f"{repository_path2}{pert_mat_file2}", with_perturbation=True)

JFField = AxisymmetricCylindricalGridField.from_matlab_file(f"{repository_path}{pert_mat_file}", with_perturbation=True)

section = CylindricalBfieldSection(JFField,phi0=1.9,R0=0.88, Z0=0)

section2 = CylindricalBfieldSection(JFField2,phi0=1.9,R0=0.88, Z0=0)

#####MANIFOLD AND CURRENT CALCULATION#########

manifold_1B  = Manifold.load(f"{repository_path}manifolds_P/OS_BO78/mf_1B.pkl")
manifold_1T  = Manifold.load(f"{repository_path}manifolds_P/OS_BO78/mf_1T.pkl")


manifold_1B_2  = Manifold.load(f"{repository_path2}manifolds_P/OS_BO78/mf_1B.pkl")

x_point1=FixedPoint(section)
x_point1= manifold_1B.fixedpoint_1
x_p= x_point1.coords[0]

manifold_1T.plot(ax=ax, stepsize_limit=0.1, linewidth=0.7, markersize=0)



T1=manifold_1B.unstable
T2=manifold_1B_2.unstable

valid_rows = np.isfinite(T1).all(axis=1)
T1 = T1[valid_rows]
T1 = T1[:-42,:]

valid_rows = np.isfinite(T2).all(axis=1)
T2 = T2[valid_rows]

niteration1=len(T1)
niteration2=len(T2)






length=np.sum(np.sqrt(np.diff(T1[:,0])**2 + np.diff(T1[:,1])**2 ))

#Aire_traj=0.03*length 

N=10000

R_Z = resample_by_arclength(T1, n_points=N, closed=False)

####### current 1 ########
Aire_traj=0.01*length #low est
j_par=-5*(10**4) #A/m2 low est.
I_surface=j_par*Aire_traj  #A
I_loop_low=I_surface/N #A

#### current 2 ######

Aire_traj_high=0.03*length #high est
j_par_high=-25*(10**4) #A/m2
I_surface_high=j_par_high*Aire_traj_high  #A
I_loop_high=I_surface_high/N  #A

####

r_p=0.84

z_p=-0.4

B_R_low_tot, B_R_high_tot = 0.0, 0.0
B_Z_low_tot, B_Z_high_tot = 0.0, 0.0

for i in range(N):

    B_R1, B_Z1, B_phi1 = B_loop_cyl(r_p, z_p, R_Z[i, 0], R_Z[i, 1], I_loop_low)
    B_R_low_tot +=B_R1
    B_Z_low_tot +=   B_Z1

    B_R2, B_Z2, B_phi2 = B_loop_cyl(r_p, z_p, R_Z[i, 0], R_Z[i, 1], I_loop_high)
    B_R_high_tot +=  B_R2
    B_Z_high_tot +=  B_Z2


ax.plot(R_Z[:, 0], R_Z[:, 1], 'x', color='green',markersize=4, markeredgewidth=0.6, label='Manifold 1B' )
plt.show()
print("B low (T) at (R=0.84m, Z=-0.4m): ", {B_R_low_tot, B_Z_low_tot})
print("B high (T) at (R=0.84m, Z=-0.4m): ", {B_R_high_tot, B_Z_high_tot})
print("B field total (T) at (R=0.84m, Z=-0.4m): ", JFField.B([r_p, 1.9, z_p]))
