from pyoculus.fields import AxisymmetricCylindricalGridField
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
import numpy as np
from matplotlib import pyplot as plt
import logging
from scipy.io import loadmat
from matplotlib.tri import Triangulation
from scipy.special import ellipk, ellipe
from Script.function.field_utils import B_loop_cyl, B_sum_of_loops, BR_BZ_PSI,zero_near_curve




fig, ax = plt.subplots(1, 1, figsize=(6, 6.6))


plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.size": 9,
    }
)


logging.basicConfig(level=logging.DEBUG)

repository_path = './Script/Jellyfisch/77062_12/normal/'
pert_mat_file = 'JF_77062_12.mat'


JFField = AxisymmetricCylindricalGridField.from_matlab_file(f"{repository_path}{pert_mat_file}", with_perturbation=True)

section = CylindricalBfieldSection(JFField,phi0=1.9,R0=0.88, Z0=0)


#####MANIFOLD AND CURRENT CALCULATION#########

manifold_1B  = Manifold.load(f"{repository_path}manifolds_P/mf_1B.pkl")
manifold_1T  = Manifold.load(f"{repository_path}manifolds_P/mf_1T.pkl")

x_point1=FixedPoint(section)
x_point1= manifold_1B.fixedpoint_1
x_p= x_point1.coords[0]

manifold_1T.plot(ax=ax, stepsize_limit=0.1, linewidth=0.7, markersize=0)

T1=manifold_1B.unstable
T1=T1[:-42,:]
T=T1[::4,:]


length=np.sum(np.sqrt(np.diff(T[:,0])**2 + np.diff(T[:,1])**2 ))


Aire_traj=0.78*length



####### current 1 ########


j_par=-7*10e1 #A/m2
I_surface=j_par*Aire_traj  #A
I_loop=I_surface/np.size(T,0)  #A

a=T[:,0]
Z0=T[:,1]

R_grid,Z_grid=np.meshgrid(JFField.R, JFField.Z)

B_R, B_Z, B_phi = B_sum_of_loops(R_grid, Z_grid, a, Z0, I_loop)

B_R_tot, B_Z_tot, psi = BR_BZ_PSI(R_grid, Z_grid, B_R, B_Z, x_p)

B_R_tot, B_Z_tot, psi = zero_near_curve(B_R_tot, B_Z_tot, psi, R_grid, Z_grid, T1, r_cut=0.01)


np.savez("./Script/Current_divertor/Current_div_B.npz", B_R_tot=B_R_tot, B_Z_tot=B_Z_tot, psi=psi)


ax.set_ylim(-0.6, 0.00)
ax.quiver(R_grid, Z_grid, B_R_tot, B_Z_tot)
ax.set_xlabel("R[m]")
ax.set_ylabel("Z[m]")
ax.set_title("Champ vectoriel 2D from current Separatrix-Strike-point")
ax.set_aspect("equal")
ax.plot(T[:,0], T[:,1], '+', color='blue', linewidth=0.7, label='Manifold 1B' )
ax.scatter(R_grid[10,10],Z_grid[10,10], color='red', marker='x')

BR = B_R_tot[10,10]
BZ = B_Z_tot[10,10]

txt = (
    f"B= ({BR:.3e}, {BZ:.3e}) T\n"
)

# positionner la boîte (coordonnées de données) — ajuste si besoin
ax.annotate(
    txt,
    xy=(R_grid[10,10], Z_grid[10,10]),   # position de la flèche/point
    xytext=(0, -12),                      # décalage (x,y) en "points"
    textcoords='offset points',
    ha='center',
    va='top',
    fontsize=9,
    color='black',
    bbox=dict(facecolor='white', alpha=1, edgecolor='black')
)


plt.savefig("./Script/Current_divertor/B_field_from_current_separatrix_strike_point.png", dpi=300)
plt.show()







