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


fig = plt.figure(figsize=(8, 4), constrained_layout=False)
gs = GridSpec(2, 2, figure=fig, width_ratios=[1, 1], height_ratios=[2, 1], wspace=0.3, hspace=0.3)
ax = fig.add_subplot(gs[:, 0])        # left column spanning both rows
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[1, 1])


plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.size": 9,
    }
)

logging.basicConfig(level=logging.DEBUG)

repository_path = './script/jellyfisch/77062_12/normal/'
repository_path2 = './script/jellyfisch/77021/1_20/backoff/'
repository_path3 = './script/jellyfish/75979_12/'

pert_mat_file = 'JF_77062_120_BO78_OS.mat'

pert_mat_file2 = 'JF_77021_120_BO78_OS.mat'

pert_mat_file3 = 'mat_files/JF_75979_120_BO78_OS.mat'


JFField = AxisymmetricCylindricalGridField.from_matlab_file(f"{repository_path}{pert_mat_file}", with_perturbation=True)

section = CylindricalBfieldSection(JFField,phi0=1.9,R0=0.88, Z0=0)


#####MANIFOLD AND CURRENT CALCULATION#########

manifold_1B  = Manifold.load(f"{repository_path}manifolds_P/OS_BO78/mf_1B.pkl")
manifold_1T  = Manifold.load(f"{repository_path}manifolds_P/OS_BO78/mf_1T.pkl")

x_point1=FixedPoint(section)
x_point1= manifold_1B.fixedpoint_1
x_p= x_point1.coords[0]

manifold_1T.plot(ax=ax, stepsize_limit=0.1, linewidth=0.7, markersize=0)

T1=manifold_1B.unstable
T1=T1[:-42,:]
T=T1[::4,:]

length=np.sum(np.sqrt(np.diff(T[:,0])**2 + np.diff(T[:,1])**2 ))

#Aire_traj=0.78*length

#Aire_traj=0.03*length 

####### current 1 ########
Aire_traj=0.01*length #low est
#j_par=-25*(10**4) #A/m2
j_par=-5*(10**4) #A/m2 low est.
I_surface=j_par*Aire_traj  #A
I_loop_low=I_surface/np.size(T,0)  #A


#### current 2 ######

Aire_traj_high=0.03*length #high est
j_par_high=-25*(10**4) #A/m2
I_surface_high=j_par_high*Aire_traj_high  #A
I_loop_high=I_surface_high/np.size(T,0)  #A

####


a=T[:,0]
Z0=T[:,1]

R_grid,Z_grid=np.meshgrid(JFField.R, JFField.Z)
B_R, B_Z, B_phi = B_sum_of_loops(R_grid, Z_grid, a, Z0, I_loop_low)

B_R_tot, B_Z_tot, psi = BR_BZ_PSI(R_grid, Z_grid, B_R, B_Z, x_p)

ax.set_ylim(-0.6, 0.00)
ax.quiver(R_grid, Z_grid, B_R_tot, B_Z_tot)
ax.set_xlabel("R[m]")
ax.set_ylabel("Z[m]")
ax.set_title("Magnetic field from divertor leg current")
ax.set_aspect("equal")
ax.plot(T[:,0], T[:,1], 'x', color='green',markersize=4, markeredgewidth=0.6, label='Manifold 1B' )


###### 

p0 = np.array([0.87, -0.33])
p1 = np.array([0.87, -0.5])  #x_point2.coords[0]


s1,B_field=plot_B_parallel_along_line(JFField, p0, p1, n=200, phi=1.9, ax=ax2, plot=True, return_data=True)

s2, B_loop =plot_B_parallel_along_line_from_loops(a, Z0, I_loop_low, p0, p1, n=200, phi=1.9, ax=ax2, plot=True, ax_line=ax, return_data=True)

s3, B_loop_high =plot_B_parallel_along_line_from_loops(a, Z0, I_loop_high, p0, p1, n=200, phi=1.9, ax=None, plot=False, ax_line=None, return_data=True)



ax2.plot(s2, B_loop_high, '-', color='C2', lw=1, label=r'B$_{||}$ loops high est.')
ax2.set_ylim(-0.004, 0.016)
B_diff = (B_field + B_loop)/B_field
B_diff_high = (B_field + B_loop_high)/B_field
ax2.legend(fontsize=8)
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax3.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))


ax.legend()
ax3.plot(s1, B_diff, '-', color='C1', lw=1,label='low est.')
ax3.plot(s1, B_diff_high, '-', color='C2', lw=1 ,label='high est.')
ax3.set_ylim(0.75, 1.00)
ax3.set_xlabel("arc length s [m]")
ax3.set_ylabel(r"$B_{\parallel corr} / B_{\parallel field}$", fontsize=14)
plt.savefig("./script/current_divertor/Div_current_both_estimation.png", dpi=300,bbox_inches="tight", pad_inches=0)
#plt.savefig("./script/current_divertor/Div_current_high_estimation.png", dpi=300,bbox_inches="tight", pad_inches=0)
plt.show()



#B_R_tot, B_Z_tot, psi = zero_near_curve(B_R_tot, B_Z_tot, psi, R_grid, Z_grid, T1, r_cut=0.01)
#B_R_tot, B_Z_tot, psi = zero_above_curve(B_R_tot, B_Z_tot, psi, R_grid, Z_grid, T1, mode='interp', fill_outside=False)
# np.savez("./Script/Current_divertor/Current_div_B.npz", B_R_tot=B_R_tot, B_Z_tot=B_Z_tot, psi=psi)

