import os, sys
proj_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))  # remonte jusqu'à "Travail master"
if proj_root not in sys.path:
    sys.path.insert(0, proj_root)
from matplotlib.ticker import FormatStrFormatter, FuncFormatter, MaxNLocator
from scipy.signal import argrelextrema
from pyoculus.fields import AxisymmetricCylindricalGridField
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
from scipy.signal import find_peaks
import numpy as np
from matplotlib import pyplot as plt
import logging
from scipy.io import loadmat
from matplotlib.tri import Triangulation , LinearTriInterpolator
from script.function.field_utils import add_constant_Bz as add_Bz
import scipy.io
from script.function.line_intersect import line_curve_MF_intersections as MF_inter
from script.function.line_intersect import sample_emi_along_line
from script.function.line_intersect import line_curve_intersections as line_intersections
from script.function.field_utils import plot_tomographic as tomo


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5), gridspec_kw={'width_ratios': [1, 1]})

plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.size": 9,
    }
)


logging.basicConfig(level=logging.DEBUG)

repository_path = './script/jellyfisch/77062_12/normal/'
pert_mat_file = 'JF_77062_120.mat'
patch_mat_file = './script/jellyfisch/77062_12/JF_77062_patch_120_C3.mat'


data = scipy.io.loadmat(f"{repository_path}{pert_mat_file}")

R=data['rr']

B_R_add,B_Z_add,psi_add= add_Bz(R, Bz_const=1.7e-2)  # Adding a constant Bz field of 17 mT

JFField = AxisymmetricCylindricalGridField.from_matlab_file_with_added_B(f"{repository_path}{pert_mat_file}", B_R_add, B_Z_add, psi_add, with_perturbation=True)

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


tpc,tri, emi= tomo(f"{patch_mat_file}",ax2,emi_vmin=0, emi_vmax=8e20)
tpc,tri, emi= tomo(f"{patch_mat_file}",ax1,emi_vmin=0, emi_vmax=8e20)


###########################. fixed point  #########################

####top fp top mf#########

# manifold_1T = Manifold(section, x_point1, x_point1,-x_point1_coord+top_o_coord, -x_point1_coord+top_o_coord)
# manifold_1T.compute(
#    eps_s=9e-6, eps_u=8e-6, nint_s=17, nint_u=16, neps_s=80, neps_u=240) #8 14 240

# manifold_1T.save('./script/current_divertor/manifolds_P/mf_1T_add_Bz.pkl')

#####top fp bottom mf#########

# manifold_1B = Manifold(section, x_point1, x_point1)
# manifold_1B.compute(
#        eps_s=9e-6, eps_u=8e-6, nint_s=10, nint_u=15, neps_s=80, neps_u=80)
# manifold_1B.save('./script/current_divertor/manifolds_P/mf_1B_add_Bz.pkl')

# #### bottom fp top mf

# manifold_2T = Manifold(section, x_point2, x_point2,-x_point2_coord+top_o_coord, -x_point2_coord+top_o_coord)
# manifold_2T.compute(
#       eps_s=9e-6, eps_u=8e-6, nint_s=32, nint_u=30, neps_s=80, neps_u=80)
# manifold_2T.save('./script/current_divertor/manifolds_P/mf_2T_add_Bz.pkl')

# #### bottom fp bottom mf

# manifold_2B = Manifold(section, x_point2, x_point2,x_point2_coord-top_o_coord, x_point2_coord-top_o_coord)
# manifold_2B.compute(
#       eps_s=9e-6, eps_u=8e-6, nint_s=32, nint_u=32, neps_s=80, neps_u=80)
# manifold_2B.save('./script/current_divertor/manifolds_P/mf_2B_add_Bz.pkl')

# #### right fp left mf

# manifold_3L = Manifold(section, x_point3, x_point3,-x_point3_coord+top_o_coord, -x_point3_coord+top_o_coord)
# manifold_3L.compute(
#       eps_s=9e-6, eps_u=8e-6, nint_s=33, nint_u=30, neps_s=80, neps_u=80)
# manifold_3L.save('./script/current_divertor/manifolds_P/mf_3L_add_Bz.pkl')

# # ##### right fp right mf

# manifold_3R = Manifold(section, x_point3, x_point3,x_point3_coord-top_o_coord, x_point3_coord-top_o_coord)
# manifold_3R.compute(
#       eps_s=9e-6, eps_u=8e-6, nint_s=20, nint_u=20, neps_s=80, neps_u=80)
# manifold_3R.save('./script/current_divertor/manifolds_P/mf_3R_add_Bz.pkl')

# ## top fp top mf 

# manifold_4T = Manifold(section, x_point4, x_point4,x_point4_coord-top_o_coord, x_point4_coord-top_o_coord)
# manifold_4T.compute(
#      eps_s=9e-6, eps_u=8e-6, nint_s=8, nint_u=8, neps_s=80, neps_u=80)
# manifold_4T.save('./script/current_divertor/manifolds_P/mf_4T_add_Bz.pkl')

# # # # ##### over the top fp bottom mf

# manifold_4B = Manifold(section, x_point4, x_point4,-x_point4_coord+top_o_coord, -x_point4_coord+top_o_coord)
# manifold_4B.compute(
#      eps_s=9e-6, eps_u=8e-6, nint_s=8, nint_u=30, neps_s=80, neps_u=80)
# manifold_4B.save('./script/current_divertor/manifolds_P/mf_4B_add_Bz.pkl')

##### loading and plotting  ###

manifold_1T  = Manifold.load('./script/current_divertor/manifolds_P/mf_1T_add_Bz.pkl')
manifold_1B  = Manifold.load('./script/current_divertor/manifolds_P/mf_1B_add_Bz.pkl')
manifold_2T  = Manifold.load('./script/current_divertor/manifolds_P/mf_2T_add_Bz.pkl')
manifold_2B  = Manifold.load('./script/current_divertor/manifolds_P/mf_2B_add_Bz.pkl')
manifold_3L  = Manifold.load('./script/current_divertor/manifolds_P/mf_3L_add_Bz.pkl')
manifold_3R  = Manifold.load('./script/current_divertor/manifolds_P/mf_3R_add_Bz.pkl')
manifold_4T  = Manifold.load('./script/current_divertor/manifolds_P/mf_4T_add_Bz.pkl')
manifold_4B  = Manifold.load('./script/current_divertor/manifolds_P/mf_4B_add_Bz.pkl')

manifold_1T.plot(stepsize_limit=0.1, ax=ax1, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"],labels=['stable MF','unstable MF'])
manifold_1B.plot(stepsize_limit=0.3,ax=ax1, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"],labels=[None,None])
manifold_2T.plot(stepsize_limit=0.1, ax=ax1, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"],labels=[None,None])
manifold_2B.plot(stepsize_limit=0.3, ax=ax1, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"],labels=[None,None])
manifold_3L.plot(stepsize_limit=0.1,ax=ax1, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"],labels=[None,None])
manifold_3R.plot(ax=ax1, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"],labels=[None,None])
manifold_4T.plot(ax=ax1, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"],labels=[None,None])
manifold_4B.plot(ax=ax1, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"],labels=[None,None])

top_o.plot(ax=ax1, marker='o', color="xkcd:white",label='O-point')
x_point1.plot(ax=ax1, marker='x', color="xkcd:white",label='X-points')
x_point2.plot(ax=ax1, marker='x', color="xkcd:white",label=None)
x_point3.plot(ax=ax1, marker='x', color="xkcd:white",label=None)
x_point4.plot(ax=ax1, marker='x', color="xkcd:white",label=None)
ratio=2.8158

p0 = np.array([0.875, -0.32])
p1 = np.array([0.875, -0.5])  #x_point2.coords[0]

pts, emi_vals,dist,perio_MANTIS,std_perio_MANTIS=sample_emi_along_line(p0, p1, tri, emi, n_samples=400,sig=10,ax=ax2)

intersections_1, dist1,perio_MF,std_perio_MF=MF_inter(p0, p1, manifold_1T, ax1=ax1,ax2=ax2)

#####figure


ax1.set_xlim(0.62, 1.15)
ax1.set_ylim(-0.75, 0.75)
ax1.set_xlim(0.7, 1.0)
ax1.set_ylim(-0.6, -0.1)
ax1.set_xlabel(r"$R[m]$")
ax1.set_ylabel(r"$Z[m]$")
ax1.set_aspect('equal') 
ax1.set_title('(77062 / 1.2s), added Bz=0.017 [T]')
# plt.savefig(f"{repository_path}figures/JF_77062_Bz_added_zoom.png", bbox_inches='tight', dpi=720)
# plt.show()


ax2.set_ylim(0, np.nanmax(emi_vals)*1.2)
ax2.set_xlim(min(dist), 0.19)
ax2.set_xlabel("distance along the line [m]")
ax2.set_ylabel("Relative Emissivity")


txt = (
    f"MANTIS (avg peak distance): {perio_MANTIS:.4f} ± {std_perio_MANTIS:.4f} [m]\n"
    f"MF (avg peak distance): {perio_MF:.4f} ± {std_perio_MF:.4f} [m]"
)

# ajouter en haut à droite de l'axe de profil
ax2.text(
    0.98, 0.98, txt,
    transform=ax2.transAxes,
    ha='right', va='top',
    fontsize=9,
    color='black',
    bbox=dict(facecolor='white', alpha=0.85, edgecolor='none')
)

ax2.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))

ax2.legend(loc='lower left')

plt.savefig(f"{repository_path}figures/jellyfisch_77062_added_Bz_120_periodicity.png", dpi=150, bbox_inches=None, facecolor=fig.get_facecolor())
plt.show()
print(f"Bz mismatch is {perio_MF/perio_MANTIS}")



#####zoomed





