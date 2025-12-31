import os, sys
proj_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..'))  # remonte jusqu'à "Travail master"
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
from script.function.line_intersect import line_curve_MF_intersections as MF_inter
from script.function.line_intersect import sample_emi_along_line
from script.function.line_intersect import line_curve_intersections as line_intersections
from script.function.field_utils import plot_tomographic as tomo

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6, 4), gridspec_kw={'width_ratios': [1, 1]})

plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.size": 9,
    }
)
logging.basicConfig(level=logging.DEBUG)

repository_path = './script/jellyfisch/77021/1_20/backoff/'
pert_mat_file = 'JF_77021_120.mat'
patch_mat_file = './script/jellyfisch/77021/1_20/JF_patch_77021_120_HeI.mat'

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

# x_point3= FixedPoint(section)
# x_point3.find(1, [1.05,-0.58], method='scipy.root')
# x_point3_coord = x_point3.coords[0]

# x_point4= FixedPoint(section)
# x_point4.find(1, [0.75,0.65], method='scipy.root')
# x_point4_coord = x_point4.coords[0]

###########################  tomographic reconstruction  #########################


tpc,tri,emi= tomo(f"{patch_mat_file}",ax=ax1,emi_vmin=0, emi_vmax=2.e19)



###########################. computation #########################



manifold_1T  = Manifold.load(f"{repository_path}manifolds_P/OS_BO78/mf_1T_n20.pkl")
manifold_1B  = Manifold.load(f"{repository_path}manifolds_P/OS_BO78/mf_1B.pkl")
manifold_2T  = Manifold.load(f"{repository_path}manifolds_P/OS_BO78/mf_2T.pkl")
manifold_2B  = Manifold.load(f"{repository_path}manifolds_P/OS_BO78/mf_2B.pkl")


manifold_1T.plot(stepsize_limit=0.2,ax=ax1, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"])
manifold_1B.plot(ax=ax1, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"])
manifold_2T.plot(stepsize_limit=0.3, ax=ax1, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"])
manifold_2B.plot(stepsize_limit=0.3, ax=ax1, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"])


top_o.plot(ax=ax1, marker='o', color="xkcd:white")
x_point1.plot(ax=ax1, marker='x', color="xkcd:white")
x_point2.plot(ax=ax1, marker='x', color="xkcd:white")

ratio=2.8158

# p0=np.array([0.85, -0.31])
# p1=x_point2.coords[0]

# intersections_1=MF_inter(p0, p1, manifold_1T,ax=ax)
#####figure


ax1.set_xlim(0.7, 1.0)
ax1.set_ylim(-0.58, -0.08)
ax1.set_xlabel(r"$R[m]$")
ax1.set_ylabel(r"$Z[m]$")
ax1.set_aspect('equal')
ax1.set_title('(77021 / 1.20 s / 1.9 rad)')
#plt.savefig(f"{repository_path}/figures/Jellyfisch_Pert_mf_patch_line.png", bbox_inches='tight', dpi=720)
# plt.show()


#####zoomed

# p0 = np.array([0.842, -0.29])
# p1 = np.array([0.842, -0.55])

p0 = np.array([0.85, -0.3])
p1 = np.array([0.822, -0.52])


pts, emi_vals,dist,perio_MANTIS,std_perio_MANTIS= sample_emi_along_line(p0, p1, tri, emi, n_samples=400,sig=8,ax=ax2)

intersections_1, dist1,perio_MF,std_perio_MF=MF_inter(p0, p1, manifold_1T,ax1=ax1,ax2=ax2)

ax2.set_ylim(0, np.nanmax(emi_vals)*1.2)
ax2.set_xlim(min(dist), 0.21)
ax2.set_xlabel("distance along the line [m]")
ax2.set_ylabel("Relative Emissivity")


txt = (
    f"MANTIS (avg $\Delta d$): {perio_MANTIS:.4f} ± {std_perio_MANTIS:.4f} [m]\n"
    f"MF (avg $\Delta d$): {perio_MF:.4f} ± {std_perio_MF:.4f} [m]"
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
ax2.legend(loc='lower left', bbox_to_anchor=(-0.1, 0.01), ncol=1, fontsize=9)


plt.savefig(f"{repository_path}figures/Jellyfisch_77021_120_periodicity.png", dpi=150, bbox_inches=None, facecolor=fig.get_facecolor())
plt.show()
print(f"Bz mismatch is {perio_MF/perio_MANTIS}")

