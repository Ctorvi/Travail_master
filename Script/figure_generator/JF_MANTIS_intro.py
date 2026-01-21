import os, sys
from pyoculus.fields import AxisymmetricCylindricalGridField
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
proj_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..','..'))  # remonte jusqu'à "Travail master"
if proj_root not in sys.path:
    sys.path.insert(0, proj_root)
import numpy as np
from matplotlib import pyplot as plt
import logging
from matplotlib.gridspec import GridSpec
from scipy.io import loadmat
from matplotlib.tri import Triangulation
from script.function.field_utils import plot_tomographic as tomo
from mpl_toolkits.axes_grid1.inset_locator import inset_axes



emi_vals_max=[3.5e20, 8e20, 2.25e19,6.e20]

subscripts = ['(a)', '(b)', '(c)']


titles=['(JF 77021 / 1.2 s / C-III)', '(JF 77062 / 1.2 s / C-III)', '(JF 75979 / 1.2 s / He-II)']

title_MANTIS = ['C-III', 'C-III',  'He-II']

fig = plt.figure(figsize=(7, 7))

gs = GridSpec(1, 3)

ax1 = fig.add_subplot(gs[0, 0])
ax1_d = fig.add_subplot(gs[0, 1])
ax1_dd = fig.add_subplot(gs[0, 2])


list_ax = [ax1, ax1_d, ax1_dd]

list_trace = [] 

plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.size": 8,
    }
)


########################. 77021, 77062,75979 ##########################

repository_path = ['./script/jellyfisch/77021/1_20/backoff/', './script/jellyfisch/77062_12/normal/',
                   './script/jellyfisch/75979_12/']


patch_mat_file = ['./script/jellyfisch/77021/1_20/JF_patch_77021_120_C3.mat', './script/jellyfisch/77062_12/JF_77062_patch_120_C3.mat',
                    './script/jellyfisch/75979_12/mat_files/JF_75979_patch_120_He2.mat']


logging.basicConfig(level=logging.DEBUG)


# tpc,tri,emi= tomo(f"{patch_mat_file[2]}",list_ax[2],emi_vmin=0, emi_vmax=emi_vals_max[2])


for i in range(len(repository_path)):

    #JFField = AxisymmetricCylindricalGridField.from_matlab_file(f"{repository_path[i]}{pert_mat_file[i]}", with_perturbation=True)

    #section = CylindricalBfieldSection(JFField,phi0=1.9,R0=0.88, Z0=0)

    #tpc,tri,emi= tomo(f"{patch_mat_file[i]}",list_ax[3*i],emi_vmin=0, emi_vmax=emi_vals_max[i])
    tpc,tri,emi= tomo(f"{patch_mat_file[i]}",list_ax[i],emi_vmin=0, emi_vmax=emi_vals_max[i])
  

    list_ax[i].set_xlim(0.62, 1.15)
    list_ax[i].set_ylim(-0.75, 0.75)

    list_ax[i].set_xlabel(r"$R[m]$")
    list_ax[i].set_aspect('equal') 
    list_ax[i].set_title(f"{titles[i]}")



for ax, lab in zip(list_ax, subscripts):
    ax.text(
        0.2, 0.95, lab,
        transform=ax.transAxes,
        fontsize=10,
        fontweight='bold',
        ha='right',
        va='top',
        color='xkcd:white',
        #bbox=dict(facecolor=bbox_face, edgecolor='none', pad=2, alpha=0.8),
        zorder=200
    )

list_ax[0].set_ylabel(r"$Z[m]$")
list_ax[1].set_yticklabels([])
list_ax[2].set_yticklabels([])




cax = inset_axes(list_ax[2], width="4%", height="80%", loc='center right',
                 bbox_to_anchor=(0.2, -0.1, 1, 1), bbox_transform=list_ax[2].transAxes, borderpad=0)


# create colorbar and hide the scientific offset text (e.g. "1e20")
cbar = fig.colorbar(tpc, cax=cax, label='Relative emissivity')
# hide the offset text that matplotlib places (the "1e20" above/beside the bar)
cbar.ax.yaxis.get_offset_text().set_visible(False)
# refresh ticks/formatter
cbar.update_ticks()
# ...existing code...

plt.savefig("./figures/MANTIS_MF_intro.png", dpi=750, bbox_inches="tight", pad_inches=0, facecolor=fig.get_facecolor())
plt.show()




