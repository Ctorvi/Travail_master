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



type_fp=['o', 'x', 'x', 'x', 'x']

emi_vals_max=[3.e20, 8e20, 2.25e19]

first_guess=[[0.9,0.18], [0.8,-0.27], [0.8,-0.57], [1.05,-0.58], [0.75,0.65]]

name_mf=['mf_1T.pkl', 'mf_1B.pkl', 'mf_2T.pkl', 'mf_2B.pkl', 'mf_3L.pkl', 'mf_3R.pkl', 'mf_4T.pkl', 'mf_4B.pkl']

titles=['(77021 / 1.2s)', '(77062 / 1.2s)', '(75979 / 1.2s)']

fig = plt.figure(figsize=(7, 8))

gs = GridSpec(2, 3, height_ratios=[1.5, 1])

ax1 = fig.add_subplot(gs[:1, 0])
ax1_d = fig.add_subplot(gs[1, 0])
ax2 = fig.add_subplot(gs[:1, 1])
ax2_d = fig.add_subplot(gs[1, 1])
ax3 = fig.add_subplot(gs[:1, 2])
ax3_d = fig.add_subplot(gs[1, 2])


list_ax = [ax1, ax1_d, ax2, ax2_d, ax3, ax3_d]

plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.size": 9,
    }
)

repository_path = ['./script/jellyfisch/77021/1_20/backoff/', './script/jellyfisch/77062_12/normal/','./script/jellyfisch/75979_12/']
pert_mat_file = ['JF_77021_120.mat', 'JF_77062_120.mat', 'JF_75979_120.mat']
patch_mat_file = ['./script/jellyfisch/77021/1_20/JF_patch_77021_120_C3.mat', './script/jellyfisch/77062_12/JF_77062_patch_120_C3.mat', './script/jellyfisch/75979_12/JF_75979_patch_120_He2.mat']

logging.basicConfig(level=logging.DEBUG)

for i in range(len(repository_path)):
    JFField = AxisymmetricCylindricalGridField.from_matlab_file(f"{repository_path[i]}{pert_mat_file[i]}", with_perturbation=True)

    section = CylindricalBfieldSection(JFField,phi0=1.9,R0=0.88, Z0=0)

    for j in range(len(type_fp)):
        fp = FixedPoint(section)
        fp.find(1, first_guess[j], method='scipy.root')
        fp.plot(ax=list_ax[2*i], marker='o' if 'o' in type_fp[j] else 'x', color="xkcd:white",zorder=10)
        fp.plot(ax=list_ax[2*i+1], marker='o' if 'o' in type_fp[j] else 'x', color="xkcd:white",zorder=10)

    tpc= tomo(f"{patch_mat_file[i]}",list_ax[2*i],emi_vmin=0, emi_vmax=emi_vals_max[i])
    tpc= tomo(f"{patch_mat_file[i]}",list_ax[2*i+1],emi_vmin=0, emi_vmax=emi_vals_max[i])

    for mf_name in name_mf:
        manifold = Manifold.load(f"{repository_path[i]}manifolds_P/{mf_name}")
        manifold.plot(stepsize_limit=0.2, ax=list_ax[2*i], markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"])
        manifold.plot(stepsize_limit=0.2, ax=list_ax[2*i+1], markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"])


    list_ax[2*i].set_xlim(0.62, 1.15)
    list_ax[2*i].set_ylim(-0.75, 0.75)

    list_ax[2*i].set_xlabel(r"$R[m]$")
    list_ax[2*i].set_aspect('equal') 
    list_ax[2*i].set_title(f"{titles[i]}")

    list_ax[2*i+1].set_xlim(0.7, 1.0)
    list_ax[2*i+1].set_ylim(-0.6, -0.1)
    list_ax[2*i+1].set_xlabel(r"$R[m]$")
    list_ax[2*i+1].set_aspect('equal')

list_ax[0].set_ylabel(r"$Z[m]$")
list_ax[1].set_ylabel(r"$Z[m]$")
list_ax[2].set_yticklabels([])
list_ax[3].set_yticklabels([])
list_ax[4].set_yticklabels([])
list_ax[5].set_yticklabels([])


plt.show()




