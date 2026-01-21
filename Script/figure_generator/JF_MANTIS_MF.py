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

type_fp=['o', 'x', 'x', 'x', 'x']

emi_vals_max=[3.5e20, 8e20, 2.25e19,6.e20]
#emi_vals_max=[3.e19, 5e18, 2.25e19,6.e20]


first_guess=[[0.9,0.18], [0.8,-0.27], [0.8,-0.57], [1.05,-0.58], [0.75,0.65]]

name_mf=['mf_1T.pkl', 'mf_1B.pkl', 'mf_2T.pkl', 'mf_2B.pkl', 'mf_3L.pkl', 'mf_3R.pkl', 'mf_4T.pkl', 'mf_4B.pkl']

allowed_mf_3 = {'mf_1T.pkl','mf_1B.pkl','mf_3L.pkl','mf_3R.pkl'}
allowed_fp_3 = [True, True, False, True, False]
label_fp = ['O-point', 'X-point', None, None, None]



titles=['(JF 77021 / 1.2 s)', '(JF 77062 / 1.2 s)', '(JF 75979 / 1.2 s)', '(JF 77060 /  1.2 s)', '(STD 80064 / 1.1 s)']

title_MANTIS = ['C-III', 'C-III',  'He-II', 'C-III', 'C-III']

#title_MANTIS = ['He-I', 'He-I',  'He-II', 'C-III', 'C-III']

fig = plt.figure(figsize=(9, 11))

gs = GridSpec(3, 5, height_ratios=[2.2, 1, 1])

ax1 = fig.add_subplot(gs[:1, 0])
ax1_d = fig.add_subplot(gs[1, 0])
ax1_dd = fig.add_subplot(gs[2, 0])

ax2 = fig.add_subplot(gs[:1, 1])
ax2_d = fig.add_subplot(gs[1, 1])
ax2_dd = fig.add_subplot(gs[2, 1])

ax3 = fig.add_subplot(gs[:1, 2])
ax3_d = fig.add_subplot(gs[1, 2])
ax3_dd = fig.add_subplot(gs[2, 2])

ax4 = fig.add_subplot(gs[:1, 3])
ax4_d = fig.add_subplot(gs[1, 3])
ax4_dd = fig.add_subplot(gs[2, 3])

ax5 = fig.add_subplot(gs[:1, 4])
ax5_d = fig.add_subplot(gs[1, 4])
ax5_dd = fig.add_subplot(gs[2, 4])

list_ax = [ax1, ax1_d, ax1_dd, ax2, ax2_d, ax2_dd, ax3, ax3_d, ax3_dd, ax4, ax4_d, ax4_dd, ax5, ax5_d, ax5_dd]

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
                   './script/jellyfisch/75979_12/', './script/jellyfisch/77060_12/']

pert_mat_file = ['JF_77021_120_BO78_OS.mat', 'JF_77062_120_BO78_OS.mat', 'mat_files/JF_75979_120_BO78_OS.mat', 'JF_77060_120_BO78_OS.mat']

patch_mat_file = ['./script/jellyfisch/77021/1_20/JF_patch_77021_120_C3.mat', './script/jellyfisch/77062_12/JF_77062_patch_120_C3.mat',
                     './script/jellyfisch/75979_12/mat_files/JF_75979_patch_120_He2.mat', './script/jellyfisch/77060_12/JF_patch_77060_120_C3.mat']

#patch_mat_file = ['./script/jellyfisch/77021/1_20/JF_patch_77021_120_HeI.mat', './script/jellyfisch/77062_12/JF_77062_patch_120_HeI_706.mat',
#                    './script/jellyfisch/75979_12/mat_files/JF_75979_patch_120_He2.mat', './script/jellyfisch/77060_12/JF_patch_77060_120_C3.mat']



pp_file =['JF_77021_pp_hits.npy', 'JF_77062_pp_hits.npy', 'JF_75979_pp_hits_P.npy', 'JF_77060_pp_hits.npy']

subscripts = ['(a1)', '(a2)', '(a3)', '(b1)', '(b2)', '(b3)', '(c1)', '(c2)', '(c3)', '(d1)', '(d2)', '(d3)', '(e1)', '(e2)', '(e3)']

logging.basicConfig(level=logging.DEBUG)

for i in range(len(repository_path)):

    JFField = AxisymmetricCylindricalGridField.from_matlab_file(f"{repository_path[i]}{pert_mat_file[i]}", with_perturbation=True)

    section = CylindricalBfieldSection(JFField,phi0=1.9,R0=0.88, Z0=0)

    Hits=np.load(f"{repository_path[i]}{pp_file[i]}")

    list_ax[3*i].scatter(Hits[:,:, 0], Hits[:,:, 1], marker='o',color="xkcd:dark grey", s=0.4, linewidths=0)
    #tpc,tri,emi= tomo(f"{patch_mat_file[i]}",list_ax[3*i],emi_vmin=0, emi_vmax=emi_vals_max[i])
    tpc,tri,emi= tomo(f"{patch_mat_file[i]}",list_ax[3*i+1],emi_vmin=0, emi_vmax=emi_vals_max[i])
    tpc,tri,emi= tomo(f"{patch_mat_file[i]}",list_ax[3*i+2],emi_vmin=0, emi_vmax=emi_vals_max[i])
    
    for mf_name in name_mf:
        if i != 3 or (i==3 and mf_name in allowed_mf_3):
         manifold = Manifold.load(f"{repository_path[i]}manifolds_P/OS_BO78/{mf_name}")
         manifold.plot(stepsize_limit=0.05, ax=list_ax[3*i], markersize=0, lw=1.3 if (mf_name=='mf_1B.pkl') else 0.5,
                       colors=["rosybrown","black"] if (mf_name=='mf_1B.pkl') else ["rosybrown", "xkcd:red"], 
                       labels=['stable MF','unstable MF'] if (mf_name=='mf_1T.pkl' and i==3) else [None ,None])
         manifold.plot(stepsize_limit=0.05, ax=list_ax[3*i+1], markersize=0, lw=0.9 if (mf_name=='mf_1B.pkl') else 0.5,
                       colors=["rosybrown","white"] if (mf_name=='mf_1B.pkl') else ["rosybrown", "xkcd:red"])
         
    for j in range(len(type_fp)):
        if i != 3 or (i==3 and allowed_fp_3[j]):
         fp = FixedPoint(section)
         fp.find(1, first_guess[j], method='scipy.root')
         fp.plot(ax=list_ax[3*i], marker='o' if 'o' in type_fp[j] else 'x', color="xkcd:dark blue",zorder=10,label=label_fp[j] if i==3 else None)
         fp.plot(ax=list_ax[3*i+1], marker='o' if 'o' in type_fp[j] else 'x', color="xkcd:white",zorder=10) 
         
         if j==1:
          list_trace.append(section.df(1,fp.coords[0]).trace())
          fp.plot(ax=list_ax[3*i+2], marker='o' if 'o' in type_fp[j] else 'x', color="xkcd:white",zorder=10)

    list_ax[3*i].set_xlim(0.62, 1.15)
    list_ax[3*i].set_ylim(-0.75, 0.75)

    list_ax[3*i].set_xlabel(r"$R[m]$")
    list_ax[3*i].set_aspect('equal') 
    list_ax[3*i].set_title(f"{titles[i]}")
    list_ax[3*i+1].set_title(f"{title_MANTIS[i]}")

    list_ax[3*i+1].set_xlim(0.7, 1.0)
    list_ax[3*i+1].set_ylim(-0.6, -0.1)
    #list_ax[3*i+1].set_xlabel(r"$R[m]$")
    list_ax[3*i+1].set_aspect('equal')

    list_ax[3*i+2].set_xlim(0.7, 1.0)
    list_ax[3*i+2].set_ylim(-0.6, -0.1)
    list_ax[3*i+2].set_xlabel(r"$R[m]$")
    list_ax[3*i+2].set_aspect('equal')


list_ax[0].set_ylabel(r"$Z[m]$")
list_ax[1].set_ylabel(r"$Z[m]$")
list_ax[2].set_ylabel(r"$Z[m]$")
list_ax[3].set_yticklabels([])
list_ax[4].set_yticklabels([])
list_ax[5].set_yticklabels([])
list_ax[6].set_yticklabels([])
list_ax[7].set_yticklabels([])
list_ax[8].set_yticklabels([])
list_ax[9].set_yticklabels([])
list_ax[10].set_yticklabels([])
list_ax[11].set_yticklabels([])






##### STANDARD SHOT 80064 BO78 OS ######

STfield = AxisymmetricCylindricalGridField.from_matlab_file('./script/standard/ST_80064_11_BO78_OS.mat', with_perturbation=True)

section = CylindricalBfieldSection(STfield,phi0=1.9,R0=0.88, Z0=0)

top_o = FixedPoint(section)
top_o.find(1, [0.89,0.3], method='scipy.root')

x_point1= FixedPoint(section)
x_point1.find(1, [0.7,-0.16], method='scipy.root')

manifold_1T  = Manifold.load('./script/standard/manifolds_P/OS_BO78/mf_1T.pkl')
manifold_1B  = Manifold.load('./script/standard/manifolds_P/OS_BO78/mf_1B.pkl')

colors=["xkcd:dark blue", "xkcd:white"]

for i in range(2):
    top_o.plot(ax=list_ax[12+i], marker='o', color=colors[i],zorder=10,label='O-point' if i==0 else None)
    x_point1.plot(ax=list_ax[12+i], marker='x', color=colors[i],zorder=10,label='X-point' if i==0 else None)
    manifold_1T.plot(stepsize_limit=0.05, ax=list_ax[12+i], markersize=0, lw=0.5,
                     colors=["rosybrown", "xkcd:red"], labels=['stable MF','unstable MF'] if i==0 else [None ,None])
    manifold_1B.plot(stepsize_limit=0.05, ax=list_ax[12+i], markersize=0, lw=0.9,
                     colors=["rosybrown", "xkcd:black"] if i==0 else ["rosybrown", "xkcd:white"],
                     labels=[None ,'divertor MF'] if i==0 else [None ,None])


x_point1.plot(ax=list_ax[14], marker='x', color=colors[1],zorder=10,label='X-point')

Hits=np.load('./script/standard/ST_80064_pp_hits.npy')

list_ax[12].scatter(Hits[:,:, 0], Hits[:,:, 1],marker='o', color="xkcd:dark grey", s=0.4, linewidths=0)

tpc,tri,emi= tomo('./script/standard/ST_80064_patch_C3.mat',ax=list_ax[13],emi_vmin=0, emi_vmax=8e20)
tpc,tri,emi= tomo('./script/standard/ST_80064_patch_C3.mat',ax=list_ax[14],emi_vmin=0, emi_vmax=8e20)


list_ax[12].set_xlim(0.62, 1.15)
list_ax[12].set_ylim(-0.75, 0.75)

list_ax[12].set_xlabel(r"$R[m]$")
list_ax[12].set_aspect('equal') 
list_ax[12].set_title(f"{titles[4]}")
list_ax[12].set_yticklabels([])

list_ax[13].set_xlim(0.7, 1.0)
list_ax[13].set_ylim(-0.6, -0.1)
list_ax[13].set_aspect('equal')
list_ax[13].set_yticklabels([])
list_ax[13].set_title(f"{title_MANTIS[4]}")

list_ax[14].set_xlim(0.7, 1.0)
list_ax[14].set_ylim(-0.6, -0.1)
list_ax[14].set_xlabel(r"$R[m]$")
list_ax[14].set_aspect('equal')
list_ax[14].set_yticklabels([])


for ax, lab in zip(list_ax, subscripts):
    if ax in (ax1, ax2, ax3, ax4, ax5):
        text_color = 'black'
    else:
        text_color = 'white'
    ax.text(
        0.95, 0.95, lab,
        transform=ax.transAxes,
        fontsize=10,
        fontweight='bold',
        ha='right',
        va='top',
        color=text_color,
        #bbox=dict(facecolor=bbox_face, edgecolor='none', pad=2, alpha=0.8),
        zorder=200
    )

cax = inset_axes(list_ax[13], width="4%", height="80%", loc='center right',
                 bbox_to_anchor=(0.2, -0.1, 1, 1), bbox_transform=list_ax[13].transAxes, borderpad=0)



list_ax[12].legend(loc='lower center', bbox_to_anchor=(0.41, 0.01), ncol=2, fontsize=7)#framealpha=0.6)

# create colorbar and hide the scientific offset text (e.g. "1e20")
cbar = fig.colorbar(tpc, cax=cax, label='Relative emission')
# hide the offset text that matplotlib places (the "1e20" above/beside the bar)
cbar.ax.yaxis.get_offset_text().set_visible(False)
# refresh ticks/formatter
cbar.update_ticks()
# ...existing code...

#plt.savefig("./figures/MANTIS_MF_v4.png", dpi=900, bbox_inches="tight", pad_inches=0, facecolor=fig.get_facecolor())
plt.show()
print(f"Trace of separatrix Jacobians: {list_trace[0]}, {list_trace[1]}, {list_trace[2]}, {list_trace[3]}")



