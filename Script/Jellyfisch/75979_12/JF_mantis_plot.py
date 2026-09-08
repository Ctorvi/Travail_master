import os, sys
proj_root = os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..','..'))  # remonte jusqu'à "Travail master"
if proj_root not in sys.path:
    sys.path.insert(0, proj_root)
from script.function.field_utils import plot_tomographic as tomo
from script.function.field_utils import plot_LIUQE as LIUQE
from pyoculus.fields import AxisymmetricCylindricalGridField
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
import numpy as np
from matplotlib import pyplot as plt
import logging
from scipy.io import loadmat
from matplotlib.lines import Line2D
from matplotlib.tri import Triangulation
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import MaxNLocator





# fig = plt.figure(figsize=(6, 8))
# # gs = GridSpec(2, 3, hspace=0.25, width_ratios=[1, 1, 0.9], height_ratios=[1, 1])
# gs = GridSpec(2, 2, hspace=0.25, width_ratios=[ 1, 0.9], height_ratios=[1, 1])

# ax_l = fig.add_subplot(gs[:, 0])    
# ax_tr = fig.add_subplot(gs[0, 1]) 
# ax_br = fig.add_subplot(gs[1, 1])


#### oral p
fig1,ax_l = plt.subplots(1,1, figsize=(6,6.6))
fig2,ax_tr = plt.subplots(1,1, figsize=(5,7))
fig3,ax_br = plt.subplots(1,1, figsize=(5,7))


pos = ax_tr.get_position()
dy = 0.07 
ax_tr.set_position([pos.x0, pos.y0 - dy, pos.width, pos.height])


# fig.tight_layout()

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

pert_mat_file = 'mat_files/JF_75979_120.mat'

patch_mat_file = './script/jellyfisch/75979_12/mat_files/JF_75979_patch_120_He2.mat'

mf_list = ['mf_1T', 'mf_1B', 'mf_2T', 'mf_2B', 'mf_3L', 'mf_3R', 'mf_4T', 'mf_4B']

JFField = AxisymmetricCylindricalGridField.from_matlab_file(f"{repository_path}{pert_mat_file}", with_perturbation=True)

section = CylindricalBfieldSection(JFField,phi0=1.9,R0=0.88, Z0=0)

top_o = FixedPoint(section)
top_o.find(1, [0.9,0.18], method='scipy.root')

x_point1= FixedPoint(section)
x_point1.find(1, [0.8,-0.27], method='scipy.root')

x_point2= FixedPoint(section)
x_point2.find(1, [0.8,-0.57], method='scipy.root')

x_point3= FixedPoint(section)
x_point3.find(1, [1.05,-0.58], method='scipy.root')

x_point4= FixedPoint(section)
x_point4.find(1, [0.75,0.65], method='scipy.root')

###########################  tomographic reconstruction  #########################

tpc1,tri,emi= tomo(f"{patch_mat_file}",ax_l,emi_vmin=0, emi_vmax=2.25e19)
tpc2,tri,emi= tomo(f"{patch_mat_file}",ax_tr,emi_vmin=0, emi_vmax=2.25e19)
tpc3,tri,emi= tomo(f"{patch_mat_file}",ax_br,emi_vmin=0, emi_vmax=2.25e19)
#tpc4,tri,emi= tomo(f"{patch_mat_file}",ax_ll,emi_vmin=0, emi_vmax=2.25e19)

##### loading and plotting  ##

manifs = {name: {name} for name in mf_list}

for i in mf_list:
       
    if i == 'mf_1T':
        lbls = ['stable MF', 'unstable MF']
    elif i == 'mf_1B':
        lbls = [None, 'divertor MF']   # troisième condition
    else:
        lbls = [None, None]

    manifs[i] = Manifold.load(f"{repository_path}{file_manifolds}{i}.pkl")
   # manifs[i].plot(ax=ax2, markersize=0, lw=0.7,labels=[None,None] if i!=0 else ['stable MF','unstable MF'])
    manifs[i].plot(stepsize_limit=0.05, ax=ax_l, markersize=0, lw=1.3 if (i=='mf_1B') else 0.7,
                       colors=["rosybrown","xkcd:white"] if (i=='mf_1B') else ["rosybrown", "xkcd:red"], 
                       labels=lbls)
    manifs[i].plot(stepsize_limit=0.05, ax=ax_tr, markersize=0, lw=1.3 if (i=='mf_1B') else 0.7,
                       colors=["rosybrown","xkcd:white"] if (i=='mf_1B') else ["rosybrown", "xkcd:red"], 
                       labels=lbls)




top_o.plot(ax=ax_l, marker='o', s=40, color="xkcd:white",label='O-point', zorder=10)
x_point1.plot(ax=ax_l, marker='x', s=50, color="xkcd:white",label='X-points', zorder=10)
x_point2.plot(ax=ax_l, marker='x', s=50, color="xkcd:white",label=None, zorder=10)
x_point3.plot(ax=ax_l, marker='x', s=50, color="xkcd:white",label=None, zorder=10)
x_point4.plot(ax=ax_l, marker='x', s=50, color="xkcd:white",label=None, zorder=10)

top_o.plot(ax=ax_tr, marker='o', s=40, color="xkcd:white",label='O-point', zorder=10)
x_point1.plot(ax=ax_tr, marker='x', s=50, color="xkcd:white",label='X-points', zorder=10)
x_point2.plot(ax=ax_tr, marker='x', s=50, color="xkcd:white",label=None, zorder=10)
x_point3.plot(ax=ax_tr, marker='x', s=50, color="xkcd:white",label=None, zorder=10)
x_point4.plot(ax=ax_tr, marker='x', s=50, color="xkcd:white",label=None, zorder=10)

# list_ax = [ax_l, ax_tr, ax_br]
# subscripts = ['(a)', '(b)', '(c)']
# text_color = 'white'


# for ax, lab in zip(list_ax, subscripts):
#     ax.text(
#         0.05, 0.9, lab,
#         transform=ax.transAxes,
#         fontsize=12,
#         fontweight='bold',
#         ha='left',
#         va='top',
#         color=text_color,
#         #bbox=dict(facecolor=bbox_face, edgecolor='none', pad=2, alpha=0.8),
#         zorder=200
#     )


# ax_ll.set_xlim(0.62, 1.15)
# ax_ll.set_ylim(-0.75, 0.75)
# ax_ll.set_xlabel(r"$R[m]$")
# ax_ll.set_aspect('equal') 
# #ax_l.set_title('Poincaré with perturbation')
# ax_ll.legend(loc='lower right', bbox_to_anchor=(1.0, 0.87), ncol=1, fontsize=9)



ax_l.set_ylabel(r"$Z[m]$")


ax_l.set_xlim(0.62, 1.15)
ax_l.set_ylim(-0.75, 0.75)
ax_l.set_xlabel(r"$R[m]$")
ax_l.set_aspect('equal') 
#ax_l.set_title('Poincaré with perturbation')
ax_l.legend(loc='lower right', bbox_to_anchor=(1.0, 0.87), ncol=1, fontsize=9)

ax_tr.set_xlim(0.65, 0.95)
ax_tr.set_ylim(-0.6, -0.25)
ax_tr.set_aspect('equal') 
#ax_tr.set_xticklabels([])

ax_tr.set_xlabel(r"$R[m]$")
ax_tr.set_ylabel(r"$Z[m]$")


ax_br.set_xlim(0.65, 0.95)
ax_br.set_ylim(-0.6, -0.25)
ax_br.set_aspect('equal') 
ax_br.set_xlabel(r"$R[m]$")

ax_br.set_ylabel(r"$Z[m]$")


cax = inset_axes(ax_tr, width="100%", height="7%", loc='lower center',
                  bbox_to_anchor=(0., 1.1, 1, 1), bbox_transform=ax_tr.transAxes, borderpad=0)
cbar = fig2.colorbar(tpc3, cax=cax, orientation='horizontal')
cbar.ax.xaxis.set_label_position('top')
cbar.ax.xaxis.tick_top()
cbar.set_label('Relative emission (He II line)', labelpad=6)  # ajuster labelpad si besoin
cbar.ax.xaxis.get_offset_text().set_visible(False)
cbar.update_ticks()



cax = inset_axes(ax_br, width="100%", height="7%", loc='lower center',
                  bbox_to_anchor=(0., 1.1, 1, 1), bbox_transform=ax_br.transAxes, borderpad=0)
cbar = fig3.colorbar(tpc3, cax=cax, orientation='horizontal')
cbar.ax.xaxis.set_label_position('top')
cbar.ax.xaxis.tick_top()
cbar.set_label('Relative emission (He II line)', labelpad=6)  # ajuster labelpad si besoin
cbar.ax.xaxis.get_offset_text().set_visible(False)
cbar.update_ticks()




# for ax in (ax_tr, ax_br):
#     # y axis on the right
#     ax.yaxis.set_label_position('right')
#     ax.yaxis.tick_right()
#     ax.yaxis.set_ticks_position('right')

#     # show right spine, hide left spine for a cleaner look
#     ax.spines['right'].set_visible(True)
#     ax.spines['left'].set_visible(False)
#     ax.yaxis.set_major_locator(MaxNLocator(nbins=4, prune='both'))



#plt.savefig(f"{repository_path}figures/JF_plot_MANTIS_chapter.png", dpi=150, bbox_inches="tight", pad_inches=0, facecolor=fig.get_facecolor())
fig2.savefig(f"{repository_path}/oral_MF_only_divertor.png", bbox_inches='tight', pad_inches=0.0, dpi=720)
#fig3.savefig(f"{repository_path}/oral_no_MF.png", bbox_inches='tight', pad_inches=0.0, dpi=720)



plt.show()




