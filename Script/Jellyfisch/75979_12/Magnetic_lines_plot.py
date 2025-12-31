import os, sys
proj_root = os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..','..'))  # remonte jusqu'à "Travail master"
if proj_root not in sys.path:
    sys.path.insert(0, proj_root)
from script.function.field_utils import plot_tomographic as tomo
from script.function.field_utils import plot_LIUQE as LIUQE
from matplotlib.ticker import FormatStrFormatter, MaxNLocator, FuncFormatter

from mpl_toolkits.axes_grid1 import make_axes_locatable
from pyoculus.fields import AxisymmetricCylindricalGridField
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
import numpy as np
from matplotlib import pyplot as plt
import logging
from scipy.io import loadmat
from matplotlib.lines import Line2D
from matplotlib.tri import Triangulation
from shapely.geometry import Polygon
from matplotlib.patches import Polygon as MplPolygon


def fill_between_polygons(ax, R_outer, Z_outer, R_inner, Z_inner, **kwargs):
    R1 = np.asarray(R_outer).ravel()
    Z1 = np.asarray(Z_outer).ravel()
    R2 = np.asarray(R_inner).ravel()
    Z2 = np.asarray(Z_inner).ravel()

    # retirer NaN éventuels
    mask1 = ~np.isnan(R1) & ~np.isnan(Z1)
    mask2 = ~np.isnan(R2) & ~np.isnan(Z2)
    R = np.concatenate([R1[mask1], R2[mask2][::-1]])
    Z = np.concatenate([Z1[mask1], Z2[mask2][::-1]])

    ax.fill(R, Z, picker=False, label='_nolegend_', **kwargs)


def fill_between_polygons_shapely(ax, R_outer, Z_outer, R_inner, Z_inner, color=(0.5,0.5,0.5)):

    outer = Polygon(np.column_stack([R_outer, Z_outer]))
    inner = Polygon(np.column_stack([R_inner, Z_inner]))
    area = outer.difference(inner)  # Polygon or MultiPolygon

    if area.is_empty:
        return

    # background pour "creuser" les trous si besoin
    bg = ax.get_facecolor()
    polys = area.geoms if hasattr(area, "geoms") else [area]

    for poly in polys:
        ext = np.array(poly.exterior.coords)
        ax.add_patch(MplPolygon(ext, closed=True, facecolor=color, edgecolor="none",picker=False, label='_nolegend_'))
        # dessiner les intérieurs (trous) avec la couleur de fond pour creuser
        for interior in poly.interiors:
            hole = np.array(interior.coords)
            ax.add_patch(MplPolygon(hole, closed=True, facecolor=bg, edgecolor="none",picker=False, label='_nolegend_'))




fig,(ax1, ax2) = plt.subplots(1, 2, figsize=(7, 6.6))


plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.size": 9,
    }
)


file_manifolds = 'manifolds_P/'


logging.basicConfig(level=logging.DEBUG)

repository_path = './script/jellyfisch/75979_12/'

pert_mat_file = '/mat_files/JF_75979_120_BO78_OS.mat'

patch_mat_file = './script/jellyfisch/75979_12/mat_files/JF_75979_patch_120_He2.mat'

JFField = AxisymmetricCylindricalGridField.from_matlab_file(f"{repository_path}{pert_mat_file}", with_perturbation=True)

section = CylindricalBfieldSection(JFField,phi0=0,R0=0.88, Z0=0)

LIUQE_mat_file = '/mat_files/Liuqe_JF_75979_120.mat'

data_LIUQE, lgd1,lgd2= LIUQE(f"{repository_path}{LIUQE_mat_file}",n_levels=15,ax=ax2,colors=["grey", "grey", "grey"],linestyles=['-','-','-'])

data_LIUQE, lgd1,lgd2= LIUQE(f"{repository_path}{LIUQE_mat_file}",n_levels=15,ax=ax1,colors=["grey", "grey", "grey"],linestyles=['-','-','-'])


data_vessel = loadmat(f"{repository_path}/mat_files/Liuqe_JF_75979_120_with_V.mat", squeeze_me=True, struct_as_record=False)

vessel=data_vessel['vessel']

R_in=vessel.R_in  
Z_in=vessel.Z_in
R_out=vessel.R_out
Z_out=vessel.Z_out
Rt=vessel.Rt
Zt=vessel.Zt



R_in_anti=R_in[::-1]
Z_in_anti=Z_in[::-1]



R_max =JFField.R.max()
R_min = JFField.R.min()
Z_max = JFField.Z.max()
Z_min = JFField.Z.min()


R_vec = np.linspace(R_min, R_max, num=56)  
Z_vec = np.linspace(Z_min, Z_max, num=130)  


BR=np.zeros((Z_vec.shape[0],R_vec.shape[0]))
BZ=np.zeros((Z_vec.shape[0],R_vec.shape[0]))


R_grid = np.tile(R_vec, (Z_vec.size, 1))          

Z_grid = np.tile(Z_vec[:, None], (1, R_vec.size)) 

R_grid2, Z_grid2 = np.meshgrid(R_vec, Z_vec, indexing='xy')
phi = 1.9

for i in range(R_vec.shape[0]):
    for j in range(Z_vec.shape[0]):
        BR[j,i] = JFField.pertfield.B_R([R_grid2[j,i], phi, Z_grid2[j,i]])
        BZ[j,i] = JFField.pertfield.B_Z([R_grid2[j,i], phi, Z_grid2[j,i]])

BR=1000*BR
BZ=1000*BZ

n_levels = 8

for ax in (ax1, ax2):
  #fill_between_polygons_shapely(ax, R_in, Z_in, R_out, Z_out)
  fill_between_polygons_shapely(ax, R_in_anti, Z_in_anti,Rt,Zt )



# s'assurer qu'il n'y a pas de NaN (ou les remplacer par min)


bz_min, bz_max = np.nanmin(BZ), np.nanmax(BZ)

if np.isfinite(bz_min) and np.isfinite(bz_max) and bz_max > bz_min:
    levels = np.linspace(bz_min, bz_max, n_levels)
else:
    levels = np.linspace(-1.0, 1.0, n_levels)

cs = ax2.contour(R_grid2, Z_grid2, BZ, levels=levels, colors='k', linewidths=1.5,linestyles='--')
#ax2.clabel(cs, fmt=lambda v: f"{v*1e3:.1f}", inline=True, fontsize=12)
ax2.clabel(cs, fmt='%.1f', inline=True, fontsize=12)

cf = ax2.contourf(R_grid2, Z_grid2, BZ, levels=levels, cmap='viridis', alpha=0.6)
# cbar = plt.colorbar(cf, ax=ax2, pad=0.02, orientation='vertical')
# cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
# cbar.ax.yaxis.get_offset_text().set_visible(False)
# cbar.update_ticks()
# cbar.set_label('B_Z [mT]', labelpad=6)


divider = make_axes_locatable(ax2)
cax = divider.append_axes("right", size="4%", pad=0.06)
nticks = n_levels
ticks = np.linspace(levels[0], levels[-1], nticks)
cbar = fig.colorbar(cf, cax=cax, orientation='vertical', ticks=ticks)
cbar.ax.set_yticklabels([f"{t:.1f}" for t in ticks])   # affiche avec 1 décimale
cbar.ax.yaxis.get_offset_text().set_visible(False)


# axis formatting
ax2.set_xlim(R_min, R_max)
ax2.set_ylim(Z_min, Z_max)
ax2.set_xlabel(r"$R\ [m]$")
ax2.set_aspect('equal')
ax2.set_title(r"$B_Z$ [mT]")



br_min, br_max = np.nanmin(BR), np.nanmax(BR)
if np.isfinite(br_min) and np.isfinite(br_max) and br_max > br_min:
    levels = np.linspace(br_min, br_max, n_levels)
else:
    levels = np.linspace(-1.0, 1.0, n_levels)


cs = ax1.contour(R_grid2, Z_grid2, BR, levels=levels, colors='k', linewidths=1.5,linestyles='--')
ax1.clabel(cs,fmt='%.1f', inline=True, fontsize=12)

cf = ax1.contourf(R_grid2, Z_grid2, BR, levels=levels, cmap='viridis', alpha=0.6)

divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="4%", pad=0.06)
nticks = n_levels
ticks = np.linspace(levels[0], levels[-1], nticks)
cbar = fig.colorbar(cf, cax=cax, orientation='vertical', ticks=ticks)
cbar.ax.set_yticklabels([f"{t:.1f}" for t in ticks])   # affiche avec 1 décimale
cbar.ax.yaxis.get_offset_text().set_visible(False)



# cbar.ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x*1e3:.1f}"))
# ticks = cbar.get_ticks()
# ticks_round = np.round(ticks, 1)
# cbar.set_ticks(ticks_round)
# #cbar.set_ticklabels([f"{t*1e3:.1f}" for t in ticks_round])
# cbar.ax.yaxis.get_offset_text().set_visible(False)
# cbar.ax.yaxis.set_major_locator(MaxNLocator(nbins=8))
# cbar.set_label('B_R [mT]', labelpad=6)

# axis formatting
ax1.set_xlim(R_min, R_max)
ax1.set_ylim(Z_min, Z_max)
ax1.set_xlabel(r"$R\ [m]$")
ax1.set_ylabel(r"$Z\ [m]$")
ax1.set_aspect('equal')
ax1.set_title(r"$B_R$ [mT]")


R_in = np.asarray(R_in).ravel()
Z_in = np.asarray(Z_in).ravel()
R_out = np.asarray(R_out).ravel()
Z_out = np.asarray(Z_out).ravel()
Rt = np.asarray(Rt).ravel()
Zt = np.asarray(Zt).ravel()


# # # first fill: fill([R_in' R_out'],[Z_in' Z_out'],'k',...)
# coords_x = np.concatenate([R_in, R_out[::-1]])
# coords_y = np.concatenate([Z_in, Z_out[::-1]])
# for ax in (ax1, ax2):
#     patch = ax.fill(coords_x, coords_y, color='k', label='_nolegend_')[0]
#     patch.set_picker(False)   # équivalent hittest='off'
#     patch.set_zorder(0)       # optionnel : placer derrière

# # second fill: fill([R_in' Rt'],[Z_in' Zt'],[0.5 0.5 0.5],...)
# coords2_x = np.concatenate([R_in, Rt[::-1]])
# coords2_y = np.concatenate([Z_in, Zt[::-1]])
# for ax in (ax1, ax2):
#     patch2 = ax.fill(coords2_x, coords2_y, color=(0.5, 0.5, 0.5), label='_nolegend_')[0]
#     patch2.set_picker(False)
#     patch2.set_zorder(1)

# plot the tip curve Rt,Zt as black line, no legend entry
for ax in (ax1, ax2):
    line = ax.plot(Rt, Zt, color='k', label='_nolegend_')[0]
    line.set_picker(False)
    line.set_zorder(2)




#plt.savefig("./figures/Field_lines_n1_field_v1_0.png", dpi=150, bbox_inches=None, facecolor=fig.get_facecolor())
plt.show()
print(Rt,Zt)