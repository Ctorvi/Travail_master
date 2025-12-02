from pyoculus.fields import AnalyticCylindricalBfield
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
from pyoculus.fields import AxisymmetricCylindricalGridField
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from jax import jit, config
import numpy as np
from pathlib import Path
import sys
import logging
from scipy.io import loadmat
from matplotlib.tri import Triangulation
config.update("jax_enable_x64", True)




plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.size": 10,
    }
)
logging.basicConfig(level=logging.DEBUG)

fig, ax = plt.subplots(1,1, figsize=(4, 6))
ax.set_aspect("equal")

def plot_perturbed_manifolds(angle, frame_index=0):

    print(f"Plotting perturbed manifolds for angle = {angle}")

    JFField = AxisymmetricCylindricalGridField.from_matlab_file('./Script/Jellyfisch/77021_105/Jellyfisch_105.mat', with_perturbation=True)

    section = CylindricalBfieldSection(JFField,phi0=angle, R0=0.88, Z0=0)


    top_o = FixedPoint(section)
    top_o.find(1, [0.9,0.18], method='scipy.root')
    top_o_coord = top_o.coords[0]
    top_o.plot(ax=ax, marker='o', color="xkcd:crimson")


    x_point= FixedPoint(section)
    x_point.find(1, [0.8,-0.2], method='scipy.root')
    x_point_coord = x_point.coords[0]
    x_point.plot(ax=ax, marker='x', color="xkcd:crimson")



    manifold2 = Manifold(section, x_point, x_point,-x_point_coord+top_o_coord, -x_point_coord+top_o_coord)
    manifold2.compute(
    eps_s=9e-6, eps_u=8e-6, nint_s=8, nint_u=14, neps_s=80, neps_u=240)#240 8 14
    manifold2.plot(stepsize_limit=0.2,ax=ax, markersize=0, lw=0.7)

    # pplot = PoincarePlot.with_segments(section, [top_o_coord, x_point_coord, x_point_coord, x_point_coord2],[20, 5],connected=False)
    # pplot.compute(1)
    # pplot.pl.  ot(ax=ax, color="xkcd:dark grey")


    data = loadmat('./Script/Jellyfisch/77021_105/JF_patch.mat', squeeze_me=True, struct_as_record=False)

    inv_grid=data['inv_grid']

    tri_x = inv_grid.tri_x
    tri_y = inv_grid.tri_y
    tri_nodes = inv_grid.tri_nodes
    R_values = inv_grid.R_values
    Z_values = inv_grid.Z_values
    R_mean = inv_grid.R_mean
    Z_mean = inv_grid.Z_mean
    Area= inv_grid.Area


    emi=data['emi']
    t=data['t']
    time=data['time']


    time_vec = np.asarray(time)
# Find index of closest time
    t_idx = np.argmin(np.abs(time_vec - t))
    emi = emi[:, t_idx]
    
  
    triangles = tri_nodes.astype(int)  # tri_nodes doit être zéro-indexé

    tri = Triangulation(tri_x, tri_y, triangles)


    
    # Plot the patch
#plt.figure()
    tpc=ax.tripcolor(tri, emi, shading='flat', edgecolors='none', vmin=0, vmax=4.5e20)
    fig.colorbar(tpc, ax=ax, label='Émission')
   

    textstr = rf"perturbed manifold for angle = {angle:.3f}$"
    props = dict(boxstyle='round', facecolor='white', edgecolor='black')
    ax.text(0.99, 0.99, textstr, transform=ax.transAxes, fontsize=6,
            verticalalignment='top', horizontalalignment='right', bbox=props)

    ax.set_xlim([0.62, 1.13])
    ax.set_ylim([-0.75, 0.75]) 
    ax.set_xlabel(r"$R[m]$")
    ax.set_ylabel(r"$Z[m]$")

    plt.savefig(f"./Script/Jellyfisch/77021_105/movie/movie_frames_annotated/pertmovie_{frame_index:03d}.png", bbox_inches='tight', dpi=720)
    
    return

    


if __name__=="__main__":
    if len(sys.argv) != 3:
        raise ValueError("Please provide the angle as argument")
    angle = float(sys.argv[1])
    figure_index = int(sys.argv[2])
    plot_perturbed_manifolds(angle, frame_index=figure_index)
    print(f"Frame {figure_index} done")
    sys.exit(0)

    