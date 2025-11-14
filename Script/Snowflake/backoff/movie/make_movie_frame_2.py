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

    SFField = AxisymmetricCylindricalGridField.from_matlab_file('./Script/Snowflake/backoff/Snowflake_09_Backoff.mat', with_perturbation=True)

    section = CylindricalBfieldSection(SFField,phi0=angle, R0=0.88, Z0=0)


    top_o = FixedPoint(section)
    top_o.find(1, [0.9,0.0], method='scipy.root')
    top_o_coord = top_o.coords[0]
    top_o.plot(ax=ax, marker='o', color="xkcd:crimson")


    x_point= FixedPoint(section)
    x_point.find(1, [0.75,-0.4], method='scipy.root')
    x_point_coord = x_point.coords[0]
    x_point.plot(ax=ax, marker='x', color="xkcd:crimson")

    x_point2= FixedPoint(section)
    x_point2.find(1, [0.8,-0.6], method='scipy.root')
    x_point_coord2 = x_point2.coords[0]
    x_point2.plot(ax=ax, marker='x', color="xkcd:crimson")

    # pplot = PoincarePlot.with_segments(section, [top_o_coord, x_point_coord, x_point_coord, x_point_coord2],[20, 5],connected=False)
    # pplot.compute(1)
    # pplot.plot(ax=ax, color="xkcd:dark grey")

    # #top fp top manifold
    manifold1 = Manifold(section, x_point, x_point)
    manifold1.compute(
      eps_s=9e-6, eps_u=8e-6, nint_s=10, nint_u=8, neps_s=160, neps_u=160)
    manifold1.plot(stepsize_limit=0.2,ax=ax, markersize=0, lw=0.7)


    textstr = rf"perturbed manifold for angle = {angle:.3f}$"
    props = dict(boxstyle='round', facecolor='white', edgecolor='black')
    ax.text(0.99, 0.99, textstr, transform=ax.transAxes, fontsize=6,
            verticalalignment='top', horizontalalignment='right', bbox=props)

    ax.set_xlim(0.62, 1.15)
    ax.set_ylim(-0.75, 0.75)
    ax.set_xlabel(r"$R[m]$")
    ax.set_ylabel(r"$Z[m]$")
    ax.set_title("Snowflake Backoff - Perturbed Manifolds")

    plt.savefig(f"./Script/Snowflake/backoff/movie/movie_frames_annotated/pertmovie_{frame_index:03d}.png", bbox_inches='tight', dpi=720)
    
    return

    


if __name__=="__main__":
    if len(sys.argv) != 3:
        raise ValueError("Please provide the angle as argument")
    angle = float(sys.argv[1])
    figure_index = int(sys.argv[2])
    plot_perturbed_manifolds(angle, frame_index=figure_index)
    print(f"Frame {figure_index} done")
    sys.exit(0)

    