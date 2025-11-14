from pyoculus.fields import AxisymmetricCylindricalGridField
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
import numpy as np
from matplotlib import pyplot as plt
import logging

fig, ax = plt.subplots(1, 1, figsize=(6, 6.6))

plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.size": 10,
    }
)
logging.basicConfig(level=logging.DEBUG)


JFField = AxisymmetricCylindricalGridField.from_matlab_file('./Script/Jellyfisch/77021_1125/Jellyfisch1125.mat', with_perturbation=False)


section = CylindricalBfieldSection(JFField,R0=0.88, Z0=0)


top_o = FixedPoint(section)
top_o.find(1, [0.9,0.18], method='scipy.root')
top_o_coord = top_o.coords[0]
top_o.plot(ax=ax, marker='o', color="xkcd:crimson")


x_point= FixedPoint(section)
x_point.find(1, [0.8,-0.2], method='scipy.root')
x_point_coord = x_point.coords[0]
x_point.plot(ax=ax, marker='x', color="xkcd:crimson")


# manifold2 = Manifold(section, x_point, x_point,-x_point_coord+top_o_coord, -x_point_coord+top_o_coord)
# manifold2.compute(
#    eps_s=9e-6, eps_u=8e-6, nint_s=10, nint_u=12, neps_s=80, neps_u=240)
# manifold2.plot(ax=ax, markersize=0, lw=1.5)


# manifold1 = Manifold(section, x_point, x_point)
# manifold1.compute(
#     eps_s=9e-6, eps_u=8e-6, nint_s=10, nint_u=10, neps_s=80, neps_u=80)
# manifold1.plot(ax=ax, markersize=0, lw=1.5)


###poincare perturbed


#pplot = PoincarePlot.with_segments(section, [top_o_coord, x_point_coord,x_point_coord,[0.9,-0.15]],[50, 30],connected=False)
pplot=PoincarePlot.with_linspace(section, top_o_coord, x_point_coord, 40)# connected=False)
pplot.compute(400)
#pplot.plot(ax=ax, color="xkcd:dark grey")




    

np.save('./Script/Jellyfisch/77021_1125/Jellyfisch_Pert.npy', pplot._hits)
#Hits=np.load('./Script/Jellyfisch/77021_1125/Jellyfisch_Pert.npy')

# ratio=2.8158

ax.scatter(Hits[:,:, 0], Hits[:,:, 1], color="xkcd:dark grey", s=1.4, linewidths=0)
ax.set_xlim(0.62, 1.15)
ax.set_ylim(-0.75, 0.75)
ax.set_xlabel(r"$R[m]$")
ax.set_ylabel(r"$Z[m]$")
ax.set_aspect('equal') 



plt.savefig('./Script/Jellyfisch/77021_1125/Jellyfisch_Pert_2.png', bbox_inches='tight', dpi=720)


plt.show()
