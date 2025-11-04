from pyoculus.fields import AxisymmetricCylindricalGridField
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
import numpy as np
from matplotlib import pyplot as plt
import logging


plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.size": 10,
    }
)
logging.basicConfig(level=logging.DEBUG)

DoubletField = AxisymmetricCylindricalGridField.from_matlab_file('./Script/Doublet/struct_for_chris.mat', with_perturbation=True)


section = CylindricalBfieldSection(DoubletField, R0=0.88, Z0=0,domain=[(-np.inf, np.inf), (-np.inf, np.inf)], rtol=1e-10)

fig, ax = plt.subplots(1, 1, figsize=(6, 6.6))

Hits=np.load('Hits/Doublet_Island.npy')

#guess_t = [0.8852,0.42]
#top_o = FixedPoint(section)
#top_o.find(1, guess_t, method='scipy.root')
#top_o_coord = top_o.coords[0]
#top_o.plot(ax=ax, color="xkcd:crimson")


#guess_b = [0.885, -0.42]
#bottom_o = FixedPoint(section)
#bottom_o.find(1, guess_b,  method='scipy.root')
#bottom_o_coord = bottom_o.coords[0]
#bottom_o.plot(ax=ax, color="xkcd:crimson")


#xpoint= FixedPoint(section)
#xpoint.find(1, [0.8865, -0.0001], method='scipy.root')
#x_point_coord = xpoint.coords[0]
#xpoint.plot(ax=ax, color="xkcd:crimson")


#pplot = PoincarePlot.with_segments(section, [top_o_coord, x_point_coord, bottom_o_coord, x_point_coord], [10, 10], connected=False)
#pplot.compute(npts=50) # nprocess=1)
#pplot.plot(ax=ax, color="xkcd:dark grey", s=1.4, linewidths=0)

#ax.set_aspect([equal])



##top island m=2
left_o1 = FixedPoint(section)
left_x1 = FixedPoint(section)
right_o1 = FixedPoint(section)
right_x1 = FixedPoint(section)


left_o1.find(2, [0.84,0.32], method='scipy.root')
left_x1.find(2, [0.82,0.48], method='scipy.root')
right_o1.find(2, [0.93,0.51], method='scipy.root')
right_x1.find(2, [0.95,0.35], method='scipy.root')


left_o1_coord = left_o1.coords[0]
right_o1_coord = right_o1.coords[0]
left_x1_coord = left_x1.coords[0]
right_x1_coord = right_x1.coords[0] 


left_o1.plot(ax=ax, marker='o', color="xkcd:blue")
left_x1.plot(ax=ax, marker='x', color="xkcd:blue")
right_o1.plot(ax=ax, marker='o', color="xkcd:blue")
right_x1.plot(ax=ax, marker='x', color="xkcd:blue")














#######. island m=3


left_o2 = FixedPoint(section)
left_x2 = FixedPoint(section)
right_o2 = FixedPoint(section)
right_x2 = FixedPoint(section)
top_o2 = FixedPoint(section)
bottom_x2 = FixedPoint(section)

#left_o2.find(3, [0.76,-0.53], method='scipy.root')
#left_x2.find(3, [0.76,-0.33], method='scipy.root')
#right_o2.find(3, [1.025,-0.47], method='scipy.root')
#right_x2.find(3, [0.96,-0.21], method='scipy.root')
#top_o2.find(3, [0.87,-0.2], method='scipy.root')
#bottom_x2.find(3, [0.9,-0.62], method='scipy.root')

#left_o2_coord = left_o2.coords[0]
#right_o2_coord = right_o2.coords[0]
#top_o2_coord = top_o2.coords[0]
#bottom_x2_coord = bottom_x2.coords[0]
#left_x2_coord = left_x2.coords[0]
#right_x2_coord = right_x2.coords[0]

#left_o2.plot(ax=ax, marker='o', color="xkcd:green")
#left_x2.plot(ax=ax, marker='x', color="xkcd:green")
#right_o2.plot(ax=ax, marker='o', color="xkcd:green")
#right_x2.plot(ax=ax, marker='x', color="xkcd:green")
#top_o2.plot(ax=ax, marker='o', color="xkcd:green")
#bottom_x2.plot(ax=ax, marker='x', color="xkcd:green")



#Hits=np.load('Script/doublet_poincare_hits_1.npy')


ax.set_xlabel(r"$R[m]$")
ax.set_ylabel(r"$Z[m]$")

ax.scatter(Hits[:,:, 0], Hits[:,:, 1], color="xkcd:dark grey", s=1.4, linewidths=0)




#perturbedmanifold = Manifold(section, xpoint, xpoint)
#perturbedmanifold.compute(
#   eps_s=9e-6, eps_u=8e-6, nint_s=8, nint_u=8, neps_s=80, neps_u=80)
#perturbedmanifold.plot(ax=ax, markersize=0, lw=1.5)

#np.save('doublet_poincare_hits', pplot._hits)

#plt.savefig('doublet_plot2.png', bbox_inches='tight', dpi=720)
plt.show()

#print(r)