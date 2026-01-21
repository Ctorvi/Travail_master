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

DoubletField = AxisymmetricCylindricalGridField.from_matlab_file('./Script/Doublet/struct_for_chris.mat', with_perturbation=False)

section = CylindricalBfieldSection(DoubletField, R0=0.88, Z0=0,domain=[(-np.inf, np.inf), (-np.inf, np.inf)], rtol=1e-10)

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(5, 8))


guess_t = [0.8852,0.42]
top_o = FixedPoint(section)
top_o.find(1, guess_t, method='scipy.root')
top_o_coord = top_o.coords[0]
top_o.plot(ax=ax, color="xkcd:dark blue",label='O-points', zorder=10)


guess_b = [0.885, -0.42]
bottom_o = FixedPoint(section)
bottom_o.find(1, guess_b,  method='scipy.root')
bottom_o_coord = bottom_o.coords[0]
bottom_o.plot(ax=ax, color="xkcd:dark blue",label=None, zorder=10)


xpoint= FixedPoint(section)
xpoint.find(1, [0.8865, -0.0001], method='scipy.root')
x_point_coord = xpoint.coords[0]
xpoint.plot(ax=ax, color="xkcd:dark blue",label='X-points', zorder=10)


pplot = PoincarePlot.with_segments(section, [top_o_coord, x_point_coord, bottom_o_coord, x_point_coord], [10, 10], connected=False)
#pplot.compute(npts=100) # nprocess=1)
#pplot.plot_q(ax=ax2,xlim=(x_point_coord[0], bottom_o_coord[0]), ylim=(x_point_coord[1], bottom_o_coord[1]))




# np.save('script/doublet/completePlot/doublet_poincare_XS_NP', pplot.xs)
# np.save('script/doublet/completePlot/doublet_poincare_IOTA_NP', pplot._iota)



xs=np.load('script/doublet/completePlot/doublet_poincare_XS_NP.npy')
iota=np.load('script/doublet/completePlot/doublet_poincare_IOTA_NP.npy')

distance_to_map_axis = xs - np.array([section.R0, section.Z0])

rho=np.linalg.norm(distance_to_map_axis, axis=1)
        

Hits=np.load('script/doublet/completePlot/doublet_poincare_hits_NP.npy')
ax.scatter(Hits[:,:, 0], Hits[:,:, 1], color="xkcd:dark grey", s=0.5, linewidths=0)

ax.set_xlabel(r"$R[m]$")
ax.set_ylabel(r"$Z[m]$")
ax.set_aspect('equal')

# manifold1 = Manifold(section, xpoint, xpoint, top_o_coord - x_point_coord, top_o_coord - x_point_coord)
# manifold1.compute(
#   eps_s=9e-6, eps_u=8e-6, nint_s=8, nint_u=8, neps_s=80, neps_u=80)
# manifold1.save('script/doublet/completePlot/doublet_mf_NP_top.pkl')
manifold1  = Manifold.load('script/doublet/completePlot/doublet_mf_NP_top.pkl')
manifold1.plot(ax=ax, markersize=0, lw=0.6)



# manifold2 = Manifold(section, xpoint, xpoint,bottom_o_coord - x_point_coord, bottom_o_coord - x_point_coord)
# manifold2.compute(
#   eps_s=9e-6, eps_u=8e-6, nint_s=8, nint_u=8, neps_s=80, neps_u=80)
# manifold2.save('script/doublet/completePlot/doublet_mf_NP_bottom.pkl')
manifold2  = Manifold.load('script/doublet/completePlot/doublet_mf_NP_bottom.pkl')
manifold2.plot(ax=ax, markersize=0, lw=0.6)



ax2.scatter(rho, iota/(2*np.pi), color="xkcd:dark red")

ax.set_xlim(0.62, 1.15)
ax.set_ylim(-0.75, 0.75)
ax.set_title('doublet config')


#plt.savefig('doublet_plot_NP.png', bbox_inches='tight', dpi=720)
plt.show()

#print(r)