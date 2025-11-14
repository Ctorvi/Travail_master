from pyoculus.fields import AxisymmetricCylindricalGridField
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
import numpy as np
from matplotlib import pyplot as plt
import logging
from scipy.io import loadmat
from matplotlib.tri import Triangulation


fig, ax = plt.subplots(1, 1, figsize=(5, 8))


plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.size": 10,
    }
)


logging.basicConfig(level=logging.DEBUG)


JFField = AxisymmetricCylindricalGridField.from_matlab_file('./Script/Jellyfisch/PC_analysis/OH coil/pert files/OH_DX_plus.mat', with_perturbation=True)


section = CylindricalBfieldSection(JFField,phi0=1.9,R0=0.88, Z0=0)

top_o = FixedPoint(section)
top_o.find(1, [0.9,0.18], method='scipy.root')
top_o_coord = top_o.coords[0]
top_o.plot(ax=ax, marker='o', color="xkcd:crimson")


x_point1= FixedPoint(section)
x_point1.find(1, [0.8,-0.27], method='scipy.root')
x_point1_coord = x_point1.coords[0]
x_point1.plot(ax=ax, marker='x', color="xkcd:crimson")


x_point2= FixedPoint(section)
x_point2.find(1, [0.8,-0.57], method='scipy.root')
x_point2_coord = x_point2.coords[0]
x_point2.plot(ax=ax, marker='x', color="xkcd:crimson")





#####save #########


manifold = Manifold(section, x_point1, x_point1,-x_point1_coord+top_o_coord, -x_point1_coord+top_o_coord)
manifold.compute(
   eps_s=9e-6, eps_u=8e-6, nint_s=8, nint_u=14, neps_s=80, neps_u=240) #8 14 240

manifold.save('./Script/Jellyfisch/PC_analysis/OH coil/manifolds/mf_OH_2.pkl')


 



### OH coil only

#manifold  = Manifold.load("./Script/Jellyfisch/PC_analysis/OH coil/manifolds/mf_OH_1.pkl") ##OY tilt #########
#manifold  = Manifold.load("./Script/Jellyfisch/PC_analysis/OH coil/manifolds/mf_OH_2.pkl") ##DX plus#########
#manifold  = Manifold.load("./Script/Jellyfisch/PC_analysis/OH coil/manifolds/mf_OH_3.pkl") ##OY tilt DX plus #########
#manifold  = Manifold.load("./Script/Jellyfisch/PC_analysis/OH coil/manifolds/mf_OH_4.pkl") ##OY tilt DY plus#########
#manifold  = Manifold.load("./Script/Jellyfisch/PC_analysis/OH coil/manifolds/mf_OH_5.pkl") ##OY tilt DX neg#########
#manifold  = Manifold.load("./Script/Jellyfisch/PC_analysis/OH coil/manifolds/mf_OH_6.pkl") ##OY tilt DY neg #########



manifold1  = Manifold.load("./Script/Jellyfisch/77021/77021_120/backoff/manifolds_P/mf_1T_NA.pkl") ##with DX all #########

####################
####################


manifold.plot(stepsize_limit=0.2, ax=ax, markersize=0, lw=0.7)
manifold1.plot(stepsize_limit=0.2, ax=ax, markersize=0, lw=0.7, colors=["yellow", "xkcd:cyan"])

ratio=2.8158






ax.set_xlim(0.62, 1.15)
ax.set_ylim(-0.75, 0.75)
ax.set_xlabel(r"$R[m]$")
ax.set_ylabel(r"$Z[m]$")
ax.set_aspect('equal') 
ax.set_title('DX plus - OH coil only')
plt.savefig('./Script/Jellyfisch/PC_analysis/OH coil/figures/OH_DX_plus.png', bbox_inches='tight', dpi=720)
plt.show()





