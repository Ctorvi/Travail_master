from pyoculus.fields import AxisymmetricCylindricalGridField
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
import numpy as np
from matplotlib import pyplot as plt
import logging
from scipy.io import loadmat
from matplotlib.tri import Triangulation


fig, axes = plt.subplots(2, 4, figsize=(8,8))


plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.size": 9,
    }
)


logging.basicConfig(level=logging.DEBUG)



m4  = Manifold.load("./Script/Jellyfisch/PC_analysis/manifolds/mf_E_1.pkl") ##OY tilt DX minus#########
m5  = Manifold.load("./Script/Jellyfisch/PC_analysis/manifolds/mf_E_2.pkl") ##OY tilt DX plus#########
m6  = Manifold.load("./Script/Jellyfisch/PC_analysis/manifolds/mf_E_3.pkl") ##OY tilt DY plus #########
m8  = Manifold.load("./Script/Jellyfisch/PC_analysis/manifolds/mf_E_4.pkl") ##OY tilt DX str plus#########
m7  = Manifold.load("./Script/Jellyfisch/PC_analysis/manifolds/mf_E_5.pkl") ##OY tilt DX str neg#########
m9  = Manifold.load("./Script/Jellyfisch/PC_analysis/manifolds/mf_E_6.pkl") ##OY tilt DY str neg #########
m1  = Manifold.load("./Script/Jellyfisch/PC_analysis/manifolds/mf_E_7.pkl") ##OY tilt  #########
m2  = Manifold.load("./Script/Jellyfisch/PC_analysis/manifolds/mf_E_8.pkl") ##DX plus  #########
m3  = Manifold.load("./Script/Jellyfisch/PC_analysis/manifolds/mf_E_9.pkl") ##DX str neg  #########

manifolds=[m1, m2, m4, m5, m3, m7, m9, m6]
title=['OY tilt', 'DX plus', 'OY tilt DX minus', 'OY tilt DX plus', 'DX str neg', 'OY tilt DX str neg', 'OY tilt DY str neg', 'OY tilt DY plus']
manifold1  = Manifold.load("./Script/Jellyfisch/77021/77021_120/backoff/manifolds_P/mf_1T_NA.pkl") ##with DX all #########

####################



for i in range(2):
    for j in range(4):
        ax = axes[i, j]
        manifolds[i * 4 + j].plot(stepsize_limit=0.2, ax=ax, markersize=0, lw=0.7)
        manifold1.plot(stepsize_limit=0.2, ax=ax, markersize=0, lw=0.7, colors=["yellow", "xkcd:cyan"])
        ax.set_xlim(0.62, 1.15)
        ax.set_ylim(-0.6, 0.6)
        ax.set_aspect('equal') 
        ax.set_title(title[i * 4 + j])

        


ratio=2.8158


plt.savefig('./Script/Jellyfisch/PC_analysis/figures/Final_E_comparison.png', bbox_inches='tight', dpi=720)
plt.show()





