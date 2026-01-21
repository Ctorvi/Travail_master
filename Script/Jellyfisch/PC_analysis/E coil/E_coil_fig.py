from pyoculus.fields import AxisymmetricCylindricalGridField
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
import numpy as np
from matplotlib import pyplot as plt
import logging
from scipy.io import loadmat
from matplotlib.tri import Triangulation


fig, axes = plt.subplots(2, 4, figsize=(10,6))


plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.size": 9,
    }
)


logging.basicConfig(level=logging.DEBUG)





m1  = Manifold.load("./script/jellyfisch/PC_analysis/E coil/manifolds/mf_E_Xp.pkl") ##DX plus  #########
m2  = Manifold.load("./script/jellyfisch/PC_analysis/E coil/manifolds/mf_E_t.pkl") ##OY tilt  #########
m3  = Manifold.load("./script/jellyfisch/PC_analysis/E coil/manifolds/mf_E_t_Xm.pkl") ##OY tilt DX minus#########
m4  = Manifold.load("./script/jellyfisch/PC_analysis/E coil/manifolds/mf_E_t_Xp.pkl") ##OY tilt DX plus#########
m5  = Manifold.load("./script/jellyfisch/PC_analysis/E coil/manifolds/mf_E_t_Yp.pkl") ##OY tilt DY plus #########
m6  = Manifold.load("./script/jellyfisch/PC_analysis/E coil/manifolds/mf_E1_to_E4.pkl") ##E1 to E4
m7  = Manifold.load("./script/jellyfisch/PC_analysis/E coil/manifolds/mf_E5_to_E8.pkl") ##E5 to E8
m8  = Manifold.load("./script/jellyfisch/PC_analysis/E coil/manifolds/mf_E4.pkl") ## E4




manifolds=[m1, m2, m3, m4, m5, m6, m7, m8]
title=[r'$\Delta X+$',r'B $\Delta \alpha_{oy}$', r'$\Delta \alpha_{oy}/ \Delta X-$', r'$\Delta \alpha_{oy}/ \Delta X+$', r'$\Delta \alpha_{oy}/ \Delta Y+$',r'E1 to E4 ($\Delta X+$)',r'E5 to E8 ($\Delta X+$)',r'E4 ($\Delta X+$)']
manifold1  = Manifold.load("./script/jellyfisch/77021/1_20/backoff/manifolds_P/OS_BO78/mf_1T_n20.pkl") ##with DX all #########

####################

subscripts = ['(a)', '(b)', '(c)','(d)', '(e)', '(f)','(g)', '(h)']


for i in range(2):
    for j in range(4):
        ax = axes[i, j]
        manifolds[i * 4 + j].plot(stepsize_limit=0.2, which="unstable", ax=ax, markersize=0, lw=0.7,labels=['stable','unstable 77021'])#labels=['stable 77021', 'unstable 77021'])"navajowhite",
        manifold1.plot(stepsize_limit=0.2, which="unstable", ax=ax, markersize=0, lw=0.7, colors=[ 'black', "xkcd:cyan"], labels=['stable','unstable E coil'])#'stable E coil',
        ax.set_xlim(0.7, 1.15)
        ax.set_ylim(-0.6, -0.1)
        ax.set_aspect('equal') 
        ax.set_title(title[i * 4 + j])
        if i == 1:
            ax.set_xlabel(r"$R[m]$")
        else:    
            ax.set_xticklabels([])
        if j == 0:
            ax.set_ylabel(r"$Z[m]$")
        else:
            ax.set_yticklabels([])  
        
        if i == 0 and j == 2:
            ax.legend(['config unst', 'LSE unst.'], loc='lower right', fontsize=7)
        ax.text(
        0.9, 0.9, subscripts[i * 4 + j],
        transform=ax.transAxes,
        fontsize=12,
        fontweight='bold',
        ha='right',
        va='top',
        color='black',
        #bbox=dict(facecolor=bbox_face, edgecolor='none', pad=2, alpha=0.8),
        zorder=200
        )
ratio=2.8158




plt.savefig('./script/jellyfisch/PC_analysis/E coil/figures/Final_E_comparison_v2.png', bbox_inches="tight", pad_inches=0, dpi=720)
plt.show()



