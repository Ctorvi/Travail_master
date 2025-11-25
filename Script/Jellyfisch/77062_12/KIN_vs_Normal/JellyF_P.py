from pyoculus.fields import AxisymmetricCylindricalGridField
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
import numpy as np
from matplotlib import pyplot as plt
import logging
from scipy.io import loadmat
from matplotlib.tri import Triangulation




#fig, ax = plt.subplots(1, 1, figsize=(5, 8))
fig, ax = plt.subplots(1, 1, figsize=(6, 6.6))


plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.size": 9,
    }
)


logging.basicConfig(level=logging.DEBUG)

repository_path = './Script/Jellyfisch/77062_12/KINETIC/'
pert_mat_file = '77062_12_KINETIC.mat'


JFField = AxisymmetricCylindricalGridField.from_matlab_file(f"{repository_path}{pert_mat_file}", with_perturbation=True)


section = CylindricalBfieldSection(JFField,phi0=1.9,R0=0.88, Z0=0)

top_o = FixedPoint(section)
top_o.find(1, [0.9,0.18], method='scipy.root')
top_o_coord = top_o.coords[0]



x_point1= FixedPoint(section)
x_point1.find(1, [0.8,-0.27], method='scipy.root')
x_point1_coord = x_point1.coords[0]



x_point2= FixedPoint(section)
x_point2.find(1, [0.8,-0.57], method='scipy.root')
x_point2_coord = x_point2.coords[0]



x_point3= FixedPoint(section)
x_point3.find(1, [1.05,-0.58], method='scipy.root')
x_point3_coord = x_point3.coords[0]

x_point4= FixedPoint(section)
x_point4.find(1, [0.75,0.65], method='scipy.root')
x_point4_coord = x_point4.coords[0]


##### loading and plotting  ###



#manifold_NK= Manifold.load('./Script/Jellyfisch/77062_12/manifolds_P/mf_1T.pkl')

manifold_1T  = Manifold.load(f"{repository_path}manifolds_P/mf_1T.pkl")
manifold_1B  = Manifold.load(f"{repository_path}manifolds_P/mf_1B.pkl")
manifold_2T  = Manifold.load(f"{repository_path}manifolds_P/mf_2T.pkl")
manifold_2B  = Manifold.load(f"{repository_path}manifolds_P/mf_2B.pkl")
manifold_3L = Manifold.load(f"{repository_path}manifolds_P/mf_3L.pkl")
manifold_3R  = Manifold.load(f"{repository_path}manifolds_P/mf_3R.pkl")
manifold_4T = Manifold.load(f"{repository_path}manifolds_P/mf_4T.pkl")
manifold_4B  = Manifold.load(f"{repository_path}manifolds_P/mf_4B.pkl")

#manifold_NK.plot(stepsize_limit=0.1, ax=ax, markersize=0, lw=0.7,colors=["grey", "xkcd:blue"])

manifold_1T.plot(stepsize_limit=0.1, ax=ax, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"],labels=['Kin stable','Kin unstable'])
manifold_1B.plot(stepsize_limit=0.3,ax=ax, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"],labels=[None,None])
manifold_2T.plot(stepsize_limit=0.1, ax=ax, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"],labels=[None,None])
manifold_2B.plot(stepsize_limit=0.3, ax=ax, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"],labels=[None,None])
manifold_3L.plot(stepsize_limit=0.1,ax=ax, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"],labels=[None,None])
manifold_3R.plot(ax=ax, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"],labels=[None,None])
manifold_4T.plot(ax=ax, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"],labels=[None,None])
manifold_4B.plot(ax=ax, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"],labels=[None,None])



top_o.plot(ax=ax, marker='o', color="red",label=None)  
x_point1.plot(ax=ax, marker='x', color="red",label=None)
x_point2.plot(ax=ax, marker='x', color="red",label='Kin fp')
x_point3.plot(ax=ax, marker='x', color="red",label=None)
x_point4.plot(ax=ax, marker='x', color="red",label=None)
ratio=2.8158








################ NORMAL





repository_path_N = './Script/Jellyfisch/77062_12/normal/'
pert_mat_file_N = 'JF_77062_12.mat'


JFField = AxisymmetricCylindricalGridField.from_matlab_file(f"{repository_path_N}{pert_mat_file_N}", with_perturbation=True)


section = CylindricalBfieldSection(JFField,phi0=1.9,R0=0.88, Z0=0)

top_o_N = FixedPoint(section)
top_o_N.find(1, [0.9,0.18], method='scipy.root')
top_o_N_coord = top_o_N.coords[0]



x_point1_N = FixedPoint(section)
x_point1_N.find(1, [0.8,-0.27], method='scipy.root')
x_point1_N_coord = x_point1_N.coords[0]



x_point2_N = FixedPoint(section)
x_point2_N.find(1, [0.8,-0.57], method='scipy.root')
x_point2_N_coord = x_point2_N.coords[0]



x_point3_N = FixedPoint(section)
x_point3_N.find(1, [1.05,-0.58], method='scipy.root')
x_point3_N_coord = x_point3_N.coords[0]

x_point4_N = FixedPoint(section)
x_point4_N.find(1, [0.75,0.65], method='scipy.root')
x_point4_N_coord = x_point4_N.coords[0]


##### loading and plotting  ###



#manifold_NK= Manifold.load('./Script/Jellyfisch/77062_12/manifolds_P/mf_1T.pkl')

manifold_1T_N  = Manifold.load(f"{repository_path_N}manifolds_P/mf_1T.pkl")
manifold_1B_N  = Manifold.load(f"{repository_path_N}manifolds_P/mf_1B.pkl")
manifold_2T_N  = Manifold.load(f"{repository_path_N}manifolds_P/mf_2T.pkl")
manifold_2B_N  = Manifold.load(f"{repository_path_N}manifolds_P/mf_2B.pkl")
manifold_3L_N = Manifold.load(f"{repository_path_N}manifolds_P/mf_3L.pkl")
manifold_3R_N  = Manifold.load(f"{repository_path_N}manifolds_P/mf_3R.pkl")
manifold_4T_N = Manifold.load(f"{repository_path_N}manifolds_P/mf_4T.pkl")
manifold_4B_N  = Manifold.load(f"{repository_path_N}manifolds_P/mf_4B.pkl")

#manifold_NK.plot(stepsize_limit=0.1, ax=ax, markersize=0, lw=0.7,colors=["grey", "xkcd:blue"])

manifold_1T_N.plot(stepsize_limit=0.1, ax=ax, markersize=0, lw=0.7,colors=["cyan", "xkcd:blue"],labels=['std stable','std unstable'])
manifold_1B_N.plot(stepsize_limit=0.3,ax=ax, markersize=0, lw=0.7,colors=["cyan", "xkcd:blue"],labels=[None,None])
manifold_2T_N.plot(stepsize_limit=0.1, ax=ax, markersize=0, lw=0.7,colors=["cyan", "xkcd:blue"],labels=[None,None])
manifold_2B_N.plot(stepsize_limit=0.3, ax=ax, markersize=0, lw=0.7,colors=["cyan", "xkcd:blue"],labels=[None,None])
manifold_3L_N.plot(stepsize_limit=0.1,ax=ax, markersize=0, lw=0.7,colors=["cyan", "xkcd:blue"],labels=[None,None])
manifold_3R_N.plot(ax=ax, markersize=0, lw=0.7,colors=["cyan", "xkcd:blue"],labels=[None,None])
manifold_4T_N.plot(ax=ax, markersize=0, lw=0.7,colors=["cyan", "xkcd:blue"],labels=[None,None])
manifold_4B_N.plot(ax=ax, markersize=0, lw=0.7,colors=["cyan", "xkcd:blue"],labels=[None,None])




top_o_N.plot(ax=ax, marker='o', color="blue",label=None)
x_point1_N.plot(ax=ax, marker='x', color="blue",label='std fp')
x_point2_N.plot(ax=ax, marker='x', color="blue",label=None)
x_point3_N.plot(ax=ax, marker='x', color="blue",label=None)
x_point4_N.plot(ax=ax, marker='x', color="blue",label=None)
ratio=2.8158






#####figure


# ax.set_xlim(0.62, 1.15)
# ax.set_ylim(-0.75, 0.75)
# ax.set_xlabel(r"$R[m]$")
# ax.set_ylabel(r"$Z[m]$")
# ax.set_aspect('equal') 
# ax.legend()
# ax.set_title('KINETIC 77062 at 1.2s, at 1.9 rad')
# #plt.savefig(f"{repository_path}figures/Jellyfisch_77062_KINETIC_BO_Pert_mf_patch.png", bbox_inches='tight', dpi=720)
# plt.show()


#####zoomed



ax.set_xlim(0.7, 1.05)
ax.set_ylim(-0.6, -0.1)
ax.set_xlabel(r"$R[m]$")
ax.set_ylabel(r"$Z[m]$")
ax.set_aspect('equal') 
ax.legend()
ax.set_title('comparison LIUQE 77062 at 1.2s')
plt.savefig('./Script/Jellyfisch/77062_12/KIN_vs_Normal/figures/Jellyfisch_77062_comparison.png', bbox_inches='tight', dpi=720)
plt.show()





