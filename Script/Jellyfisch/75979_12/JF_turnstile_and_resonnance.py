import os, sys
proj_root = os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..','..'))  # remonte jusqu'à "Travail master"
if proj_root not in sys.path:
    sys.path.insert(0, proj_root)
from script.function.field_utils import plot_tomographic as tomo
from script.function.field_utils import plot_LIUQE as LIUQE


from pyoculus.fields import AxisymmetricCylindricalGridField
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
import numpy as np
from matplotlib import pyplot as plt
import logging
from scipy.io import loadmat
from matplotlib.lines import Line2D
from matplotlib.tri import Triangulation

#fig, ax2 = plt.subplots(1, 1, figsize=(5, 8))
fig, ax2 = plt.subplots(1, 1, figsize=(20, 20))

plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.size": 9,
    }
)


file_manifolds = 'manifolds_P/OS_BO78/'


logging.basicConfig(level=logging.DEBUG)

repository_path = './script/jellyfisch/75979_12/'

pert_mat_file = 'mat_files/JF_75979_120_BO78_OS.mat'

mf_list = ['mf_1T_high_RES', 'mf_1B', 'mf_2T', 'mf_2B', 'mf_3L', 'mf_3R', 'mf_4T', 'mf_4B']

patch_mat_file = './script/jellyfisch/75979_12/mat_files/JF_75979_patch_120_He2.mat'

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


###top fp top mf#########

# manifold_1T = Manifold(section, x_point1, x_point1,-x_point1_coord+top_o_coord, -x_point1_coord+top_o_coord)
# manifold_1T.compute(
#    eps_s=9e-6, eps_u=8e-6, nint_s=30, nint_u=22, neps_s=480, neps_u=480) #8 14 240

# manifold_1T.save(f"{repository_path}{file_manifolds}mf_1T_high_RES.pkl")


# #####top fp bottom mf#########

# manifold_1B = Manifold(section, x_point1, x_point1)
# manifold_1B.compute(
#        eps_s=9e-6, eps_u=8e-6, nint_s=17, nint_u=18, neps_s=80, neps_u=80)
# manifold_1B.save(f"{repository_path}{file_manifolds}mf_1B.pkl")


# #### bottom fp top mf


# manifold_2T = Manifold(section, x_point2, x_point2,-x_point2_coord+top_o_coord, -x_point2_coord+top_o_coord)
# manifold_2T.compute(
#       eps_s=9e-6, eps_u=8e-6, nint_s=32, nint_u=30, neps_s=80, neps_u=80)
# manifold_2T.save(f"{repository_path}{file_manifolds}mf_2T.pkl")

# #### bottom fp bottom mf


# manifold_2B = Manifold(section, x_point2, x_point2,x_point2_coord-top_o_coord, x_point2_coord-top_o_coord)
# manifold_2B.compute(
#       eps_s=9e-6, eps_u=8e-6, nint_s=20, nint_u=20, neps_s=80, neps_u=80)
# manifold_2B.save(f"{repository_path}{file_manifolds}mf_2B.pkl")

# #### right fp left mf

# manifold_3L = Manifold(section, x_point3, x_point3,-x_point3_coord+top_o_coord, -x_point3_coord+top_o_coord)
# manifold_3L.compute(
#       eps_s=9e-6, eps_u=8e-6, nint_s=33, nint_u=30, neps_s=80, neps_u=80)
# manifold_3L.save(f"{repository_path}{file_manifolds}mf_3L.pkl")

# # ##### right fp right mf

# manifold_3R = Manifold(section, x_point3, x_point3,x_point3_coord-top_o_coord, x_point3_coord-top_o_coord)
# manifold_3R.compute(
#       eps_s=9e-6, eps_u=8e-6, nint_s=20, nint_u=30, neps_s=80, neps_u=80)
# manifold_3R.save(f"{repository_path}{file_manifolds}mf_3R.pkl")

# ## top fp top mf 

# manifold_4T = Manifold(section, x_point4, x_point4,x_point4_coord-top_o_coord, x_point4_coord-top_o_coord)
# manifold_4T.compute(
#      eps_s=9e-6, eps_u=8e-6, nint_s=8, nint_u=8, neps_s=80, neps_u=80)
# manifold_4T.save(f"{repository_path}{file_manifolds}mf_4T.pkl")

# # # # ##### over the top fp bottom mf


# manifold_4B = Manifold(section, x_point4, x_point4,-x_point4_coord+top_o_coord, -x_point4_coord+top_o_coord)
# manifold_4B.compute(
#      eps_s=9e-6, eps_u=8e-6, nint_s=8, nint_u=30, neps_s=80, neps_u=80)
# manifold_4B.save(f"{repository_path}{file_manifolds}mf_4B.pkl")

###################.  poincaré plot #####################

# pplot = PoincarePlot.with_linspace(section, top_o_coord, x_point1_coord,80)
# pplot.compute(400)
# np.save('./script/jellyfisch/75979_12/JF_75979_pp_hits_n80.npy', pplot._hits)

################### tomographic plot #####################

#tpc= tomo(f"{patch_mat_file}",ax2,emi_vmin=0, emi_vmax=2.25e19)


##########################


Hits=np.load('./script/jellyfisch/75979_12/JF_75979_pp_hits_n80.npy')
ax2.scatter(Hits[:,:, 0], Hits[:,:, 1], color="xkcd:dark grey", s=0.8, linewidths=0)
# ax1.scatter(Hits[:,:, 0], Hits[:,:, 1], color="xkcd:dark grey", s=0.6, linewidths=0)
eps_s1 = 8.2e-07
eps_u1 = 7.4e-07


manifs = {name: {name} for name in mf_list}

for i in mf_list:
     manifs[i] = Manifold.load(f"{repository_path}{file_manifolds}{i}.pkl")
    # manifs[i].plot(ax=ax2, markersize=0, lw=0.7,labels=[None,None] if i!=0 else ['stable MF','unstable MF'])
    # if i=='mf_1T_high_RES':
          #manifs[i].find_clinics(first_guess_eps_s=eps_s1, first_guess_eps_u=eps_u1,n_points=2)
          #A=manifs[i].compute_turnstile_areas()
          
          #fig, ax2=manifs[i].plot_clinics(ax=ax2,label='homoclinics',s=15,color='xkcd:dark blue')
          #manifs[i].plot_filled_lobe(ax=ax2, lobe_number=21,which_section=2,alpha=0.7,color='red')
          #manifs[i].plot_filled_lobe(ax=ax2, lobe_number=13,which_section=2,alpha=0.7,color='blue')
          #manifs[i].save(f"{repository_path}{file_manifolds}{i}.pkl")
     manifs[i].plot(stepsize_limit=0.05, ax=ax2, markersize=0, lw=0.5,colors=["xkcd:royal blue", "xkcd:red"],labels=['stable MF','unstable MF'] if i=='mf_1T' else [None,None])
#stepsize_limit=0.05

top_o.plot(ax=ax2, marker='o',s=40, color="xkcd:dark blue",label=None, zorder=10)
x_point1.plot(ax=ax2, marker='x', s=80, color="xkcd:dark blue",label='X-points', zorder=10)
x_point2.plot(ax=ax2, marker='x', s=50, color="xkcd:dark blue",label=None, zorder=10)
x_point3.plot(ax=ax2, marker='x', s=50, color="xkcd:dark blue",label=None, zorder=10)
x_point4.plot(ax=ax2, marker='x', s=50, color="xkcd:dark blue",label=None, zorder=10)


# top_o.plot(ax=ax2, marker='o', s=50, color="xkcd:green",label='O-point', zorder=10)
# x_point1.plot(ax=ax2, marker='x', s=70, color="xkcd:green",label='X-points', zorder=10)
# x_point2.plot(ax=ax2, marker='x', s=50, color="xkcd:green",label=None, zorder=10)
# x_point3.plot(ax=ax2, marker='x', s=50, color="xkcd:green",label=None, zorder=10)
# x_point4.plot(ax=ax2, marker='x', s=50, color="xkcd:green",label=None, zorder=10)
# ratio=2.8158

# ax2.set_xlim(0.62, 1.15)
# ax2.set_ylim(-0.45, -0.0)

ax2.set_xlim(0.62, 1.15)
ax2.set_ylim(-0.75, 0.75)



ax2.set_xlabel(r"$R[m]$")
ax2.set_ylabel(r"$Z[m]$")
ax2.set_aspect('equal') 


# x_text, y_text = 0.89, -0.4
# # safe extraction de la valeur A[0] et formatage en mathtext (3 décimales)
# val = float(A[1]) #if (isinstance(A, (list, tuple, np.ndarray)) and len(A) > 0) else float(A)
# label_text = rf"$\phi_{{turn}}={val:.4f}$"

# ax2.text(
#     x_text-0.01, y_text-0.025, label_text,
#     color='black', fontsize=10, fontstyle='italic',  # black on (white) background, italic via fontstyle
#     ha='left', va='center',
#     transform=ax2.transData,
#     zorder=200
# )
# label_text1 = rf"$f^{8}(E)$"
# ax2.text(
#     x_text, y_text, label_text1,
#     color='red', fontsize=13, fontstyle='italic',  # black on (white) background, italic via fontstyle
#     ha='left', va='center',
#     transform=ax2.transData,
#     zorder=200
# )

# label_text2 = rf"$E$"
# ax2.text(
#     0.78, -0.25, label_text2,
#     color='blue', fontsize=13, fontstyle='italic',  # black on (white) background, italic via fontstyle
#     ha='left', va='center',
#     transform=ax2.transData,
#     zorder=200
# )




#ax2.set_title('Poincaré perturbed')

# ax2.legend(loc='lower right', bbox_to_anchor=(1., 0.0), ncol=2, fontsize=9)

plt.savefig(f"{repository_path}figures/JF_75979_turnstile_ORAL.png",  bbox_inches="tight",dpi=720, pad_inches=0, facecolor=fig.get_facecolor())
plt.show()
# print(f"Turnstile area : {A} m^2")




########### resonnance zone #############



# top_o = FixedPoint(section)
# top_o.find(1, [0.9,0.18], method='scipy.root')
# top_o_coord = top_o.coords[0]


# o_point_island= FixedPoint(section)
# o_point_island.find(1, [0.99,0.199], method='scipy.root')
# o_point_island_coord = o_point_island.coords[0]

# # x_point_island_top= FixedPoint(section)
# # x_point_island_top.find(1, [0.85,0.3], method='scipy.root')
# # x_point_island_top_coord = x_point_island_top.coords[0]

# # x_point_island_bottom= FixedPoint(section)
# # x_point_island_bottom.find(1, [0.85,0.06], method='scipy.root')
# # x_point_island_bottom_coord = x_point_island_bottom.coords[0]



# x_point_island_n= FixedPoint(section)
# x_point_island_n.find(1, [0.79,0.18], method='scipy.root')
# x_point_island_n_coord = x_point_island_n.coords[0]


# o_point_island.plot(ax=ax2, marker='o',s=40, color="xkcd:dark blue",label='O-points', zorder=10)
# # x_point_island_top.plot(ax=ax2, marker='x', s=50, color="xkcd:dark blue",label='X-points', zorder=10)
# x_point_island_n.plot(ax=ax2, marker='x', s=50, color="xkcd:dark blue",label=None, zorder=10)
# #x_point_island_bottom.plot(ax=ax2, marker='x', s=50, color="xkcd:dark blue",label=None, zorder=10)


# #####top fp bottom mf#########

# manifold_1T = Manifold(section, x_point_island_n, x_point_island_n, [0,1],[0,1] )
# manifold_1T.compute(
#        eps_s=9e-6, eps_u=8e-6, nint_s=80, nint_u=80, neps_s=40, neps_u=40)



# manifold_1T = Manifold.load(f"{repository_path}{file_manifolds}mf_1_island_n.pkl")


# manifold_1T.find_clinics(first_guess_eps_s=8e-6, first_guess_eps_u=9e-8)


# # manifold_1T.save(f"{repository_path}{file_manifolds}mf_1_island_n.pkl")



# manifold_1T.plot(ax=ax2, stepsize_limit=0.1, linewidth=0.7, markersize=0)

# # manifold_1B = Manifold(section, x_point_island_n, x_point_island_n, [0,-1],[0,-1] )
# # manifold_1B.compute(
# #        eps_s=9e-6, eps_u=8e-6, nint_s=80, nint_u=80, neps_s=40, neps_u=40)

# manifold_1B = Manifold.load(f"{repository_path}{file_manifolds}mf_1_island_n2.pkl")

# manifold_1B.find_clinics(first_guess_eps_s=9e-5, first_guess_eps_u=8e-5)




# # manifold_1B.save(f"{repository_path}{file_manifolds}mf_1_island_n2.pkl")

# manifold_1B.plot(ax=ax2, stepsize_limit=0.1, linewidth=0.7, markersize=0)

# manifold_1B = Manifold(section, x_point_island_bottom, x_point_island_bottom,-x_point_island_bottom.coords[0]+o_point_island.coords[0], -x_point_island_bottom.coords[0]+o_point_island.coords[0])
# manifold_1B.compute(
#        eps_s=9e-6, eps_u=8e-6, nint_s=17, nint_u=18, neps_s=80, neps_u=80)
# manifold_1B.save(f"{repository_path}{file_manifolds}mf_2_island.pkl")
# manifold_1B.plot(ax=ax2, stepsize_limit=0.1, linewidth=0.7, markersize=0, color="xkcd:red", label='manifold bottom FP')



# ax2.set_xlim(0.75, 1.0)
# ax2.set_ylim(0.0, 0.35)
# ax2.set_xlabel(r"$R[m]$")
# ax2.set_ylabel(r"$Z[m]$")
# ax2.set_aspect('equal') 
# ax2.set_title('Poincaré perturbed')

# ax2.legend(loc='lower right', bbox_to_anchor=(1., 0.86), ncol=1, fontsize=9)

#plt.savefig(f"{repository_path}figures/JF_75979_P_pp.png", bbox_inches='tight', dpi=720)
# plt.show()


