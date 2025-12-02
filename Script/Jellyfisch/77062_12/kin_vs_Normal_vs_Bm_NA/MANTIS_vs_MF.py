from pyoculus.fields import AxisymmetricCylindricalGridField
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
import numpy as np
from matplotlib import pyplot as plt
import logging
from scipy.io import loadmat
from matplotlib.tri import Triangulation




#fig, ax = plt.subplots(2, 1, figsize=(5, 8))
#fig, ax = plt.subplots(1, 1, figsize=(6, 6.6))
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 8), gridspec_kw={'width_ratios': [1, 1]})

plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.size": 9,
    }
)


logging.basicConfig(level=logging.DEBUG)

#repository_path = './script/jellyfisch/77062_12//'

pert_mat_file_KIN = './script/jellyfisch/77062_12/kinetic/JF_77062_120_KINETIC.mat'
pert_mat_file_norm = './script/jellyfisch/77062_12/normal/JF_77062_120.mat'
pert_mat_file_Bm_NA = './script/jellyfisch/77062_12/bm_corrected/JF_77062_120_Bm_corrected.mat'

patch_mat_file = './script/jellyfisch/77062_12/JF_77062_patch_12_C3.mat'

JFField_KIN = AxisymmetricCylindricalGridField.from_matlab_file(f"{pert_mat_file_KIN}", with_perturbation=True)
JFField_norm = AxisymmetricCylindricalGridField.from_matlab_file(f"{pert_mat_file_norm}", with_perturbation=True)
JFField_Bm_NA = AxisymmetricCylindricalGridField.from_matlab_file(f"{pert_mat_file_Bm_NA}", with_perturbation=True)

section_KIN = CylindricalBfieldSection(JFField_KIN,phi0=1.9,R0=0.88, Z0=0)
section_norm = CylindricalBfieldSection(JFField_norm,phi0=1.9,R0=0.88, Z0=0)
section_Bm_NA = CylindricalBfieldSection(JFField_Bm_NA,phi0=1.9,R0=0.88, Z0=0)





############ fp dec #######################



top_o_B= FixedPoint(section_Bm_NA)
top_o_B.find(1, [0.9,0.18], method='scipy.root')
#top_o_coord = top_o.coords[0]
 

x_point1_B= FixedPoint(section_Bm_NA)
x_point1_B.find(1, [0.8,-0.27], method='scipy.root')
#x_point1_coord = x_point1.coords[0]

x_point2_B= FixedPoint(section_Bm_NA)
x_point2_B.find(1, [0.8,-0.57], method='scipy.root')
#x_point2_coord = x_point2.coords[0]

x_point3_B= FixedPoint(section_Bm_NA)
x_point3_B.find(1, [1.05,-0.58], method='scipy.root')
#x_point3_coord = x_point3.coords[0]

x_point4_B= FixedPoint(section_Bm_NA)
x_point4_B.find(1, [0.75,0.65], method='scipy.root')
#x_point4_coord = x_point4.coords[0]




############ fp inc #######################



top_o_K = FixedPoint(section_KIN)
top_o_K.find(1, [0.9,0.18], method='scipy.root')
#top_o_coord = top_o_i.coords[0]

x_point1_K= FixedPoint(section_KIN)
x_point1_K.find(1, [0.8,-0.27], method='scipy.root')
#x_point1_coord = x_point1_i.coords[0]

x_point2_K= FixedPoint(section_KIN)
x_point2_K.find(1, [0.8,-0.57], method='scipy.root')
#x_point2_coord = x_point2_i.coords[0]

x_point3_K= FixedPoint(section_KIN)
x_point3_K.find(1, [1.05,-0.58], method='scipy.root')
#x_point3_coord = x_point3_i.coords[0]

x_point4_K= FixedPoint(section_KIN)
x_point4_K.find(1, [0.75,0.65], method='scipy.root')
#x_point4_coord = x_point4_i.coords[0]




###########################  tomographic reconstruction  #########################





data = loadmat(f"{patch_mat_file}", squeeze_me=True, struct_as_record=False)

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
tpc=ax1.tripcolor(tri, emi, shading='flat', edgecolors='none', vmin=0, vmax=8e20)
tpc=ax2.tripcolor(tri, emi, shading='flat', edgecolors='none', vmin=0, vmax=8e20)
#fig.colorbar(tpc, ax=ax, label='Émission')
   




##### loading and plotting  ###

#########################################.  kinetic  ############################


manifold_1T_K  = Manifold.load('./script/jellyfisch/77062_12/kinetic/manifolds_P/mf_1T.pkl')
manifold_1B_K  = Manifold.load('./script/jellyfisch/77062_12/kinetic/manifolds_P/mf_1B.pkl')
manifold_2T_K  = Manifold.load('./script/jellyfisch/77062_12/kinetic/manifolds_P/mf_2T.pkl')
manifold_2B_K  = Manifold.load('./script/jellyfisch/77062_12/kinetic/manifolds_P/mf_2B.pkl')
manifold_3L_K  = Manifold.load('./script/jellyfisch/77062_12/kinetic/manifolds_P/mf_3L.pkl')
manifold_3R_K  = Manifold.load('./script/jellyfisch/77062_12/kinetic/manifolds_P/mf_3R.pkl')
manifold_4T_K  = Manifold.load('./script/jellyfisch/77062_12/kinetic/manifolds_P/mf_4T.pkl')
manifold_4B_K  = Manifold.load('./script/jellyfisch/77062_12/kinetic/manifolds_P/mf_4B.pkl')

manifold_1T_K.plot(stepsize_limit=0.1, ax=ax1, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"],labels=['kin stable','kin unstable'])
manifold_1B_K.plot(stepsize_limit=0.3,ax=ax1, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"],labels=[None,None])
manifold_2T_K.plot(stepsize_limit=0.1, ax=ax1, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"],labels=[None,None])
manifold_2B_K.plot(stepsize_limit=0.3, ax=ax1, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"],labels=[None,None])
manifold_3L_K.plot(stepsize_limit=0.1,ax=ax1, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"],labels=[None,None])
manifold_3R_K.plot(ax=ax1, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"],labels=[None,None])
manifold_4T_K.plot(ax=ax1, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"],labels=[None,None])
manifold_4B_K.plot(ax=ax1, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"],labels=[None,None])

#########################################.  Bm corrected  ############################

manifold_1T_B  = Manifold.load('./script/jellyfisch/77062_12/bm_corrected/manifolds_P/mf_1T.pkl') 
manifold_1B_B  = Manifold.load('./script/jellyfisch/77062_12/bm_corrected/manifolds_P/mf_1B.pkl')
manifold_2T_B  = Manifold.load('./script/jellyfisch/77062_12/bm_corrected/manifolds_P/mf_2T.pkl')
manifold_2B_B  = Manifold.load('./script/jellyfisch/77062_12/bm_corrected/manifolds_P/mf_2B.pkl')
manifold_3L_B  = Manifold.load('./script/jellyfisch/77062_12/bm_corrected/manifolds_P/mf_3L.pkl')
manifold_3R_B  = Manifold.load('./script/jellyfisch/77062_12/bm_corrected/manifolds_P/mf_3R.pkl')
manifold_4T_B  = Manifold.load('./script/jellyfisch/77062_12/bm_corrected/manifolds_P/mf_4T.pkl')
manifold_4B_B  = Manifold.load('./script/jellyfisch/77062_12/bm_corrected/manifolds_P/mf_4B.pkl')

manifold_1T_B.plot(stepsize_limit=0.1, ax=ax2, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"],labels=['Bm corr stable','Bm corr unstable'])
manifold_1B_B.plot(stepsize_limit=0.3,ax=ax2, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"],labels=[None,None])
manifold_2T_B.plot(stepsize_limit=0.1, ax=ax2, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"],labels=[None,None])
manifold_2B_B.plot(stepsize_limit=0.3, ax=ax2, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"],labels=[None,None])
manifold_3L_B.plot(stepsize_limit=0.1,ax=ax2, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"],labels=[None,None])
manifold_3R_B.plot(ax=ax2, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"],labels=[None,None])
manifold_4T_B.plot(ax=ax2, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"],labels=[None,None])
manifold_4B_B.plot(ax=ax2, markersize=0, lw=0.7,colors=["rosybrown", "xkcd:red"],labels=[None,None])






top_o_K.plot(ax=ax1, marker='o', color="xkcd:white",label=None,zorder=10)
x_point1_K.plot(ax=ax1, marker='s', color="xkcd:white",label='X-p kin',zorder=10)
x_point2_K.plot(ax=ax1, marker='s', color="xkcd:white",label=None,zorder=10)
x_point3_K.plot(ax=ax1, marker='s', color="xkcd:white",label=None,zorder=10)
x_point4_K.plot(ax=ax1, marker='s', color="xkcd:white",label=None,zorder=10)


top_o_B.plot(ax=ax2, marker='o', color="xkcd:white",label=None,zorder=10)
x_point1_B.plot(ax=ax2, marker='s', color="xkcd:white",label='X-p Bm corr',zorder=10)
x_point2_B.plot(ax=ax2, marker='s', color="xkcd:white",label=None,zorder=10)
x_point3_B.plot(ax=ax2, marker='s', color="xkcd:white",label=None,zorder=10)
x_point4_B.plot(ax=ax2, marker='s', color="xkcd:white",label=None,zorder=10)




#####figure


# ax.set_xlim(0.62, 1.15)
# ax.set_ylim(-0.75, 0.75)
# ax.set_xlabel(r"$R[m]$")
# ax.set_ylabel(r"$Z[m]$")
# ax.set_aspect('equal') 
# ax.set_title('77062 at 1.2s, at 1.9 rad')
# #plt.savefig(f"{repository_path}figures/Jellyfisch_77062_Ip_Bm_dec.png", bbox_inches='tight', dpi=720)
# plt.show()


#####.  zoomed comparison



ax1.set_xlim(0.7, 1.0)
ax1.set_ylim(-0.6, -0.1)
ax1.set_xlabel(r"$R[m]$")
ax1.set_ylabel(r"$Z[m]$")
ax1.set_aspect('equal') 
ax1.legend(loc='center right', bbox_to_anchor=(1.03, 0.5), ncol=1, fontsize=9)
ax1.set_title('kinetic vs normal')

ax2.set_xlim(0.7, 1.0)
ax2.set_ylim(-0.6, -0.1)
ax2.set_xlabel(r"$R[m]$")
#ax2.set_ylabel(r"$Z[m]$")
ax2.set_aspect('equal') 
ax2.legend(loc='center right', bbox_to_anchor=(1., 0.5), ncol=1, fontsize=9)
ax2.set_title('Bm corrected vs normal')




plt.savefig('./script/jellyfisch/77062_12/kin_vs_Normal_vs_Bm_NA/figures/comparison_kin_Bm_corr_MANTIS.png', bbox_inches='tight', dpi=720)
plt.show()





