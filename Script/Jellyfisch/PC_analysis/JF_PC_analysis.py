from pyoculus.fields import AxisymmetricCylindricalGridField
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
import numpy as np
from matplotlib import pyplot as plt
import logging
from scipy.io import loadmat
from matplotlib.tri import Triangulation

fig, ax = plt.subplots(1, 1, figsize=(5, 6))

plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.size": 11,
    }
)


logging.basicConfig(level=logging.DEBUG)

JFField = AxisymmetricCylindricalGridField.from_matlab_file('./script/jellyfisch/PC_analysis/E coil/pert_files/BO78_no_E4.mat', with_perturbation=True)

section = CylindricalBfieldSection(JFField,phi0=1.9,R0=0.88, Z0=0)

# top_o = FixedPoint(section)
# top_o.find(1, [0.9,0.18], method='scipy.root')
# top_o_coord = top_o.coords[0]
# #top_o.plot(ax=ax, marker='o', color="xkcd:crimson")


# x_point1= FixedPoint(section)
# x_point1.find(1, [0.8,-0.27], method='scipy.root')
# x_point1_coord = x_point1.coords[0]
# #x_point1.plot(ax=ax, marker='x', color="xkcd:crimson")


# x_point2= FixedPoint(section)
# x_point2.find(1, [0.8,-0.57], method='scipy.root')
# x_point2_coord = x_point2.coords[0]
# #x_point2.plot(ax=ax, marker='x', color="xkcd:crimson")


# x_point3= FixedPoint(section)
# x_point3.find(1, [1.05,-0.58], method='scipy.root')
# x_point3_coord = x_point3.coords[0]
# #x_point3.plot(ax=ax, marker='x', color="xkcd:crimson")

# x_point4= FixedPoint(section)
# x_point4.find(1, [0.75,0.65], method='scipy.root')
# x_point4_coord = x_point4.coords[0]
# x_point4.plot(ax=ax, marker='x', color="xkcd:crimson")

# #####save #########

# manifold = Manifold(section, x_point1, x_point1,-x_point1_coord+top_o_coord, -x_point1_coord+top_o_coord)
# manifold.compute(
#    eps_s=9e-6, eps_u=8e-6, nint_s=8, nint_u=14, neps_s=80, neps_u=240) #8 14 240

# manifold.save('./script/jellyfisch/PC_analysis/E coil/manifolds/mf_E_BO78_no_E4.pkl')

###random try

#manifold  = Manifold.load("./Script/Jellyfisch/PC_analysis/manifolds/mf_1.pkl") ##no DX E and OH #########
#manifold  = Manifold.load("./Script/Jellyfisch/PC_analysis/manifolds/mf_2.pkl") ##no tilt E and OH #########
#manifold  = Manifold.load("./Script/Jellyfisch/PC_analysis/manifolds/mf_3.pkl") ##no tilt all #########
#manifold  = Manifold.load("./Script/Jellyfisch/PC_analysis/manifolds/mf_4.pkl") ##no DX all #########
#manifold  = Manifold.load("./Script/Jellyfisch/PC_analysis/manifolds/mf_5.pkl") ##E tilt AMP plus #########
#manifold  = Manifold.load("./Script/Jellyfisch/PC_analysis/manifolds/mf_6.pkl") ##F tilt AMP plus minus #########

### E coil only

#manifold  = Manifold.load("./Script/Jellyfisch/PC_analysis/manifolds/mf_E_t_Xm.pkl") ##OY tilt DX minus#########
#manifold  = Manifold.load("./Script/Jellyfisch/PC_analysis/manifolds/mf_E_t_Xp.pkl") ##OY tilt DX plus#########
#manifold  = Manifold.load("./Script/Jellyfisch/PC_analysis/manifolds/mf_E_t_Yp.pkl") ##OY tilt DY plus #########
#manifold  = Manifold.load("./Script/Jellyfisch/PC_analysis/manifolds/mf_E_t_Xsp.pkl") ##OY tilt DX str plus#########
#manifold  = Manifold.load("./Script/Jellyfisch/PC_analysis/manifolds/mf_E_t_Xsn.pkl") ##OY tilt DX str neg#########
#manifold  = Manifold.load("./Script/Jellyfisch/PC_analysis/manifolds/mf_E_t_Ysn.pkl") ##OY tilt DY str neg #########
#manifold  = Manifold.load("./Script/Jellyfisch/PC_analysis/manifolds/mf_E_t.pkl") ##OY tilt  #########
#manifold  = Manifold.load("./Script/Jellyfisch/PC_analysis/E coil/manifolds/mf_E_Xp.pkl") ##DX plus  #########
#manifold  = Manifold.load("./Script/Jellyfisch/PC_analysis/manifolds/mf_E_Xsn.pkl") ##DX str neg  #########
#manifold  = Manifold.load("./Script/Jellyfisch/PC_analysis/manifolds/mf_E_Xsp.pkl") ##DX str pos  #########


#manifold  = Manifold.load("./Script/Jellyfisch/PC_analysis/manifolds/mf_E_11.pkl") ##OY tilt str pos  #########
#manifold  = Manifold.load("./Script/Jellyfisch/PC_analysis/manifolds/mf_E_12.pkl") ##OY tilt str neg  #########


#E coil isolated

#manifold1  = Manifold.load("./Script/Jellyfisch/PC_analysis/E coil/manifolds/mf_E1_to_E4.pkl") ##Only E1 to E4  #########
#manifold2 = Manifold.load("./Script/Jellyfisch/PC_analysis/E coil/manifolds/mf_E3_to_E4.pkl") ##Only E3 to E4  #########
#manifold1 = Manifold.load("./Script/Jellyfisch/PC_analysis/E coil/manifolds/mf_E5_to_E8.pkl") ##Only E5 to E8  #########
#manifold2 = Manifold.load("./Script/Jellyfisch/PC_analysis/E coil/manifolds/mf_E3.pkl") ##Only E3 #########
#manifold2 = Manifold.load("./Script/Jellyfisch/PC_analysis/E coil/manifolds/mf_E4.pkl") ##Only E4 #########

eps_s1 = 8.2e-07
eps_u1 = 7.4e-07


manifold1 = Manifold.load("./script/jellyfisch/PC_analysis/E coil/manifolds/mf_E_BO78_only_E4_tilt.pkl") ##Only E4 #########
# manifold1.find_clinics(first_guess_eps_s=eps_s1, first_guess_eps_u=eps_u1,n_points=2)
# manifold1.save('./script/jellyfisch/PC_analysis/E coil/manifolds/mf_E_BO78_E4_tilt.pkl')

#manifold1.compute_turnstile_areas()
#manifold1.save('./script/jellyfisch/PC_analysis/E coil/manifolds/mf_E_BO78_E4_tilt.pkl')


# A1=manifold1.turnstile_areas[0]
# A2=manifold1.turnstile_areas[1]

#manifold2  = Manifold.load("./script/jellyfisch/PC_analysis/E coil/manifolds/mf_E_BO78_no_E4.pkl")
# manifold2.find_clinics(first_guess_eps_s=eps_s1, first_guess_eps_u=eps_u1,n_points=2)
# manifold2.save('./script/jellyfisch/PC_analysis/E coil/manifolds/mf_E_BO78_no_E4.pkl')

# manifold2.compute_turnstile_areas()
# manifold2.save('./script/jellyfisch/PC_analysis/E coil/manifolds/mf_E_BO78_no_E4.pkl')

# B1=manifold2.turnstile_areas[0]
# B2=manifold2.turnstile_areas[1]


manifold3  = Manifold.load("./script/jellyfisch/77021/1_20/backoff/manifolds_P/OS_BO78/mf_1T_n20.pkl") ##with DX all #########
#manifold3.find_clinics(first_guess_eps_s=eps_s1, first_guess_eps_u=eps_u1,n_points=2)
# manifold3.save("./script/jellyfisch/77021/1_20/backoff/manifolds_P/OS_BO78/mf_1T_n20.pkl")
# manifold3.compute_turnstile_areas()
# manifold3.save("./script/jellyfisch/77021/1_20/backoff/manifolds_P/OS_BO78/mf_1T_n20.pkl")

# C1=manifold3.turnstile_areas[0]
# C2=manifold3.turnstile_areas[1]

#manifold1  = Manifold.load("./Script/Jellyfisch/77021/77021_120/backoff/manifolds_P/mf_1T_NA.pkl") ##with DX all #########

#manifold = Manifold.load("./script/jellyfisch/PC_analysis/E coil/manifolds/mf_E_BO78_no_E4.pkl") ##Only E4 #########


####################
####################


# manifold.plot(stepsize_limit=0.2, ax=ax, markersize=0, lw=0.7)



manifold1.plot(stepsize_limit=0.2, which="unstable", ax=ax, markersize=0, lw=0.7, colors=["green", "xkcd:red"])
manifold3.plot(stepsize_limit=0.2, which="unstable", ax=ax, markersize=0, lw=0.7, colors=["navajowhite", "xkcd:lightblue"])
#manifold3.plot(stepsize_limit=0.2, which="unstable", ax=ax, markersize=0, lw=0.7, colors=["orange", "blue"])
#ax.legend(['E4 only st.','E4 only unst.','no E4 unst.','no E4 st.','LSE st.', 'LSE unst.'], loc='lower right', fontsize=7)
ax.legend(['E4 only unst.', 'LSE unst.'], loc='lower right', fontsize=11)

ratio=2.8158



# texte = (
#     r"$\phi_{only E4}=$" + f"{A1:.4f}\n"
#     r"$\phi_{no E4}=$" + f"{B1:.4f}\n"
#     r"$\phi_{all}=$" + f"{C1:.4f}"
# )

# ax.text(0.9, -0.46, texte, fontsize=10, color='blue')

ax.set_xlim(0.7, 1.)
ax.set_ylim(-0.6, -0.1)
ax.set_xlabel(r"$R[m]$")
ax.set_ylabel(r"$Z[m]$")
ax.set_aspect('equal') 
ax.set_title('E4 isolated from LSE results')
plt.savefig('./script/jellyfisch/PC_analysis/E coil/figures/E4_isolated.png', bbox_inches='tight', pad_inches=0.0, dpi=720)
plt.show()





