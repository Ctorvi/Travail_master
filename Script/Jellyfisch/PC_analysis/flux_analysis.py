from pyoculus.fields import AxisymmetricCylindricalGridField
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
import numpy as np
from matplotlib import pyplot as plt
import logging
from scipy.io import loadmat
from matplotlib.tri import Triangulation

fig, ax = plt.subplots(1, 1, figsize=(6, 8))

eps_s1 = 8.2e-07
eps_u1 = 7.4e-07

plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.size": 10,
    }
)

logging.basicConfig(level=logging.DEBUG)

JFField1 = AxisymmetricCylindricalGridField.from_matlab_file('./script/jellyfisch/PC_analysis/E coil/pert_files/BO78_no_E4.mat', with_perturbation=True)

section = CylindricalBfieldSection(JFField1,phi0=1.9,R0=0.88, Z0=0)

top_o = FixedPoint(section)
top_o.find(1, [0.9,0.18], method='scipy.root')
top_o_coord = top_o.coords[0]
#top_o.plot(ax=ax, marker='o', color="xkcd:crimson")

x_point1= FixedPoint(section)
x_point1.find(1, [0.8,-0.27], method='scipy.root')
x_point1_coord = x_point1.coords[0]
#x_point1.plot(ax=ax, marker='x', color="xkcd:crimson")

# manifold1 = Manifold(section, x_point1, x_point1,-x_point1_coord+top_o_coord, -x_point1_coord+top_o_coord)
# manifold1.compute(eps_s=9e-6, eps_u=8e-6, nint_s=8, nint_u=14, neps_s=80, neps_u=240) #8 14 240
# manifold1.find_clinics(first_guess_eps_s=eps_s1, first_guess_eps_u=eps_u1,n_points=2)
# manifold1.compute_turnstile_areas()

#manifold1.save('./script/jellyfisch/PC_analysis/E coil/manifolds/mf_E_BO78_no_E4.pkl')
manifold1 = Manifold.load("./script/jellyfisch/PC_analysis/E coil/manifolds/mf_E_BO78_no_E4.pkl") ##Only E4 #########

A1=manifold1.turnstile_areas[0]
A2=manifold1.turnstile_areas[1]




####Second turnstile #########




JFField2 = AxisymmetricCylindricalGridField.from_matlab_file('./script/jellyfisch/PC_analysis/E coil/pert_files/BO78_only_E4_tilt.mat', with_perturbation=True)

section2 = CylindricalBfieldSection(JFField2,phi0=1.9,R0=0.88, Z0=0)

top_o_2 = FixedPoint(section2)
top_o_2.find(1, [0.9,0.18], method='scipy.root')
top_o_coord_2 = top_o_2.coords[0]
#top_o.plot(ax=ax, marker='o', color="xkcd:crimson")

x_point2= FixedPoint(section2)
x_point2.find(1, [0.8,-0.27], method='scipy.root')
x_point2_coord = x_point2.coords[0]
#x_point1.plot(ax=ax, marker='x', color="xkcd:crimson")

# manifold2 = Manifold(section2, x_point2, x_point2,-x_point2_coord+top_o_coord_2, -x_point2_coord+top_o_coord_2)
# manifold2.compute(eps_s=9e-6, eps_u=8e-6, nint_s=8, nint_u=14, neps_s=80, neps_u=240) #8 14 240
# manifold2.find_clinics(first_guess_eps_s=eps_s1, first_guess_eps_u=eps_u1,n_points=2)
# manifold2.compute_turnstile_areas()
# manifold2.save('./script/jellyfisch/PC_analysis/E coil/manifolds/mf_E_BO78_only_E4_tilt.pkl')


manifold2 = Manifold.load("./script/jellyfisch/PC_analysis/E coil/manifolds/mf_E_BO78_only_E4_tilt.pkl") ##Only E4 #########

B1=manifold2.turnstile_areas[0]
B2=manifold2.turnstile_areas[1]


### third turnstile


manifold3  = Manifold.load("./script/jellyfisch/77021/1_20/backoff/manifolds_P/OS_BO78/mf_1T_n20.pkl") ##with DX all #########
# #manifold3.find_clinics(first_guess_eps_s=eps_s1, first_guess_eps_u=eps_u1,n_points=2)
# # manifold3.save("./script/jellyfisch/77021/1_20/backoff/manifolds_P/OS_BO78/mf_1T_n20.pkl")
manifold3.compute_turnstile_areas()
# manifold3.save("./script/jellyfisch/77021/1_20/backoff/manifolds_P/OS_BO78/mf_1T_n20.pkl")

C1=manifold3.turnstile_areas[0]
C2=manifold3.turnstile_areas[1]

#manifold1  = Manifold.load("./Script/Jellyfisch/77021/77021_120/backoff/manifolds_P/mf_1T_NA.pkl") ##with DX all #########
#manifold = Manifold.load("./script/jellyfisch/PC_analysis/E coil/manifolds/mf_E_BO78_no_E4.pkl") ##Only E4 #########

####################
####################

# manifold.plot(stepsize_limit=0.2, ax=ax, markersize=0, lw=0.7)
manifold1.plot_clinics(ax=ax)
manifold2.plot_clinics(ax=ax)
#manifold3.plot_clinics(ax=ax, color="xkcd:orange")

manifold1.plot(stepsize_limit=0.2, ax=ax, markersize=0, lw=0.7, colors=["green", "xkcd:red"])
manifold2.plot(stepsize_limit=0.2, ax=ax, markersize=0, lw=0.7, colors=["navajowhite", "xkcd:cyan"])
#manifold3.plot(stepsize_limit=0.2, ax=ax, markersize=0, lw=0.7, colors=["orange", "blue"])
#ax.legend(['E4 only st.','E4 only unst.','no E4 unst.','no E4 st.','LSE st.', 'LSE unst.'], loc='lower right', fontsize=7)

texte = (
    r"$\phi_{only E4}=$" + f"{A1:.4f}\n"
    r"$\phi_{no E4}=$" + f"{B1:.4f}\n"
    r"$\phi_{all}=$" + f"{C1:.4f}"
)

ax.text(0.9, -0.46, texte, fontsize=10, color='blue')

ax.set_xlim(0.7, 1.)
ax.set_ylim(-0.6, -0.1)
ax.set_xlabel(r"$R[m]$")
ax.set_ylabel(r"$Z[m]$")
ax.set_aspect('equal') 
ax.set_title('E4 isolated from LSE results')
#plt.savefig('./script/jellyfisch/PC_analysis/E coil/figures/E4_vs_noE4_turnstile_flux.png', bbox_inches='tight', dpi=720)
plt.show()





