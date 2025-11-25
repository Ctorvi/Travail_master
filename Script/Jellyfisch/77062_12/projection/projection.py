import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from pyoculus.fields import AxisymmetricCylindricalGridField
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
import logging
from scipy.io import loadmat
from matplotlib.tri import Triangulation



fig, ax = plt.subplots(figsize=(6,6))
# -----------------------------
# 1) Camera parameters
# -----------------------------
C = np.array([1.127, -0.287, -0.016])  # camera position (m)
view_dir = np.array([-0.914, -0.401, -0.059])  # view direction

def camera_rotation_from_direction(d):
    d = d / np.linalg.norm(d)
    up = np.array([0,0,1])
    if abs(np.dot(d, up)) > 0.99:
        up = np.array([0,1,0])
    x_cam = np.cross(up, d)
    x_cam /= np.linalg.norm(x_cam)
    y_cam = np.cross(d, x_cam)
    y_cam /= np.linalg.norm(y_cam)
    R = np.vstack([x_cam, y_cam, d])
    return R

Rcam = camera_rotation_from_direction(view_dir)



def toroidal_slice_to_xyz(R, Z, phi):
    X = R * np.cos(phi)
    Y = R * np.sin(phi)
    return np.vstack([X.ravel(), Y.ravel(), Z.ravel()]).T

def project_points_camera_frame(P, Rcam, C):
    Pc = (Rcam @ (P - C).T).T
    valid = Pc[:,2] > 0
    u = Pc[:,0] / Pc[:,2]
    v = Pc[:,1] / Pc[:,2]
    return u, v, valid




manifold1=Manifold.load('./Script/Jellyfisch/77021/77021_120/backoff/manifolds_P/mf_1T_NA.pkl')

#manifold1.plot(ax=ax,stepsize_limit=0.1)
# manifold1.plot_clinics(ax=ax)
# manifold1.plot_filled_lobe(ax=ax, lobe_number=15,alpha=0.3)
# manifold1.plot_filled_lobe(ax=ax, lobe_number=14,alpha=0.3)
# manifold1.plot_filled_lobe(ax=ax, lobe_number=13,alpha=0.3)
# manifold1.plot_filled_lobe(ax=ax, lobe_number=12,alpha=0.3)
# manifold1.plot_filled_lobe(ax=ax, lobe_number=11,alpha=0.3)

#LB_s=[]
# LB1_s,LB1_u=manifold1.get_lobe_boundary(lobe_number=11)
# LB2_s,LB2_u=manifold1.get_lobe_boundary(lobe_number=12)
# LB3_s,LB3_u=manifold1.get_lobe_boundary(lobe_number=13)
# LB4_s,LB4_u=manifold1.get_lobe_boundary(lobe_number=14)
# LB5_s,LB5_u=manifold1.get_lobe_boundary(lobe_number=15)


# LB1=np.concatenate([LB1_u,LB1_s])
# LB2=np.concatenate([LB2_u,LB2_s])
# LB3=np.concatenate([LB3_u,LB3_s])
# LB4=np.concatenate([LB4_u,LB4_s])
# LB5=np.concatenate([LB5_u,LB5_s])


data = np.load("./Script/Jellyfisch/77062_12/projection/Lobe_boundary.npz")
LB1 = data["LB1"]
LB2 = data["LB2"]
LB3 = data["LB3"]
LB4 = data["LB4"]
LB5 = data["LB5"]


DB = {
    "LB1": LB1,
    "LB2": LB2,
    "LB3": LB3,
    "LB4": LB4,
    "LB5": LB5
}

LB = ["LB1", "LB2", "LB3", "LB4", "LB5"]


#np.savez("./Script/Jellyfisch/77062_12/projection/Lobe_boundary.npz", LB1=LB1, LB2=LB2, LB3=LB3, LB4=LB4, LB5=LB5)


scat = ax.scatter([], [], s=0.7, color='black')
ax.set_xlabel("u [m] (camera plane)")
ax.set_ylabel("v [m] (camera plane)")
ax.set_title("Projection of manifold/legs from MANTIS pov (285-295 deg)")


ax.set_aspect('equal')



# -----------------------------
# 3) Animation setup
# -----------------------------


FOVy_deg = 69.1  # horizontal
FOVx_deg = 52.4  # vertical

FOVx = np.deg2rad(FOVx_deg)
FOVy = np.deg2rad(FOVy_deg)

# xlim = (-np.tan(FOVx/2), np.tan(FOVx/2))
# ylim = (-np.tan(FOVy/2), np.tan(FOVy/2))

xlim = (-0.2, 0.2)
ylim = (-0.6, 0.)

ax.set_xlim(xlim)
ax.set_ylim(ylim)




# 4) Animation function
# -----------------------------


theta_int = np.linspace(285, 295, 150) 

for i in range(len(LB)):
  LBi=DB[LB[i]]

  angle= np.deg2rad(292.5)
  P_plane = toroidal_slice_to_xyz(LBi[:,0], LBi[:,1], angle)
  u, v, valid = project_points_camera_frame(P_plane, Rcam, C)
  u_valid = u[valid]
  v_valid = v[valid]
  ax.plot(u_valid, v_valid,'-',linewidth=0.7,color='black')


  for theta in theta_int:
        phi_plane = np.deg2rad(theta)  # 0->360 degrees
        P_plane = toroidal_slice_to_xyz(LBi[:,0], LBi[:,1], phi_plane)
        u, v, valid = project_points_camera_frame(P_plane, Rcam, C)
        u_valid = u[valid]
        v_valid = v[valid]
        #ax.scatter(u_valid, v_valid, s=0.5)
        ax.fill(u_valid, v_valid, alpha=0.05) # color='black')



plt.savefig("./Script/Jellyfisch/77062_12/projection/Lobe_boundary_projection.png", dpi=300)
plt.show()