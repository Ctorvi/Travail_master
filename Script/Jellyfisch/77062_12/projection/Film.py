import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from pyoculus.fields import AxisymmetricCylindricalGridField
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
import logging
from scipy.io import loadmat
from matplotlib.tri import Triangulation
from mpl_toolkits.mplot3d import Axes3D  



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
manifold2=Manifold.load('./Script/Jellyfisch/77021/77021_120/backoff/manifolds_P/mf_1B_NA.pkl')
manifold3=Manifold.load('./Script/Jellyfisch/77021/77021_120/backoff/manifolds_P/mf_2T_NA.pkl')
manifold4=Manifold.load('./Script/Jellyfisch/77021/77021_120/backoff/manifolds_P/mf_2B_NA.pkl')



manifold1_points_u = manifold1.unstable
manifold1_points_s = manifold1.stable

manifold2_points_u = manifold2.unstable
manifold2_points_s = manifold2.stable

manifold3_points_u = manifold3.unstable
manifold3_points_s = manifold3.stable

manifold4_points_u = manifold4.unstable
manifold4_points_s = manifold4.stable


Rgrid1_u = manifold1_points_u[:,0].reshape(-1, 1)
Zgrid1_u = manifold1_points_u[:,1].reshape(-1, 1)

Rgrid1_s = manifold1_points_s[:,0].reshape(-1, 1)
Zgrid1_s = manifold1_points_s[:,1].reshape(-1, 1)

Rgrid2_u = manifold2_points_u[:,0].reshape(-1, 1)
Zgrid2_u = manifold2_points_u[:,1].reshape(-1, 1)

Rgrid2_s = manifold2_points_s[:,0].reshape(-1, 1)
Zgrid2_s = manifold2_points_s[:,1].reshape(-1, 1)

Rgrid3_u = manifold3_points_u[:,0].reshape(-1, 1)
Zgrid3_u = manifold3_points_u[:,1].reshape(-1, 1)

Rgrid3_s = manifold3_points_s[:,0].reshape(-1, 1)
Zgrid3_s = manifold3_points_s[:,1].reshape(-1, 1)   

Rgrid4_u = manifold4_points_u[:,0].reshape(-1, 1)
Zgrid4_u = manifold4_points_u[:,1].reshape(-1, 1)

Rgrid4_s = manifold4_points_s[:,0].reshape(-1, 1)
Zgrid4_s = manifold4_points_s[:,1].reshape(-1, 1)



#Rgrid, Zgrid = np.meshgrid(Rgrid, Zgrid)



scat = ax.scatter([], [], s=0.7, color='red')
scat1 = ax.scatter([], [], s=0.7, color='blue')

scat2 = ax.scatter([], [], s=0.7, color='red')
scat3 = ax.scatter([], [], s=0.7, color='blue')

scat4 = ax.scatter([], [], s=0.7, color='red')
scat5 = ax.scatter([], [], s=0.7, color='blue')

scat6 = ax.scatter([], [], s=0.7, color='red')
scat7 = ax.scatter([], [], s=0.7, color='blue')

ax.set_xlabel("u [m] (camera plane)")
ax.set_ylabel("v [m] (camera plane)")
ax.set_title("Poloidal plane moving around torus")


ax.set_aspect('equal')




# -----------------------------
# 3) Animation setup
# -----------------------------


FOVy_deg = 69.1  # horizontal
FOVx_deg = 52.4  # vertical

FOVx = np.deg2rad(FOVx_deg)
FOVy = np.deg2rad(FOVy_deg)

xlim = (-np.tan(FOVx/2), np.tan(FOVx/2))
ylim = (-np.tan(FOVy/2), np.tan(FOVy/2))

ax.set_xlim(xlim)
ax.set_ylim(ylim)


print("xlim =", xlim)
print("ylim =", ylim)


# 4) Animation function
# -----------------------------

def update(frame):
    phi_plane = np.deg2rad(frame)  # 0->360 degrees

    P_plane1_u = toroidal_slice_to_xyz(Rgrid1_u, Zgrid1_u, phi_plane)
    P_plane1_s = toroidal_slice_to_xyz(Rgrid1_s, Zgrid1_s, phi_plane)

    P_plane2_u = toroidal_slice_to_xyz(Rgrid2_u, Zgrid2_u, phi_plane)
    P_plane2_s = toroidal_slice_to_xyz(Rgrid2_s, Zgrid2_s, phi_plane)    

    P_plane3_u = toroidal_slice_to_xyz(Rgrid3_u, Zgrid3_u, phi_plane)
    P_plane3_s = toroidal_slice_to_xyz(Rgrid3_s, Zgrid3_s, phi_plane)

    P_plane4_u = toroidal_slice_to_xyz(Rgrid4_u, Zgrid4_u, phi_plane)
    P_plane4_s = toroidal_slice_to_xyz(Rgrid4_s, Zgrid4_s, phi_plane)

    u1_u, v1_u, valid1_u = project_points_camera_frame(P_plane1_u, Rcam, C)
    u1_s, v1_s, valid1_s = project_points_camera_frame(P_plane1_s, Rcam, C)
    u2_u, v2_u, valid2_u = project_points_camera_frame(P_plane2_u, Rcam, C)
    u2_s, v2_s, valid2_s = project_points_camera_frame(P_plane2_s, Rcam, C)
    u3_u, v3_u, valid3_u = project_points_camera_frame(P_plane3_u, Rcam, C)
    u3_s, v3_s, valid3_s = project_points_camera_frame(P_plane3_s, Rcam, C)
    u4_u, v4_u, valid4_u = project_points_camera_frame(P_plane4_u, Rcam, C)
    u4_s, v4_s, valid4_s = project_points_camera_frame(P_plane4_s, Rcam, C)



    u_valid1_u = u1_u[valid1_u]
    v_valid1_u = v1_u[valid1_u]
    u_valid1_s = u1_s[valid1_s]
    v_valid1_s = v1_s[valid1_s]

    u_valid2_u = u2_u[valid2_u]
    v_valid2_u = v2_u[valid2_u]
    u_valid2_s = u2_s[valid2_s]
    v_valid2_s = v2_s[valid2_s]

    u_valid3_u = u3_u[valid3_u]
    v_valid3_u = v3_u[valid3_u]
    u_valid3_s = u3_s[valid3_s]
    v_valid3_s = v3_s[valid3_s]

    u_valid4_u = u4_u[valid4_u]
    v_valid4_u = v4_u[valid4_u]
    u_valid4_s = u4_s[valid4_s]
    v_valid4_s = v4_s[valid4_s]



    scat.set_offsets(np.c_[u_valid1_u, v_valid1_u])
    scat1.set_offsets(np.c_[u_valid1_s, v_valid1_s])

    scat2.set_offsets(np.c_[u_valid2_u, v_valid2_u])
    scat3.set_offsets(np.c_[u_valid2_s, v_valid2_s])

    scat4.set_offsets(np.c_[u_valid3_u, v_valid3_u])
    scat5.set_offsets(np.c_[u_valid3_s, v_valid3_s])

    scat6.set_offsets(np.c_[u_valid4_u, v_valid4_u])
    scat7.set_offsets(np.c_[u_valid4_s, v_valid4_s])

    ax.set_title(f"Poloidal plane φ={frame}°")

    return scat,scat1, scat2, scat3, scat4, scat5, scat6, scat7

# -----------------------------
# 5) Run animation
# -----------------------------
anim = FuncAnimation(fig, update, frames=np.arange(0, 360, 2), interval=1)


#ax.plot()

plt.show()
