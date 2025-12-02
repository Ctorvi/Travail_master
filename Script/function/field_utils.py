import numpy as np
from scipy.special import ellipk, ellipe
from scipy.spatial import cKDTree

mu0 = 4e-7 * np.pi
def mask_within_distance_to_curve(R_grid, Z_grid, curve_points, r_cut):
    """
    Retourne un masque booléen de la même forme que R_grid/Z_grid,
    True si la cellule est à distance <= r_cut de la courbe (curve_points Nx2).
    """
    # forme (N,2)
    pts = np.vstack((R_grid.ravel(), Z_grid.ravel())).T
    tree = cKDTree(curve_points)
    dists, _ = tree.query(pts, distance_upper_bound=r_cut)
    mask = (dists <= r_cut)
    return mask.reshape(R_grid.shape)


def zero_near_curve(B_R, B_Z, psi, R_grid, Z_grid, curve_points, r_cut=0.05):
    """
    Met B_R, B_Z et psi à zéro pour tous les points à distance <= r_cut de la courbe.
    Modifie les arrays in-place et les renvoie.
    """
    mask = mask_within_distance_to_curve(R_grid, Z_grid, curve_points, r_cut)
    B_R[mask] = 0.0
    B_Z[mask] = 0.0
    if psi is not None:
        psi[mask] = 0.0
    return B_R, B_Z, psi



def B_loop_cyl(R, Z, a, Z0, I):
    R = np.asarray(R, dtype=float)
    Z = np.asarray(Z, dtype=float)
    z = Z - Z0
    r = np.maximum(R, 1e-12)
    Delta = (a + r)**2 + z**2
    k2 = 4 * a * r / Delta
    k2 = np.clip(k2, 0.0, 1.0 - 1e-12)
    K = ellipk(k2)
    E = ellipe(k2)
    denom = np.sqrt(Delta)
    D = (a - r)**2 + z**2
    Br = mu0 * I * z / (2 * np.pi * r * denom) * (
         -K + (a**2 + r**2 + z**2) / D * E
    )
    Bz = mu0 * I / (2 * np.pi * denom) * (
          K + (a**2 - r**2 - z**2) / D * E
    )
    Bphi = np.zeros_like(Br)
    return Br, Bz, Bphi

def B_sum_of_loops(Rg, Zg, R_loops, Z_loops, I_loops):
    B_R_tot = np.zeros_like(Rg)
    B_Z_tot = np.zeros_like(Rg)
    B_phi_tot = np.zeros_like(Rg)
    for i,j in np.ndindex(Rg.shape):
        for k in range(len(R_loops)):
            Br, Bz, Bphi = B_loop_cyl(Rg[i,j], Zg[i,j], R_loops[k], Z_loops[k], I_loops)
            B_R_tot[i,j] += Br
            B_Z_tot[i,j] += Bz
            B_phi_tot[i,j] += Bphi
    return B_R_tot, B_Z_tot, B_phi_tot

def BR_BZ_PSI(R_grid, Z_grid, B_R_tot, B_Z_tot, x_p, r_cut=0.05):
    psi = np.zeros_like(R_grid)
    psi[0,0] = 0.0
    NR = R_grid.shape[1]
    NZ = R_grid.shape[0]
    for i in range(1, NR):
        R_mid = 0.5*(R_grid[0,i] + R_grid[0,i-1])
        BZ_mid = 0.5*(B_Z_tot[0,i] + B_Z_tot[0,i-1])
        dR = R_grid[0,i] - R_grid[0,i-1]
        psi[0,i] = psi[0,i-1] + R_mid * BZ_mid * dR
    for i in range(NR):
        for j in range(1, NZ):
            R_here = R_grid[j,i]
            BR_mid = 0.5*(B_R_tot[j,i] + B_R_tot[j-1,i])
            dZ = Z_grid[j,i] - Z_grid[j-1,i]
            psi[j,i] = psi[j-1,i] - R_here * BR_mid * dZ
    psi = 2*np.pi*psi

    Rc, Zc = x_p[0], x_p[1]
    dist = np.hypot(R_grid - Rc, Z_grid - Zc)
    mask = dist <= r_cut
    B_R_tot[mask] = 0.0
    B_Z_tot[mask] = 0.0
    psi[mask] = 0.0



    return B_R_tot, B_Z_tot, psi







# psi = np.zeros_like(R_grid)

# # 1. Intégration le long de Z = Z0 (ligne de référence, par ex j=0)
# psi[0,0] = 0.0  # condition de référence

# NR = R_grid.shape[1]
# NZ = R_grid.shape[0]

# for i in range(1, NR):
#     R_mid = 0.5*(R_grid[0,i] + R_grid[0,i-1])
#     BZ_mid = 0.5*(B_Z_tot[0,i] + B_Z_tot[0,i-1])
#     dR = R_grid[0,i] - R_grid[0,i-1]
#     psi[0,i] = psi[0,i-1] + R_mid * BZ_mid * dR  # ψ(R_i,Z0)

# # 2. Intégration en Z pour chaque colonne R_i
# for i in range(NR):
#     for j in range(1, NZ):
#         R_here = R_grid[j,i]
#         BR_mid = 0.5*(B_R_tot[j,i] + B_R_tot[j-1,i])
#         dZ = Z_grid[j,i] - Z_grid[j-1,i]
#         psi[j,i] = psi[j-1,i] - R_here * BR_mid * dZ



# psi=2*np.pi*psi  # en Weber

# ## zero les éléments DANS un disque de rayon r_cut autour de (Rc, Zc)
# Rc, Zc = x_p[0],x_p[1]   # centre du disque (à adapter)
# r_cut = 0.05          # rayon (à adapter)
# dist = np.hypot(R_grid - Rc, Z_grid - Zc)  # distance euclidienne sur la grille
# mask = dist <= r_cut   # True pour les points à l'intérieur du disque

# B_R_tot[mask] = 0.0
# B_Z_tot[mask] = 0.0
# B_phi_tot[mask] = 0.0