from scipy.ndimage import gaussian_filter1d

from scipy.signal import argrelextrema, find_peaks
from pyoculus.fields import AxisymmetricCylindricalGridField
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
import numpy as np
from matplotlib import pyplot as plt
import logging
from scipy.io import loadmat
from matplotlib.tri import Triangulation


def line_curve_intersections(p0, p1, curves, tol=1e-9, ax=None, plot_kwargs=None):
    """
    Trouve les intersections entre la droite (p0->p1) et UNE courbe (polyligne).
    Retour: (pts Nx2, t_params, u_params) -- tous numpy arrays (vide si aucune intersection).
    """
    P = np.asarray(p0, dtype=float)
    Q = np.asarray(p1, dtype=float)
    s = Q - P

    def cross(u, v):
        return u[0]*v[1] - u[1]*v[0]

    if plot_kwargs is None:
        plot_kwargs = dict(color='k', marker='x', markersize=6, mew=1.5, linestyle='none')

    # accepter soit directement une array, soit une liste/tuple contenant une seule curve
    if isinstance(curves, (list, tuple)):
        if len(curves) == 0:
            return np.empty((0,2)), np.empty((0,)), np.empty((0,))
        curve = curves[0]
    else:
        curve = curves

    pts_arr = np.asarray(curve, dtype=float)
    if pts_arr.ndim != 2 or pts_arr.shape[0] < 2:
        return np.empty((0,2)), np.empty((0,)), np.empty((0,))

    hits = []
    t_list = []
    u_list = []
    for k in range(len(pts_arr) - 1):
        A = pts_arr[k]
        B = pts_arr[k + 1]
        r = B - A
        denom = cross(r, s)
        if abs(denom) < 1e-15:
            continue
        PA = P - A
        t = cross(PA, s) / denom
        u = cross(PA, r) / denom
        if -tol <= t <= 1 + tol and -tol <= u <= 1 + tol:
            pt = A + t * r
            hits.append(pt)
            t_list.append(t)
            u_list.append(u)
            if ax is not None:
                ax.plot(pt[0], pt[1], **plot_kwargs)
                ax.plot([P[0], Q[0]], [P[1], Q[1]], color='black', linestyle='--', linewidth=0.8, alpha=0.5)
                ax.scatter(P[0], P[1], color='black', s=10,marker='s')
                ax.scatter(Q[0], Q[1], color='black', s=10,marker='s')

    if len(hits) == 0:
        return np.empty((0,2)), np.empty((0,)), np.empty((0,))

    return np.vstack(hits), np.asarray(t_list, dtype=float), np.asarray(u_list, dtype=float)
# ...existing code...

def line_curve_MF_intersections(p0, p1, manifolds, unstable=True, tol=1e-9, ax1=None, ax2=None, plot_kwargs1=None, plot_kwargs2=None):
    """
    Wrapper qui accepte soit un objet manifold (avec .unstable/.stable) soit directement
    une polyligne (Nx2). Retourne (pts Nx2, dist_along_line).

    pts: intersection point coordinates
    dists: distance from p0 along the line 
    ax: to plot pts if needed
    """
    P = np.asarray(p0, dtype=float)
    Q = np.asarray(p1, dtype=float)
    total_dist = np.linalg.norm(Q - P)
    if total_dist == 0:
        total_dist = 1.0  # éviter division par zéro si p0==p1

    # extraire la courbe (polyligne) depuis l'objet manifold si nécessaire
    if hasattr(manifolds, 'unstable') or hasattr(manifolds, 'stable'):
        C = manifolds.unstable if unstable else manifolds.stable
        curve = np.asarray(C, dtype=float)
    else:
        curve = np.asarray(manifolds, dtype=float)

    pts, t_params, u_params = line_curve_intersections(p0, p1, curve, tol=tol, ax=ax1, plot_kwargs=plot_kwargs1)
    if pts.size == 0:
        return np.empty((0,2)), np.empty((0,))
    
    dists = np.linalg.norm(pts - P, axis=1)
    middle_zone_MF=[]
    if ax2 is not None:
             for i, (a, b) in enumerate(zip(dists[::2], dists[1::2])):
                    lbl = 'Leg zone (MF intersection)' if i == 0 else None
                    lbl2= 'Half width of leg zone' if i == 0 else None

                    ax2.axvspan(a, b, color='lightgrey', alpha=0.3, label=lbl)
                    ax2.axvline(0.5*(a + b), color='red', linestyle='--', linewidth=0.8, alpha=0.8, label=lbl2)
                    middle_zone_MF.append(0.5*(a + b))
                     
           # else:
           #     ax2.axvline(x=x, **plot_kwargs2)
    periodicity = np.mean(np.diff(middle_zone_MF))
    std_periodicity = np.std(np.diff(middle_zone_MF))
    return pts, dists, periodicity, std_periodicity
# ...existing code...


def sample_emi_along_line(p0, p1, tri, emi, n_samples=300, limit_find_peaks=0.2,sig=5.0,ax=None, plot_kwargs=None):# plot_profile=True):
    """
    Échantillonne l'émissivité `emi` le long de la droite p0->p1.
    - p0, p1 : (R,Z)
    - tri : matplotlib.tri.Triangulation
    - emi : valeurs par noeud (len == n_points) ou par triangle (len == n_tri)
    Retour : 
    - pts (n,2): coordinates of points along the line 
    - emi_vals (n,): emissivity values at these points (NaN if outside triangulation)
    - dist (n,): distance from p0 along the line 
    """

    p0 = np.asarray(p0, dtype=float)
    p1 = np.asarray(p1, dtype=float)
    t = np.linspace(0.0, 1.0, int(n_samples))
    pts = p0[None, :] + np.outer(t, (p1 - p0))          # (n,2)
    x = pts[:, 0]; y = pts[:, 1]
    dist = np.hypot(x - x[0], y - y[0])

    ntri = tri.triangles.shape[0]
        # interpolation linéaire depuis les noeuds
    trif = tri.get_trifinder()
    tri_idx = trif(x, y)           # -1 si hors triangulation
    emi_vals = np.full_like(x, np.nan, dtype=float)
    mask = tri_idx >= 0
    if np.any(mask):
         emi_vals[mask] = emi[tri_idx[mask]]

    
    y_gauss = gaussian_filter1d(emi_vals, sigma=sig)
         
    if ax is not None:
        # safe use of plot_kwargs
        if isinstance(plot_kwargs, dict):
            ax.plot(dist, emi_vals, label='Emissivity along line', **plot_kwargs)
            ax.plot(dist, y_gauss, color='green', linewidth=1)#label='Gaussian smoothed',

        else:
            ax.plot(dist, emi_vals,color='royalblue', label='Emissivity along line')
            ax.plot(dist, y_gauss, color='royalblue', alpha=0.5, linewidth=1)#label='Gaussian smoothed',

    idx_end=np.argmin(np.abs(dist - limit_find_peaks))
    emi_vals_end=emi_vals[:idx_end]
    y_gauss_end=y_gauss[:idx_end]

    max_idx, props = find_peaks(y_gauss_end, distance=10, prominence=1e-3)  # ajuster distance/prominence
    max_dist = dist[max_idx]

    ax.scatter(max_dist, y_gauss_end[max_idx], color='blue', zorder=10)

    for x in np.atleast_1d(max_dist):
        lbl = 'peak emissivity' if x == np.atleast_1d(max_dist)[0] else None
        ax.axvline(float(x), color='blue', linestyle='--', linewidth=0.8, alpha=0.8,label=lbl)


    periodicity = np.mean(np.diff(max_dist))
    std_periodicity = np.std(np.diff(max_dist))


    
    return pts, emi_vals, dist, periodicity,std_periodicity




   


    