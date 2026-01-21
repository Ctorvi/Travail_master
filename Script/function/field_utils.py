import numpy as np
from scipy.special import ellipk, ellipe
from scipy.spatial import cKDTree
from scipy.io import loadmat
from matplotlib.lines import Line2D
from matplotlib.tri import Triangulation
from pyoculus.utils.plot import create_canvas
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
from matplotlib.patches import Polygon as MplPolygon






def fill_between_polygons(ax, R_outer, Z_outer, R_inner, Z_inner, **kwargs):
    R1 = np.asarray(R_outer).ravel()
    Z1 = np.asarray(Z_outer).ravel()
    R2 = np.asarray(R_inner).ravel()
    Z2 = np.asarray(Z_inner).ravel()

    # retirer NaN éventuels
    mask1 = ~np.isnan(R1) & ~np.isnan(Z1)
    mask2 = ~np.isnan(R2) & ~np.isnan(Z2)
    R = np.concatenate([R1[mask1], R2[mask2][::-1]])
    Z = np.concatenate([Z1[mask1], Z2[mask2][::-1]])

    ax.fill(R, Z, picker=False, label='_nolegend_', **kwargs)


def fill_between_polygons_shapely(ax, R_outer, Z_outer, R_inner, Z_inner, color=(0.5,0.5,0.5)):

    outer = Polygon(np.column_stack([R_outer, Z_outer]))
    inner = Polygon(np.column_stack([R_inner, Z_inner]))
    area = outer.difference(inner)  # Polygon or MultiPolygon

    if area.is_empty:
        return

    # background pour "creuser" les trous si besoin
    bg = ax.get_facecolor()
    polys = area.geoms if hasattr(area, "geoms") else [area]

    for poly in polys:
        ext = np.array(poly.exterior.coords)
        ax.add_patch(MplPolygon(ext, closed=True, facecolor=color, edgecolor="none",picker=False, label='_nolegend_'))
        # dessiner les intérieurs (trous) avec la couleur de fond pour creuser
        for interior in poly.interiors:
            hole = np.array(interior.coords)
            ax.add_patch(MplPolygon(hole, closed=True, facecolor=bg, edgecolor="none",picker=False, label='_nolegend_'))



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



# ...existing code...
def zero_above_curve(B_R, B_Z, psi, R_grid, Z_grid, curve_points, mode='interp', fill_outside=False):
    """
    Met à zéro B_R, B_Z et psi pour tous les points de la grille situés *au-dessus*
    de la courbe `curve_points` (Nx2 array of [R,Z]).
    Args:
      B_R, B_Z, psi : arrays shape like R_grid/Z_grid (psi can be None)
      R_grid, Z_grid : meshgrid-like arrays of coordinates
      curve_points : array-like (N,2) of [R,Z] points describing la courbe
      mode : 'interp'|'nearest'
        - 'interp' : on interpole Z_curve(R) et on considère "au‑dessus" si Z_grid > Z_curve(R)
        - 'nearest' : on prend le Z du point de la courbe le plus proche de chaque cellule
      fill_outside : bool, si True on considère aussi points hors domaine R de la courbe
    Returns:
      B_R, B_Z, psi  (modifiés in-place puis renvoyés)
    """
    import numpy as np
    from scipy.spatial import cKDTree
    from scipy.interpolate import interp1d

    curve = np.asarray(curve_points, dtype=float)
    if curve.size == 0:
        return B_R, B_Z, psi

    R_flat = R_grid.ravel()
    Z_flat = Z_grid.ravel()

    if mode == 'interp':
        # construire une interpolation Z(R) : nécessite R strictement croissant -> trier et unique
        Rc = curve[:, 0]
        Zc = curve[:, 1]
        # trier par R
        order = np.argsort(Rc)
        Rc_s = Rc[order]
        Zc_s = Zc[order]
        # enlever doublons en R (garder la première occurrence)
        uniq_idx = np.concatenate(([0], np.where(np.diff(Rc_s) > 0)[0] + 1))
        Rc_u = Rc_s[uniq_idx]
        Zc_u = Zc_s[uniq_idx]

        if Rc_u.size < 2:
            # fallback au nearest si interpolation impossible
            mode = 'nearest'
        else:
            interp = interp1d(Rc_u, Zc_u, bounds_error=False, fill_value=np.nan)
            Z_curve_at_R = interp(R_flat)  # NaN pour R hors intervalle
            if not fill_outside:
                mask_outside = np.isnan(Z_curve_at_R)
                Z_curve_at_R[mask_outside] = -np.inf  # outside -> never > => not flagged
            else:
                # remplir par extrapolation via nearest end value
                left = Zc_u[0]; right = Zc_u[-1]
                Z_curve_at_R = np.where(np.isnan(Z_curve_at_R),
                                        np.where(R_flat < Rc_u[0], left, right),
                                        Z_curve_at_R)
            mask = Z_flat > Z_curve_at_R

    if mode == 'nearest':
        tree = cKDTree(curve)
        _, idx = tree.query(np.vstack((R_flat, Z_flat)).T)
        Z_curve_near = curve[idx, 1]
        mask = Z_flat > Z_curve_near

    # reshape mask
    mask2 = mask.reshape(R_grid.shape)

    # apply zeroing
    B_R[mask2] = 0.0
    B_Z[mask2] = 0.0
    if psi is not None:
        psi[mask2] = 0.0

    return B_R, B_Z, psi




# ...existing code...

# def zero_above_curve(B_R, B_Z, psi, R_grid, Z_grid, curve_points_to_cut):

#     Z_max=
    
#     B_R[mask] = 0.0
#     B_Z[mask] = 0.0
#     if psi is not None:
#         psi[mask] = 0.0

#     return B_R, B_Z, psi

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

# def B_sum_of_loops(Rg, Zg, R_loops, Z_loops, I_loops):
#     B_R_tot = np.zeros_like(Rg)
#     B_Z_tot = np.zeros_like(Rg)
#     B_phi_tot = np.zeros_like(Rg)
#     for i,j in np.ndindex(Rg.shape):
#         for k in range(len(R_loops)):
#             Br, Bz, Bphi = B_loop_cyl(Rg[i,j], Zg[i,j], R_loops[k], Z_loops[k], I_loops)
#             B_R_tot[i,j] += Br
#             B_Z_tot[i,j] += Bz
#             B_phi_tot[i,j] += Bphi
#     return B_R_tot, B_Z_tot, B_phi_tot


# def B_sum_of_loops(Rg, Zg, R_loops, Z_loops, I_loops):
#     """
#     Sum field of multiple circular loops.
#     Accepts Rg/Zg either as 2D meshgrid arrays or as 1D vectors (Rg = R axis, Zg = Z axis).
#     Returns arrays shaped like the evaluated grid.
#     """
#     Rg_a = np.asarray(Rg)
#     Zg_a = np.asarray(Zg)

#     # build R_grid, Z_grid depending on input shape
#     if Rg_a.ndim == 1 and Zg_a.ndim == 1:
#         R_grid, Z_grid = np.meshgrid(Rg_a, Zg_a, indexing='xy')
#     elif Rg_a.ndim == 2 and Zg_a.ndim == 2:
#         if Rg_a.shape != Zg_a.shape:
#             raise ValueError("Rg and Zg 2D arrays must have the same shape")
#         R_grid, Z_grid = Rg_a, Zg_a
#     else:
#         # try to broadcast (e.g. one is column vector, other row vector)
#         try:
#             R_grid, Z_grid = np.broadcast_arrays(Rg_a, Zg_a)
#             R_grid = R_grid.copy()
#             Z_grid = Z_grid.copy()
#         except Exception:
#             raise ValueError("Rg and Zg must be both 1D (axes) or both 2D (meshgrid) or broadcastable to the same shape")

#     B_R_tot = np.zeros_like(R_grid, dtype=float)
#     B_Z_tot = np.zeros_like(R_grid, dtype=float)
#     B_phi_tot = np.zeros_like(R_grid, dtype=float)

#     R_loops = np.asarray(R_loops).ravel()
#     Z_loops = np.asarray(Z_loops).ravel()
#     I_loops = np.asarray(I_loops)

#     n_loops = len(R_loops)
#     if len(Z_loops) != n_loops:
#         raise ValueError("R_loops and Z_loops must have the same length")
#     # allow scalar I_loops or array of same length as loops
#     if I_loops.size not in (1, n_loops):
#         raise ValueError("I_loops must be scalar or same length as R_loops")

#     for k in range(n_loops):
#         I_k = float(I_loops[k]) if I_loops.size == n_loops else float(I_loops)
#         Br, Bz, Bphi = B_loop_cyl(R_grid, Z_grid, R_loops[k], Z_loops[k], I_k)
#         B_R_tot += Br
#         B_Z_tot += Bz
#         B_phi_tot += Bphi
#     return B_R_tot, B_Z_tot, B_phi_tot
# ...existing code...
def B_sum_of_loops(Rg, Zg, R_loops, Z_loops, I_loops):
    """
    Sum field of multiple circular loops.
    Accepts:
      - Rg, Zg as 2D meshgrid arrays -> returns arrays shaped like the grid
      - Rg, Zg as 1D axis vectors -> returns full meshgrid result (len(Zg), len(Rg))
      - Rg, Zg as 1D coordinate arrays of same length -> interpreted as pointwise pairs
        and returns 1D arrays of length N (pointwise evaluation).
    """
    import numpy as np

    Rg_a = np.asarray(Rg)
    Zg_a = np.asarray(Zg)

    # detect pointwise pairs: both 1D and same length -> evaluate at N points
    if Rg_a.ndim == 1 and Zg_a.ndim == 1 and Rg_a.size == Zg_a.size:
        pointwise = True
        R_pts = Rg_a.ravel()
        Z_pts = Zg_a.ravel()
    else:
        pointwise = False

    # build R_grid, Z_grid for grid evaluation
    if not pointwise:
        if Rg_a.ndim == 1 and Zg_a.ndim == 1:
            R_grid, Z_grid = np.meshgrid(Rg_a, Zg_a, indexing='xy')
        elif Rg_a.ndim == 2 and Zg_a.ndim == 2:
            if Rg_a.shape != Zg_a.shape:
                raise ValueError("Rg and Zg 2D arrays must have the same shape")
            R_grid, Z_grid = Rg_a, Zg_a
        else:
            try:
                R_grid, Z_grid = np.broadcast_arrays(Rg_a, Zg_a)
                R_grid = R_grid.copy()
                Z_grid = Z_grid.copy()
            except Exception:
                raise ValueError("Rg and Zg must be both 1D (axes) or both 2D (meshgrid) or broadcastable to the same shape")

    # prepare outputs
    if pointwise:
        out_shape = (R_pts.size,)
        B_R_tot = np.zeros(out_shape, dtype=float)
        B_Z_tot = np.zeros(out_shape, dtype=float)
        B_phi_tot = np.zeros(out_shape, dtype=float)
    else:
        out_shape = R_grid.shape
        B_R_tot = np.zeros(out_shape, dtype=float)
        B_Z_tot = np.zeros(out_shape, dtype=float)
        B_phi_tot = np.zeros(out_shape, dtype=float)

    R_loops = np.asarray(R_loops).ravel()
    Z_loops = np.asarray(Z_loops).ravel()
    I_loops = np.asarray(I_loops)

    n_loops = len(R_loops)
    if len(Z_loops) != n_loops:
        raise ValueError("R_loops and Z_loops must have the same length")
    if I_loops.size not in (1, n_loops):
        raise ValueError("I_loops must be scalar or same length as R_loops")

    # sum contributions
    for k in range(n_loops):
        I_k = float(I_loops[k]) if I_loops.size == n_loops else float(I_loops)
        if pointwise:
            # try vectorized call on arrays of points
            try:
                Br, Bz, Bphi = B_loop_cyl(R_pts, Z_pts, R_loops[k], Z_loops[k], I_k)
                Br = np.asarray(Br).ravel()
                Bz = np.asarray(Bz).ravel()
                Bphi = np.asarray(Bphi).ravel()
            except Exception:
                # fallback to per-point scalar calls
                Br = np.zeros(R_pts.size, dtype=float)
                Bz = np.zeros(R_pts.size, dtype=float)
                Bphi = np.zeros(R_pts.size, dtype=float)
                for i in range(R_pts.size):
                    br, bz, bp = B_loop_cyl(float(R_pts[i]), float(Z_pts[i]), R_loops[k], Z_loops[k], I_k)
                    Br[i] = br; Bz[i] = bz; Bphi[i] = bp
            B_R_tot += Br
            B_Z_tot += Bz
            B_phi_tot += Bphi
        else:
            # grid evaluation (prefer vectorized B_loop_cyl)
            Br, Bz, Bphi = B_loop_cyl(R_grid, Z_grid, R_loops[k], Z_loops[k], I_k)
            B_R_tot += Br
            B_Z_tot += Bz
            B_phi_tot += Bphi

    return B_R_tot, B_Z_tot, B_phi_tot
# ...existing code...
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


def add_constant_Bz(R_grid, Bz_const):
    B_z_added = np.zeros_like(R_grid)
    B_r_added = np.zeros_like(R_grid)
    psi_added = np.zeros_like(R_grid)
    psi_added[:, :] = 1 / 2 * Bz_const * R_grid[:, :] ** 2
    return B_r_added, B_z_added, psi_added


def plot_tomographic(path,ax,emi_vmin=0, emi_vmax=3.e20):
    
    data = loadmat(path, squeeze_me=True, struct_as_record=False)

    inv_grid=data['inv_grid']

    tri_x = inv_grid.tri_x
    tri_y = inv_grid.tri_y
    tri_nodes = inv_grid.tri_nodes
    R_values = inv_grid.R_values
    Z_values = inv_grid.Z_values
    R_mean = inv_grid.R_mean
    Z_mean = inv_grid.Z_mean
    Area = inv_grid.Area

    emi = data['emi']
    t = data['t']
    time = data['time']
    time_vec = np.asarray(time)
    t_idx = np.argmin(np.abs(time_vec - t))
    emi = emi[:, t_idx]
    triangles = tri_nodes.astype(int)  # tri_nodes doit être zéro-indexé
    tri = Triangulation(tri_x, tri_y, triangles)

    tpc = ax.tripcolor(tri, emi, shading='flat', edgecolors='none', vmin=emi_vmin, vmax=emi_vmax)
    
    return tpc,tri,emi


def plot_LIUQE(path,n_levels=15,**kwargs):


   fig, ax, kwargs = create_canvas(**kwargs)

   colors = kwargs.pop("colors", ["xkcd:cyan", "xkcd:cyan", "black"])

   labels = kwargs.pop("labels", ["LIUQE", None, "LCFS"])

   lw = kwargs.pop("lw", [0.8, 0.8, 1])

   linestyles = kwargs.pop("linestyles", ['-', '-', '-'])

   data = loadmat(path, squeeze_me=True, struct_as_record=False)

   flux=data['flux']
   r = flux.r
   z = flux.z
   psi = flux.psi.T
   psi_r = flux.psi_r
   psi_z = flux.psi_z.T
   psi_axis = flux.psi_axis

   levels1 = np.linspace(-flux.psi_axis, 0, n_levels)
   levels2 = np.linspace(flux.psi_axis, 0, n_levels)

   levels1.sort()
   levels2.sort()

   ax.contour(psi_r, psi_z, psi, levels=levels1,colors=colors[0], linestyles=linestyles[0], linewidths=lw[0], label=labels[0],**kwargs)
   ax.contour(psi_r, psi_z, psi,  levels=levels2, colors=colors[1], linestyles=linestyles[1], linewidths=lw[1], label=labels[1],**kwargs)
   ax.contour(psi_r, psi_z, psi,   levels=[0],   colors=colors[2], linestyles=linestyles[2], linewidths=lw[2], label=labels[2],**kwargs)

   lgd1 = Line2D([0], [0], color=colors[2], lw=lw[2], linestyle=linestyles[2])   
   lgd2 = Line2D([0], [0], color=colors[1], lw=lw[0], linestyle=linestyles[0])

   return data, lgd1, lgd2


def plot_Br_Bz(JFField, n_levels=8, ax_pair=None, R_count=56, Z_count=130, phi=0,
               vessel=None, Rt=None, Zt=None, cmap='viridis', alpha_map=0.7):
    """
    Trace B_R et B_Z (en mT) côte à côte pour un AxisymmetricCylindricalGridField.
    Args:
      JFField : AxisymmetricCylindricalGridField (doit avoir .R .Z et .pertfield with B_R/B_Z)
      n_levels: int, nombre de niveaux de contour
      ax_pair : tuple(ax1,ax2) ou None -> crée une figure si None
      R_count, Z_count : résolution du maillage pour le tracé
      phi : valeur de phi (scalare) utilisée pour l'évaluation du champ
      vessel, Rt, Zt : (optionnel) géométrie à dessiner par dessus
      cmap : colormap pour contourf
    Returns:
      fig, (ax1, ax2), (BR, BZ)
    """


    # create axes if not provided
    if ax_pair is None:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.0, 6.6))
    else:
        ax1, ax2 = ax_pair
        fig = ax1.figure

    # grid extents from field
    R_min, R_max = float(np.min(JFField.R)), float(np.max(JFField.R))
    Z_min, Z_max = float(np.min(JFField.Z)), float(np.max(JFField.Z))

    R_vec = np.linspace(R_min, R_max, num=R_count)
    Z_vec = np.linspace(Z_min, Z_max, num=Z_count)

    R_grid2, Z_grid2 = np.meshgrid(R_vec, Z_vec, indexing='xy')

    # evaluate field robustly: try vectorized then fallback to loops
    BR = np.zeros(R_grid2.shape)
    BZ = np.zeros(R_grid2.shape)

    pts = np.column_stack((R_grid2.ravel(), np.full(R_grid2.size, phi), Z_grid2.ravel()))
    try:
        BR_flat = JFField.pertfield.B_R(pts)
        BZ_flat = JFField.pertfield.B_Z(pts)
        BR = np.asarray(BR_flat).reshape(R_grid2.shape)
        BZ = np.asarray(BZ_flat).reshape(Z_grid2.shape)
    except Exception:
        # fallback point by point
        for j in range(Z_grid2.shape[0]):
            for i in range(R_grid2.shape[1]):
                Rv = float(R_grid2[j, i]); Zv = float(Z_grid2[j, i])
                try:
                    BR[j, i] = float(JFField.pertfield.B_R([Rv, phi, Zv]))
                    BZ[j, i] = float(JFField.pertfield.B_Z([Rv, phi, Zv]))
                except Exception:
                    # try single-point tuple if other signatures required
                    BR[j, i] = float(JFField.pertfield.B_R([Rv, phi, Zv]))
                    BZ[j, i] = float(JFField.pertfield.B_Z([Rv, phi, Zv]))

    # convert to mT for display
    BR = 1000.0 * BR
    BZ = 1000.0 * BZ

    # optional draw vessel / tips if provided
    if vessel is not None or (Rt is not None and Zt is not None):
        # support passing vessel struct or arrays Rt,Zt
        if vessel is not None:
            R_in = np.asarray(vessel.R_in).ravel()
            Z_in = np.asarray(vessel.Z_in).ravel()
            Rt_ = np.asarray(vessel.Rt).ravel()
            Zt_ = np.asarray(vessel.Zt).ravel()
        else:
            R_in = np.asarray(Rt).ravel() * 0  # empty outer ring
            Z_in = np.asarray(Zt).ravel() * 0
            Rt_ = np.asarray(Rt).ravel()
            Zt_ = np.asarray(Zt).ravel()

        # draw tip curve on both axes
        # for ax in (ax1, ax2):
        #     ax.plot(Rt_, Zt_, color='k', label='_nolegend_')
         
        

        R_in_anti=R_in[::-1]
        Z_in_anti=Z_in[::-1]

        for ax in (ax1, ax2):
            fill_between_polygons_shapely(ax, R_in_anti, Z_in_anti,Rt_,Zt_) 
            line = ax.plot(Rt_, Zt_, color='k', label='_nolegend_')[0]
            line.set_picker(False)
            line.set_zorder(2)
        # # create a proper ring (outer boundary then inner boundary reversed) to avoid covering whole interior
        #     inner = np.column_stack((R_in, Z_in))
        #     outer = np.column_stack((Rt_, Zt_))
        # # detect which curve is outer by mean radius (fallback to given order)
        #     try:
        #       if np.mean(outer[:, 0]) < np.mean(inner[:, 0]):
        #          outer, inner = inner, outer
        #     except Exception:
        #        pass
        # # build polygon verts: outer (ccw), then inner reversed (cw) to form a ring
        #     verts = np.vstack([outer, inner[::-1]])
        #     patch = MplPolygon(verts, closed=True, facecolor=(0.5, 0.5, 0.5), edgecolor='none', label='_nolegend_', zorder=0)
        #     patch.set_picker(False)
        #     ax.add_patch(patch)

   
    # plotting helper
    def _plot_field(ax, data, alpha=alpha_map):# title):
        bz_min, bz_max = np.nanmin(data), np.nanmax(data)
        if np.isfinite(bz_min) and np.isfinite(bz_max) and bz_max > bz_min:
            levels = np.linspace(bz_min, bz_max, n_levels)
        else:
            levels = np.linspace(-1.0, 1.0, n_levels)

        cs = ax.contour(R_grid2, Z_grid2, data, levels=levels, colors='k',
                        linewidths=1.2, linestyles='--')
        ax.clabel(cs, fmt='%.1f', inline=True, fontsize=12)
        cf = ax.contourf(R_grid2, Z_grid2, data, levels=levels, cmap=cmap, alpha=alpha_map)

        # colorbar on the right
        #divider = make_axes_locatable(ax)
        #cax = divider.append_axes("right", size="3%", pad=0.06)
        #ticks = np.linspace(levels[0], levels[-1], len(levels))
        #cbar = fig.colorbar(cf, cax=cax, orientation='vertical', ticks=ticks)
        #cbar.ax.set_yticklabels([f"{t:.0f}" for t in ticks])
        #cbar.ax.yaxis.get_offset_text().set_visible(False)
        #cbar.set_label(f"{title} [mT]")

        ax.set_xlim(R_min, R_max)
        ax.set_ylim(Z_min, Z_max)
        ax.set_aspect('equal')

    _plot_field(ax2, BZ, alpha=alpha_map)# None)
    _plot_field(ax1, BR, alpha=alpha_map)# None)

    ax1.set_xlabel(r"$R\ [m]$",fontsize=11); ax1.set_ylabel(r"$Z\ [m]$",fontsize=11)
    ax2.set_xlabel(r"$R\ [m]$",fontsize=11)
    ax1.set_title("$B_R$[mT]",fontsize=11); ax2.set_title("$B_Z$[mT]",fontsize=11)

    fig.tight_layout()
    return fig, (ax1, ax2), (BR, BZ)
# ...existing code...
# ...existing code...
def plot_B_parallel_along_line_from_loops(a, Z0, I_loop, p0, p1, n=100, phi=1.9,
                                          ax=None, plot=True, return_data=False,ax_line=None):
    """
    Compute and (optionally) plot the magnetic field component parallel to the straight
    line from p0 -> p1 using the field produced by B_sum_of_loops(a, Z0, I_loop).
    Args:
      a, Z0     : arrays of loop positions (same format as used by B_sum_of_loops)
      I_loop    : current per loop (scalar or array)
      p0, p1    : (R,Z) tuples or lists in data coordinates
      n         : number of sample points along the segment
      phi       : toroidal angle (not used by B_sum_of_loops but kept for API parity)
      ax        : matplotlib Axes to plot the scalar vs arc-length (if None a new fig/ax is created)
      plot      : if True draw the plot on ax
      return_data: if True return (s, B_par) arrays
    Returns:
      if return_data: (s, B_par) where s is arc length along segment and B_par is B·t (in T)
      otherwise returns (fig, ax) if plot created a new figure, or None if plotted on provided ax.
    """


    # build sampling points along the segment in (R,Z)
    R0, Z0p = float(p0[0]), float(p0[1])
    R1, Z1p = float(p1[0]), float(p1[1])
    Rs = np.linspace(R0, R1, n)
    Zs = np.linspace(Z0p, Z1p, n)

    # arc-length parameter s
    ds = np.sqrt((np.diff(Rs))**2 + (np.diff(Zs))**2)
    s = np.concatenate([[0.0], np.cumsum(ds)])

    # try vectorized evaluation of the loop field
    BR = None
    BZ = None
    try:
        # B_sum_of_loops accepts grid-like R, Z. Passing 1D arrays should return arrays.
        BR_try, BZ_try, Bphi_try = B_sum_of_loops(Rs, Zs, a, Z0, I_loop)
        BR = np.asarray(BR_try).ravel()
        BZ = np.asarray(BZ_try).ravel()
        if BR.size != n or BZ.size != n:
            # try reshape if returned shape (n,1) or (1,n)
            BR = BR.reshape(-1)[:n]
            BZ = BZ.reshape(-1)[:n]
    except Exception:
        # fallback: evaluate point-by-point
        BR = np.zeros(n)
        BZ = np.zeros(n)
        for i in range(n):
            try:
                br, bz, bp = B_sum_of_loops(np.array([Rs[i]]), np.array([Zs[i]]), a, Z0, I_loop)
                BR[i] = np.asarray(br).ravel()[0]
                BZ[i] = np.asarray(bz).ravel()[0]
            except Exception as e_point:
                raise RuntimeError(f"Failed to evaluate B_sum_of_loops at point {i} (R={Rs[i]},Z={Zs[i]}): {e_point}")

    # tangent unit vector in RZ plane (same for whole straight line)
    t = np.array([R1 - R0, Z1p - Z0p], dtype=float)
    norm_t = np.hypot(t[0], t[1])
    if norm_t == 0:
        raise ValueError("p0 and p1 are identical; tangent undefined.")
    t_unit = t / norm_t

    # compute B_par = B · t_unit using R and Z components (BR, BZ)
    B_par = BR * t_unit[0] + BZ * t_unit[1]  # in Tesla

    # plotting
    created_fig_ax = None
    if plot:
        if ax is None:
            fig, axp = plt.subplots(1, 1, figsize=(6, 3.5))
            created_fig_ax = (fig, axp)
        else:
            axp = ax
        axp.plot(s, B_par, '-', color='C1', lw=1.5, label='B$_{||}$ loops')
        #axp.set_xlabel("arc length s [m]")
        axp.set_ylabel(r"$B_{\parallel}$ [T]")
        axp.grid(True, ls='--', alpha=0.4)
    if ax_line is not None:
        ax_line.plot([R0, R1], [Z0p, Z1p], color='black', linewidth=2, zorder=10)
    if return_data:
        return s, B_par
    if created_fig_ax:
        return created_fig_ax
    return None
# ...existing code...

# ...existing code...
def plot_B_parallel_along_line(JFField, p0, p1, n=100, phi=1.9, ax=None, plot=True, return_data=False):
    """
    Compute and (optionally) plot the magnetic field component parallel to the straight
    line from p0 -> p1. Points sampled = n.
    Args:
      JFField : field object (AxisymmetricCylindricalGridField)
      p0, p1  : (R,Z) tuples or lists in data coordinates
      n       : number of sample points along the segment
      phi     : toroidal angle used to evaluate the field
      ax      : matplotlib Axes to plot the scalar vs arc-length (if None a new fig/ax is created)
      plot    : if True draw the plot on ax
      return_data : if True return (s, Bpar) arrays
    Returns:
      if return_data: (s, Bpar) where s is arc length along segment and Bpar is B·t (in T)
      otherwise returns (fig, ax) if plot created a new figure, or None if plotted on provided ax.
    Notes:
      - tries to call JFField.B(pts) (vectorized) first. Falls back to component calls
        (B_R, B_phi, B_Z) on JFField or JFField.pertfield if necessary.
      - uses only R and Z components (assumes vector ordering [BR, Bphi, BZ]).
    """
    import numpy as np
    import matplotlib.pyplot as plt

    # build sampling points along the segment in (R,Z)
    R0, Z0 = float(p0[0]), float(p0[1])
    R1, Z1 = float(p1[0]), float(p1[1])
    Rs = np.linspace(R0, R1, n)
    Zs = np.linspace(Z0, Z1, n)

    # arc-length parameter s
    ds = np.sqrt((np.diff(Rs))**2 + (np.diff(Zs))**2)
    s = np.concatenate([[0.0], np.cumsum(ds)])

    # prepare points for field calls: shape (n,3) as (R,phi,Z)
    pts = np.column_stack((Rs, np.full(n, phi), Zs))
    
    

    B = np.zeros_like(pts)
    for i in range(n):
        B[i,:]= JFField.B(pts[i,:])
    

    # now B is (n,3). tangent unit vector in RZ plane (same for whole straight line)
    t = np.array([R1 - R0, Z1 - Z0], dtype=float)
    norm_t = np.hypot(t[0], t[1])
    if norm_t == 0:
        raise ValueError("p0 and p1 are identical; tangent undefined.")
    t_unit = t / norm_t
    # compute B_par = B · t_unit using R and Z components (B[:,0], B[:,2])
    B_par = B[:, 0] * t_unit[0] + B[:, 2] * t_unit[1]  # in Tesla

    # plotting
    created_fig_ax = None
    if plot:
        if ax is None:
            fig, axp = plt.subplots(1, 1, figsize=(6, 3.5))
            created_fig_ax = (fig, axp)
        else:
            axp = ax

        axp.plot(s, B_par, '-', color='C0', lw=1.5,label='B$_{||}$ from field ')
        #axp.set_xlabel("arc length s [m]")
       # axp.set_ylabel(r"$B_{\parallel}$ [T]")
        axp.grid(True, ls='--', alpha=0.4)

    if return_data:
        return s, B_par
    if created_fig_ax:
        return created_fig_ax
    return None
# ...existing code...
    


def plot_vessel(ax):
        data_vessel = loadmat('./script/jellyfisch/75979_12/mat_files/Liuqe_JF_75979_120_with_V.mat', squeeze_me=True, struct_as_record=False)

        vessel=data_vessel['vessel']
        R_in = np.asarray(vessel.R_in).ravel()
        Z_in = np.asarray(vessel.Z_in).ravel()
        Rt = np.asarray(vessel.Rt).ravel()
        Zt = np.asarray(vessel.Zt).ravel()

        R_in_anti=R_in[::-1]
        Z_in_anti=Z_in[::-1]

        fill_between_polygons_shapely(ax, R_in_anti, Z_in_anti,Rt,Zt) 
        line = ax.plot(Rt, Zt, color='k', label='_nolegend_')[0]
        line.set_picker(False)
        line.set_zorder(2)
        