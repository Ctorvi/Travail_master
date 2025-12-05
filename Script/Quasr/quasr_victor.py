from pyoculus.fields import SimsoptBfield
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
import numpy as np
import logging
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
from simsopt import load


plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.size": 10,
    }
)
logging.basicConfig(level=logging.DEBUG)

fig = plt.figure(figsize=(12, 8))

gs = GridSpec(1, 3, wspace=0.2, width_ratios=[1, 1, 1])
ax3 = fig.add_subplot(gs[0, 2])
ax2 = fig.add_subplot(gs[0, :2])



QUASR_ID = 229079
IOTA_ITS = 40
RECOMPUTE_PPLOT = False
COMPUTE_TURNSTILE_AREAS_FILE = False

logging.basicConfig(level=logging.DEBUG)

jsonname = 'script/quasr/serial0229079.json'
surfaces, coils = load(jsonname)


### AX2: Poincare plot

AXISGUESS = [.868, 0.]
simsoptfield = SimsoptBfield.from_coils(coils, Nfp=3, interpolate=False)  #, surf=s, n=50)
fieldlinemap = CylindricalBfieldSection.without_axis(simsoptfield, guess=AXISGUESS, rtol=1e-10, nsteps=1e5)
fieldlinemap_poincare = CylindricalBfieldSection.without_axis(simsoptfield, guess=AXISGUESS, rtol=1e-7, nsteps=1e5)



#pplot = PoincarePlot.with_horizontal(fieldlinemap_poincare, 0.5, 60)
#pplot.compute(400)
#np.save('quasr_hits', pplot._hits)
Hits=np.load('quasr_hits.npy', allow_pickle=True)

ax2.scatter(Hits[:,:, 0], Hits[:,:, 1], color="xkcd:dark grey", s=1.4, linewidths=0)
#ax2.set_aspect('equal')
# fp8xa = FixedPoint(fieldlinemap)
# fp8xa.find(t=8, guess=[1.14374773, 0.0203871])
# fp8xa.m = 8
# fp8xa.plot(ax=ax2, zorder=15, s=7, color='xkcd:blue')

# fp8xb = FixedPoint(fieldlinemap)
# fp8xb.find(t=8, guess=[1.13535758, 0.07687874])
# fp8xb.m = 8
# fp8xb.plot(ax=ax2, zorder=15, s=7, color='xkcd:light blue')




# colors8a = ['xkcd:light blue', 'xkcd:royal blue']
# manifold8a_in = Manifold(fieldlinemap, fp8xa, fp8xb, first_stable=True)
# #manifold8a_in.compute(eps_s=1e-4, eps_u=1e-4, nint_s=11, nint_u=11, neps_s=10, neps_u=10)
# #manifold8a_in.save_mf_quasr('./script/quasr/manifolds/8a_in.pkl')

# manifold8a_in.load_mf_quasr("./script/quasr/manifolds/8a_in.pkl")

# manifold8a_in.plot(ax=ax2, markersize=0, lw=0.5, colors=colors8a)
# manifold8a_in.plot_manifold_copies(ax=ax2, markersize=0, lw=0.5, colors=colors8a)




# # print('calculating manifold8a_out')
# manifold8a_out = Manifold(fieldlinemap, fp8xa, fp8xb, first_stable=False)
# #manifold8a_out.compute(eps_s=1e-4, eps_u=1e-4, nint_s=10, nint_u=10, neps_s=10, neps_u=10)
# #manifold8a_out.save_mf_quasr('./script/quasr/manifolds/8a_out.pkl')


# manifold8a_out.load_mf_quasr("./script/quasr/manifolds/8a_out.pkl")

# manifold8a_out.plot(ax=ax2, markersize=0, lw=0.5, colors=colors8a)
# manifold8a_out.plot_manifold_copies(ax=ax2, markersize=0, lw=0.5, colors=colors8a)

# fp8xb2 = FixedPoint(fieldlinemap)
# fp8xb2.find(t=8, guess=fp8xb.coords[1])


# #### NEED TO SWITCH TO Z=O OCCURRENCE AND RE-FIDDLE MANIFOLDS
# colors8b = ['xkcd:lavender', 'xkcd:bluey purple']
# print('calculating manifold8b_in')
# manifold8b_in = Manifold(fieldlinemap, fp8xa, fp8xb2, first_stable=False)


# #manifold8b_in.compute(eps_s=1e-4, eps_u=1e-4, nint_s=11, nint_u=11, neps_s=10, neps_u=10)
# #manifold8b_in.save_mf_quasr('./script/quasr/manifolds/8b_in.pkl')


# manifold8b_in.load_mf_quasr("./script/quasr/manifolds/8b_in.pkl")

# manifold8b_in.plot(ax=ax2, markersize=0, lw=0.5, colors=colors8b)
# manifold8b_in.plot_manifold_copies(ax=ax2, markersize=0, lw=0.5, colors=colors8b)


# print('calculating manifold8b_out')


# manifold8b_out = Manifold(fieldlinemap, fp8xa, fp8xb2, first_stable=True)


# #manifold8b_out.compute(eps_s=1e-4, eps_u=1e-4, nint_s=11, nint_u=10, neps_s=10, neps_u=11)
# #manifold8b_out.save_mf_quasr('./script/quasr/manifolds/8b_out.pkl')


# manifold8b_out.load_mf_quasr("./script/quasr/manifolds/8b_out.pkl")

# manifold8b_out.plot(ax=ax2, markersize=0, lw=0.5, colors=colors8b)
# manifold8b_out.plot_manifold_copies(ax=ax2, markersize=0, lw=0.5, colors=colors8b)

# fp9x = FixedPoint(fieldlinemap)
# fp9x.find(t=9, guess=[1.203, .05])
# fp9x.m = 9
# #fp9x.plot(ax=ax2, zorder=15, s=7, color='xkcd:green')

# fp9x2 = FixedPoint(fieldlinemap)
# fp9x2.find(t=9, guess=fp9x.coords[1])

# colors9 = ['xkcd:bright green', 'xkcd:dark green']






# print('calculating manifold9_in')

# manifold9_in = Manifold(fieldlinemap, fp9x, fp9x2, first_stable=False)
# #manifold9_in.compute(eps_s=1e-4, eps_u=1e-4, nint_s=11, nint_u=11, neps_s=8, neps_u=10)

# #manifold9_in.save_mf_quasr('./script/quasr/manifolds/9_in.pkl')


# manifold9_in.load_mf_quasr("./script/quasr/manifolds/9_in.pkl")

# manifold9_in.plot(ax=ax2, markersize=0, lw=0.5, colors=colors9)
# manifold9_in.plot_manifold_copies(ax=ax2, markersize=0, lw=0.5, colors=colors9)





# print('calculating manifold9_out')

# manifold9_out = Manifold(fieldlinemap, fp9x, fp9x2, first_stable=True)

# #manifold9_out.compute(eps_s=1e-4, eps_u=1e-4, nint_s=12, nint_u=11, neps_s=10, neps_u=10)
# #manifold9_out.save_mf_quasr('./script/quasr/manifolds/9_out.pkl')


# manifold9_out.load_mf_quasr("./script/quasr/manifolds/9_out.pkl")

# manifold9_out.plot(ax=ax2, markersize=0, lw=0.5, colors=colors9)

# manifold9_out.plot_manifold_copies(ax=ax2, markersize=0, lw=0.5, colors=colors9)




# #fp9o = FixedPoint(fieldlinemap)
# #fp9o.find(t=9, guess=[ 1.203, .00])
# #fp9o.plot(ax=ax2, zorder=15, s=10, color='xkcd:light orange')
# fp10x = FixedPoint(fieldlinemap)
# fp10x.find(t=10, guess=[1.253, .00])
# #fp10x.plot(ax=ax2, zorder=15, s=7, color='xkcd:goldenrod')
# fp10x.m = 10
# fp10x2 = FixedPoint(fieldlinemap)
# fp10x2.find(t=10, guess=fp10x.coords[1])

# colors10 = ['xkcd:goldenrod', 'xkcd:coral']






# manifold10_in = Manifold(fieldlinemap, fp10x, fp10x2, first_stable=False)
# #manifold10_in.compute(nint_s=4, nint_u=4, neps_s=40, neps_u=40, eps_s=2e-4, eps_u=2e-4)
# #manifold10_in.save_mf_quasr('./script/quasr/manifolds/10_in.pkl')

# manifold10_in.load_mf_quasr("./script/quasr/manifolds/10_in.pkl")
# manifold10_in.plot(ax=ax2, markersize=0, lw=0.5, colors=colors10)
# manifold10_in.plot_manifold_copies(ax=ax2, markersize=0, lw=0.5, colors=colors10)







# manifold10_out = Manifold(fieldlinemap, fp10x, fp10x2, first_stable=True)
# #manifold10_out.compute(nint_s=4, nint_u=4, neps_s=40, neps_u=40, eps_s=3e-4, eps_u=3e-4)


# #manifold10_out.save_mf_quasr('./script/quasr/manifolds/10_out.pkl')
# manifold10_out.load_mf_quasr("./script/quasr/manifolds/10_out.pkl")

# manifold10_out.plot(ax=ax2, markersize=0, lw=0.5, colors=colors10)
# manifold10_out.plot_manifold_copies(ax=ax2, markersize=0, lw=0.5, colors=colors10)






# fp11 = FixedPoint(fieldlinemap)
# fp11.find(t=11, guess=[1.29093278, 0.04675814])
# #fp11.plot(ax=ax2, zorder=15, s=7, color='xkcd:hot pink')
# fp11.m = 11
# fp11b = FixedPoint(fieldlinemap)
# fp11b.find(t=11, guess=fp11.coords[1])





# colors11 = ['xkcd:pink', 'xkcd:fuchsia']



# manifold11_in = Manifold(fieldlinemap, fp11, fp11b, first_stable=False)
# #manifold11_in.compute(nint_s=4, nint_u=4, neps_s=20, neps_u=20, eps_s=3e-5, eps_u=3e-5)
# #manifold11_in.save_mf_quasr('./script/quasr/manifolds/11_in.pkl')


# manifold11_in.load_mf_quasr("./script/quasr/manifolds/11_in.pkl")

# manifold11_in.plot(ax=ax2, markersize=0, lw=0.5, colors=colors11)
# manifold11_in.plot_manifold_copies(ax=ax2, markersize=0, lw=0.5, colors=colors11)







# manifold11_out = Manifold(fieldlinemap, fp11, fp11b, first_stable=True)
# #manifold11_out.compute(nint_s=4, nint_u=4, neps_s=20, neps_u=20, eps_s=3e-5, eps_u=3e-5)
# #manifold11_out.save_mf_quasr('./script/quasr/manifolds/11_out.pkl')

# manifold11_out.load_mf_quasr("./script/quasr/manifolds/11_out.pkl")

# manifold11_out.plot(ax=ax2, markersize=0, lw=0.5, colors=colors11)
# manifold11_out.plot_manifold_copies(ax=ax2, markersize=0, lw=0.5, colors=colors11)

# fpcrazy12o = FixedPoint(fieldlinemap)
# fpcrazy12o.find(t=12, guess=[1.495, .00])
# #fpcrazy12o.plot(ax=ax2, zorder=15, s=10, color='xkcd:burgundy')

# fpcrazy12x = FixedPoint(fieldlinemap)
# fpcrazy12x.find(t=12, guess=[1.48976075, 0.09298623])
# fpcrazy12x.m = 12
# #fpcrazy12x.plot(ax=ax2, zorder=15, s=7, color='xkcd:burgundy')

# fpcrazy12x2 = FixedPoint(fieldlinemap)
# fpcrazy12x2.find(t=12, guess=fpcrazy12x.coords[1])




# crazymanifold = Manifold(fieldlinemap, fpcrazy12x, fpcrazy12x2, first_stable=True)
# #crazymanifold.compute(nint_s=5, nint_u=5, neps_s=20, neps_u=20, eps_s=1e-4, eps_u=1e-4)
# crazymanifold.save_mf_quasr('./script/quasr/manifolds/crazymanifold.pkl')

# manifold11_out.load_mf_quasr("./script/quasr/manifolds/crazymanifold.pkl")

# crazymanifold.plot(ax=ax3, markersize=0, lw=0.5)
# crazymanifold.plot_manifold_copies(ax=ax3, markersize=0, lw=0.5)

# #pplot.plot(ax=ax2, color='xkcd:grey', s=0.7, linewidths=0)
ax2.set_aspect('equal')
ax2.set_xlim(.55, 1.6)
ax2.set_ylim(-.3, .3)
ax2.set_xlabel(r'$R [m]$')
ax2.set_ylabel(r'$Z [m]$')
# # do fixed points



# ############ inset ax3 #####3
# manifolds = {
#     "manifold8a_in": manifold8a_in,
#     "manifold8a_out": manifold8a_out,
#     "manifold8b_in": manifold8b_in,
#     "manifold8b_out": manifold8b_out,
#     "manifold9_in": manifold9_in,
#     "manifold9_out": manifold9_out,
#     "manifold10_in": manifold10_in,
#     "manifold10_out": manifold10_out,
#     "manifold11_in": manifold11_in,
#     "manifold11_out": manifold11_out,
# }
# greycolors = ['xkcd:grey', 'xkcd:light grey']
# #pplot.plot(ax=ax3, color='xkcd:grey', s=0.7, linewidths=0)
# manifold8a_in.plot(ax=ax3, markersize=0, lw=0.5, colors=colors8a)
# manifold8a_out.plot(ax=ax3, markersize=0, lw=0.5, colors=colors8a)
# manifold8b_in.plot(ax=ax3, markersize=0, lw=0.5, colors=colors8b)
# manifold8b_out.plot(ax=ax3, markersize=0, lw=0.5, colors=colors8b)
# manifold9_in.plot(ax=ax3, markersize=0, lw=0.5, colors=colors9)
# manifold9_out.plot(ax=ax3, markersize=0, lw=0.5, colors=colors9)
# manifold10_in.plot(ax=ax3, markersize=0, lw=0.5, colors=colors10)
# manifold10_out.plot(ax=ax3, markersize=0, lw=0.5, colors=colors10)
# manifold11_in.plot(ax=ax3, markersize=0, lw=0.5, colors=colors11)
# manifold11_out.plot(ax=ax3, markersize=0, lw=0.5,  colors=colors11)
# manifold8a_in.plot_manifold_copies(ax=ax3, markersize=0, lw=0.5, colors=greycolors)
# manifold8a_out.plot_manifold_copies(ax=ax3, markersize=0, lw=0.5, colors=greycolors)
# manifold8b_in.plot_manifold_copies(ax=ax3, markersize=0, lw=0.5, colors=greycolors)
# manifold8b_out.plot_manifold_copies(ax=ax3, markersize=0, lw=0.5, colors=greycolors)
# manifold9_in.plot_manifold_copies(ax=ax3, markersize=0, lw=0.5, colors=greycolors)
# manifold9_out.plot_manifold_copies(ax=ax3, markersize=0, lw=0.5, colors=greycolors)
# manifold10_in.plot_manifold_copies(ax=ax3, markersize=0, lw=0.5, colors=greycolors)
# manifold10_out.plot_manifold_copies(ax=ax3, markersize=0, lw=0.5, colors=greycolors)
# manifold11_in.plot_manifold_copies(ax=ax3, markersize=0, lw=0.5, colors=greycolors)
# manifold11_out.plot_manifold_copies(ax=ax3, markersize=0, lw=0.5,  colors=greycolors)
# fp8xa.plot(ax=ax3, zorder=15, s=10, color='xkcd:blue')
# fp8xb.plot(ax=ax3, zorder=15, s=10, color='xkcd:light blue')
# fp9x.plot(ax=ax3, zorder=15, s=10, color='xkcd:green')
# fp10x.plot(ax=ax3, zorder=15, s=10, color='xkcd:goldenrod')
# fp11.plot(ax=ax3, zorder=15, s=10, color='xkcd:hot pink')
# #inmanifold.plot(ax=ax3, markersize=0, lw=0.5)
# #outmanifold.plot(ax=ax3, markersize=0, lw=0.5)
#ax3.set_xlim(.827, 1.152)
ax3.set_aspect('equal')
ax3.set_xlabel(r'$R [m]$')
ax3.set_ylabel(r'$Z [m]$')
ax3.set_xlim(1.1, 1.3)
ax3.set_ylim(-.1, .1)

# # #

# # ######### calculate rotational transform ######

# #iotapplot = PoincarePlot.with_horizontal(fieldlinemap_poincare, 0.7, 60)
# #iotacompute_or_load(iotapplot, f'iotaplot_{QUASR_ID}.npy', 20, recompute_anyway=RECOMPUTE_PPLOT)

# #ax4.plot(iotapplot.rho, -1*iotapplot.iota)
# ax4.set_xlabel(r'$R-R_0 [m]$')
# ax4.set_ylabel(r'$\imath$')


# # add axhlines at 3/8, 3/9, 3/10, 3/11, 3/12:
# ax4.axhline(3/8, color='xkcd:light blue', linestyle='--', lw=0.5)
# ax4.axhline(3/9, color='xkcd:green', linestyle='--', lw=0.5)
# ax4.axhline(3/10, color='xkcd:goldenrod', linestyle='--', lw=0.5)
# ax4.axhline(3/11, color='xkcd:fuchsia', linestyle='--', lw=0.5)
# ax4.axhline(3/12, color='xkcd:burgundy', linestyle='--', lw=0.5)

# # add text for each line on ax4:
# ax4.text(0.53, 3/8, r'$\frac{3}{8}$', va='center')
# ax4.text(0.55, 3/9, r'$\frac{3}{9}$', va='center')
# ax4.text(0.58, 3/10, r'$\frac{3}{10}$', va='center')
# ax4.text(0.61, 3/11, r'$\frac{3}{11}$', va='center')
# ax4.text(0.63, 3/12, r'$\frac{3}{12}$', va='center')

# plt.savefig('quasr_manyisland.png', bbox_inches='tight', dpi=720)
# #plt.ion()
plt.show()

# if COMPUTE_TURNSTILE_AREAS_FILE:

#     # Open a file to write the LaTeX table
#     with open("quasr_manyisland_manifold_turnstile_areas.tex", "w") as f:
#         # Write the table header
#         f.write("\\begin{tabular}{lcc}\n")
#         f.write("\\hline\n")
#         f.write("Manifold & Turnstile Area (in) & Turnstile Area (out) \\\\\n")
#         f.write("\\hline\n")

#         # Loop over the manifold pairs
#         for name, manifold in manifolds.items():
#             if name.endswith("_in"):
#                 out_name = name.replace("_in", "_out")
#                 manifold_in = manifold
#                 manifold_out = manifolds[out_name]

#                 if not len(manifold_in.clinics)==2:
#                     # Compute turnstile areas
#                     manifold_in.find_clinic_single(guess_eps_s=1e-4, guess_eps_u=1e-4, root_args={'options': {'factor': 1e-3}}, nretry=3, reset_clinics=True)
#                     manifold_in.find_other_clinic(shift_in_stable=0.5, nretry=3)
#                 area_in = manifold_in.compute_turnstile_areas()
                
#                 if not len(manifold_out.clinics)==2:
#                     manifold_out.find_clinic_single(guess_eps_s=1e-4, guess_eps_u=1e-4, root_args={'options': {'factor': 1e-3}}, nretry=3, reset_clinics=True)
#                     manifold_out.find_other_clinic(shift_in_stable=0.5, nretry=3)
#                 area_out = manifold_out.compute_turnstile_areas()

#                 # Write the results to the LaTeX table
#                 f.write(f"{name[:-3]} & ${area_in[0]:.4e}$ & ${area_out[0]:.4e}$ \\\\\n")

#         # Write the table footer
#         f.write("\\hline\n")
#         f.write("\\end{tabular}\n")

#     logging.info("LaTeX table written to manifold_turnstile_areas.tex")
