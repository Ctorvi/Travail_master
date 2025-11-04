from pyoculus.fields import AxisymmetricCylindricalGridField
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
import numpy as np
from matplotlib import pyplot as plt
import logging
from scipy.io import loadmat

fig, ax = plt.subplots(1, 1, figsize=(6, 6.6))

plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.size": 10,
    }
)
logging.basicConfig(level=logging.DEBUG)

JFField = AxisymmetricCylindricalGridField.from_matlab_file('./Script/Snowflake/Snowflake_09.mat', with_perturbation=False)



#### Liuqe Snowflake data load ####






section = CylindricalBfieldSection(JFField,R0=0.88, Z0=0)



top_o = FixedPoint(section)
top_o.find(1, [0.9,0.0], method='scipy.root')
top_o_coord = top_o.coords[0]
top_o.plot(ax=ax, marker='o', color="xkcd:crimson")


x_point= FixedPoint(section)
x_point.find(1, [0.75,-0.4], method='scipy.root')
x_point_coord = x_point.coords[0]
x_point.plot(ax=ax, marker='x', color="xkcd:crimson")

x_point2= FixedPoint(section)
x_point2.find(1, [0.8,-0.6], method='scipy.root')
x_point_coord2 = x_point2.coords[0]
x_point2.plot(ax=ax, marker='x', color="xkcd:crimson")


# pplot = PoincarePlot.with_segments(section, [top_o_coord, x_point_coord, x_point_coord, x_point_coord2],[30, 10],connected=False)
# pplot.compute(400)


############################################################### Liuqe ############################################################### #########################



data = loadmat('./Script/Snowflake/Liuqe_SF_09.mat', squeeze_me=True, struct_as_record=False)

flux=data['flux']
r = flux.r
z = flux.z
psi = flux.psi.T
psi_r = flux.psi_r
psi_z = flux.psi_z.T
psi_axis = flux.psi_axis



levels1 = np.linspace(-flux.psi_axis, 0, 15)
levels2 = np.linspace(flux.psi_axis, 0, 15)



# Tri pour éviter l'erreur
levels1.sort()
levels2.sort()


# Premier contour
ax.contour(
    psi_r, psi_z, psi,
    levels=levels1,
    colors='cyan', linestyles='-',
    linewidths=0.5

)

# Deuxième contour
ax.contour(
    psi_r, psi_z, psi,
    levels=levels2,
    colors='cyan', linestyles='-',
    linewidths=0.5
)

# Troisième contour : niveau 0, trait plus épais
ax.contour(
    psi_r, psi_z, psi,
    levels=[0],
    colors='k',
    linewidths=2
)




#######################################################################  NO PERT SCRIPT #################################################################### 



#np.save('./Script/Snowflake/Snowflake_NoPert.npy', pplot._hits)
Hits=np.load('./Script/Snowflake/Snowflake_NoPert.npy')



#top fp top manifold

# manifold1 = Manifold(section, x_point, x_point)
# manifold1.compute(
#      eps_s=9e-6, eps_u=8e-6, nint_s=8, nint_u=8, neps_s=80, neps_u=80)
# manifold1.save('./Script/Snowflake/manifolds_NP/mf_1.pkl')


#manifold1  = Manifold.load("./Script/Snowflake/manifolds_NP/mf_1.pkl")

# manifold1.plot(ax=ax, markersize=0, lw=1.5)


# np.save('./Script/Snowflake/manifolds_NP/SF_top_mfs.npy', manifold1.stable)
# np.save('./Script/Snowflake/manifolds_NP/SF_top_mfu.npy', manifold1.unstable)



mfs=np.load('./Script/Snowflake/manifolds_NP/SF_top_mfs.npy')
mfu=np.load('./Script/Snowflake/manifolds_NP/SF_top_mfu.npy')

ax.plot(mfs[:,0],mfs[:,1], color="xkcd:royal blue",markersize=0, lw=1.)
ax.plot(mfu[:,0],mfu[:,1], color="xkcd:magenta",markersize=0, lw=1.)







## top fp bottom manifold
# manifold2 = Manifold(section, x_point, x_point,x_point_coord-top_o_coord, x_point_coord-top_o_coord)
# manifold2.compute(
#    eps_s=9e-6, eps_u=8e-6, nint_s=8, nint_u=8, neps_s=80, neps_u=80)
# manifold2.save('./Script/Snowflake/manifolds_NP/mf_2.pkl')


#manifold2  = Manifold.load("./Script/Snowflake/manifolds_NP/mf_2.pkl")

# manifold2.plot(ax=ax, markersize=0, lw=1.5)

# np.save('./Script/Snowflake/manifolds_NP/SF_top_mfs1.npy', manifold2.stable)
# np.save('./Script/Snowflake/manifolds_NP/SF_top_mfu1.npy', manifold2.unstable)



mfs1=np.load('./Script/Snowflake/manifolds_NP/SF_top_mfs1.npy')
mfu1=np.load('./Script/Snowflake/manifolds_NP/SF_top_mfu1.npy')

ax.plot(mfs1[:,0],mfs1[:,1], color="xkcd:royal blue",markersize=0, lw=1.)
ax.plot(mfu1[:,0],mfu1[:,1], color="xkcd:magenta",markersize=0, lw=1.)








## bottom fp bottom manifolds
# manifold3 = Manifold(section, x_point2, x_point2,x_point_coord2-x_point_coord,x_point_coord2-x_point_coord)
# manifold3.compute(
#    eps_s=9e-6, eps_u=8e-6, nint_s=8, nint_u=8, neps_s=80, neps_u=80)
# manifold3.save('./Script/Snowflake/manifolds_NP/mf_3.pkl')


#manifold3  = Manifold.load("./Script/Snowflake/manifolds_NP/mf_3.pkl")

# manifold3.plot(ax=ax, markersize=0, lw=1.5)

# np.save('./Script/Snowflake/manifolds_NP/SF_down_mfs2.npy', manifold3.stable)
# np.save('./Script/Snowflake/manifolds_NP/SF_down_mfu2.npy', manifold3.unstable)


mfs2=np.load('./Script/Snowflake/manifolds_NP/SF_down_mfs2.npy')
mfu2=np.load('./Script/Snowflake/manifolds_NP/SF_down_mfu2.npy')

ax.plot(mfs2[:,0],mfs2[:,1], color="xkcd:royal blue",markersize=0, lw=1.)
ax.plot(mfu2[:,0],mfu2[:,1], color="xkcd:magenta",markersize=0, lw=1.)





## bottom fp top manifolds
# manifold4 = Manifold(section, x_point2, x_point2,-x_point_coord2+x_point_coord,-x_point_coord2+x_point_coord)
# manifold4.compute(
#    eps_s=9e-6, eps_u=8e-6, nint_s=8, nint_u=12, neps_s=80, neps_u=80)
# manifold4.save('./Script/Snowflake/manifolds_NP/mf_4.pkl')

#manifold4  = Manifold.load("./Script/Snowflake/manifolds_NP/mf_4.pkl")


# manifold4.plot(ax=ax, markersize=0, lw=1.5)

# np.save('./Script/Snowflake/manifolds_NP/SF_down_mfs1.npy', manifold4.stable)
# np.save('./Script/Snowflake/manifolds_NP/SF_down_mfu1.npy', manifold4.unstable)


mfs3=np.load('./Script/Snowflake/manifolds_NP/SF_down_mfs1.npy')
mfu3=np.load('./Script/Snowflake/manifolds_NP/SF_down_mfu1.npy')

ax.plot(mfs3[:,0],mfs3[:,1], color="xkcd:royal blue",markersize=0, lw=1.)
ax.plot(mfu3[:,0],mfu3[:,1], color="xkcd:magenta",markersize=0, lw=1.)



ratio=2.8158

ax.scatter(Hits[:,:, 0], Hits[:,:, 1], color="xkcd:dark grey", s=1., linewidths=0)
ax.set_xlim(0.62, 1.15)
ax.set_ylim(-0.75, 0.75)
ax.set_xlabel(r"$R[m]$")
ax.set_ylabel(r"$Z[m]$")
ax.set_aspect('equal') 


plt.savefig('./Script/Snowflake/Snowflake_NoPert_Liuqe.png', bbox_inches='tight', dpi=720)
plt.show()




















######################################################################################### PERTURBED SCRIPT ######################################################################################




# #np.save('./Script/Snowflake/Snowflake_NoPert.npy', pplot._hits)
# Hits=np.load('./Script/Snowflake/Snowflake_Pert.npy')



# #top fp top manifold
# # manifold1 = Manifold(section, x_point, x_point)
# # manifold1.compute(
# #      eps_s=9e-6, eps_u=8e-6, nint_s=10, nint_u=10, neps_s=160, neps_u=160)
# manifold1.save('./Script/Snowflake/manifolds_P/mf_1.pkl')


# manifold1  = Manifold.load("./Script/Snowflake/manifolds_P/mf_1.pkl")

# # manifold1.plot(stepsize_limit=0.2,ax=ax, markersize=0, lw=1.5)



# #np.save('./Script/Snowflake/manifolds_P/SF_top_mfs.npy', manifold1.stable)
# #np.save('./Script/Snowflake/manifolds_P/SF_top_mfu.npy', manifold1.unstable)



# mfs=np.load('./Script/Snowflake/manifolds_P/SF_top_mfs.npy')
# mfu=np.load('./Script/Snowflake/manifolds_P/SF_top_mfu.npy')

# ax.plot(mfs[:,0],mfs[:,1], color="xkcd:royal blue",markersize=0, lw=1.)
# ax.plot(mfu[:,0],mfu[:,1], color="xkcd:magenta",markersize=0, lw=1.)








# ##top fp bottom manifold


# # manifold2 = Manifold(section, x_point, x_point,x_point_coord-top_o_coord, x_point_coord-top_o_coord)
# # manifold2.compute(
# #    eps_s=9e-6, eps_u=8e-6, nint_s=8, nint_u=8, neps_s=80, neps_u=80)
# manifold2.save('./Script/Snowflake/manifolds_P/mf_2.pkl')


# manifold2 = Manifold.load("./Script/Snowflake/manifolds_P/mf_2.pkl")
# # manifold2.plot(ax=ax, markersize=0, lw=1.5,colors=["cyan", "tomato"])

# # np.save('./Script/Snowflake/manifolds_P/SF_top_mfs1.npy', manifold2.stable)
# # np.save('./Script/Snowflake/manifolds_P/SF_top_mfu1.npy', manifold2.unstable)




# mfs1=np.load('./Script/Snowflake/manifolds_P/SF_top_mfs1.npy')
# mfu1=np.load('./Script/Snowflake/manifolds_P/SF_top_mfu1.npy')

# ax.plot(mfs1[:,0],mfs1[:,1], color="cyan",markersize=0, lw=1.)
# ax.plot(mfu1[:,0],mfu1[:,1], color="fuchsia",markersize=0, lw=1.)








# ## bottom fp bottom manifolds


# # manifold3 = Manifold(section, x_point2, x_point2,x_point_coord2-x_point_coord,x_point_coord2-x_point_coord)
# # manifold3.compute(
# #    eps_s=9e-6, eps_u=8e-6, nint_s=8, nint_u=8, neps_s=80, neps_u=80)
# manifold3.save('./Script/Snowflake/manifolds_P/mf_3.pkl')


# manifold3 = Manifold.load("./Script/Snowflake/manifolds_P/mf_3.pkl")
# # manifold3.plot(ax=ax, markersize=0, lw=1.5,colors=["cornflowerblue", "lightcoral"])


# # np.save('./Script/Snowflake/manifolds_P/SF_down_mfs2.npy', manifold3.stable)
# # np.save('./Script/Snowflake/manifolds_P/SF_down_mfu2.npy', manifold3.unstable)




# mfs2=np.load('./Script/Snowflake/manifolds_P/SF_down_mfs2.npy')
# mfu2=np.load('./Script/Snowflake/manifolds_P/SF_down_mfu2.npy')

# ax.plot(mfs2[:,0],mfs2[:,1], color="cornflowerblue",markersize=0, lw=1.)
# ax.plot(mfu2[:,0],mfu2[:,1], color="lightcoral",markersize=0, lw=1.)








# ## bottom fp top manifolds


# # manifold4 = Manifold(section, x_point2, x_point2,-x_point_coord2+x_point_coord,-x_point_coord2+x_point_coord)
# # manifold4.compute(
# #    eps_s=9e-6, eps_u=8e-6, nint_s=8, nint_u=10, neps_s=80, neps_u=80)
# manifold4.save('./Script/Snowflake/manifolds_P/mf_4.pkl')


# manifold4 = Manifold.load("./Script/Snowflake/manifolds_P/mf_4.pkl")
# # manifold4.plot(ax=ax, markersize=0, lw=1.5,colors=["cornflowerblue", "lightcoral"])

# # np.save('./Script/Snowflake/manifolds_P/SF_down_mfs1.npy', manifold4.stable)
# # np.save('./Script/Snowflake/manifolds_P/SF_down_mfu1.npy', manifold4.unstable)



# mfs3=np.load('./Script/Snowflake/manifolds_P/SF_down_mfs1.npy')
# mfu3=np.load('./Script/Snowflake/manifolds_P/SF_down_mfu1.npy')

# ax.plot(mfs3[:,0],mfs3[:,1], color="cornflowerblue",markersize=0, lw=1.)
# ax.plot(mfu3[:,0],mfu3[:,1], color="lightcoral",markersize=0, lw=1.)







# ratio=2.8158

# ax.scatter(Hits[:,:, 0], Hits[:,:, 1], color="xkcd:dark grey", s=1.4, linewidths=0)
# ax.set_xlim(0.62, 1.15)
# ax.set_ylim(-0.75, 0.75)
# ax.set_xlabel(r"$R[m]$")
# ax.set_ylabel(r"$Z[m]$")
# ax.set_aspect('equal') 



# plt.savefig('./Script/Snowflake/Snowflake_Pert_mf.png', bbox_inches='tight', dpi=720)
# plt.show()






