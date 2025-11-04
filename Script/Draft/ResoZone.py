
from pyoculus.fields import AnalyticCylindricalBfield
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
from matplotlib import pyplot as plt
import numpy as np
import logging
from matplotlib.path import Path


def ComputeShoelace(xy): 
    """
    Shoelace formula for calculating the area of a polygon. 
    xy: list of tuples or arrays of x,y coordinates of the polygon vertices. 
    """

    x = np.array([p[0] for p in xy])
    y = np.array([p[1] for p in xy])
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))



plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.size": 10,
    }
)
logging.basicConfig(level=logging.DEBUG)

POINCARE_ITS = 500
SHEAR =  0.4
SF = 0.8875

maxwellboltzmann = {
    "m": 3,
    "n": 2,
    "d": 1.75/np.sqrt(2),
    "type": "maxwell-boltzmann",
    "amplitude": 0.1,
}

separatrix = {"type": "circular-current-loop", "amplitude": -10, "R": 6, "Z": -5.5}

fig, ax = plt.subplots(1, 1, figsize=(6, 6.6))



perturbedfield = AnalyticCylindricalBfield.with_new_axis(
    R=6, Z=0, sf=SF, shear=SHEAR, perturbations_args=[maxwellboltzmann]
)
perturbedmap = CylindricalBfieldSection.without_axis(perturbedfield, guess=[6., 0])
axiscoords = np.array((perturbedmap.R0, perturbedmap.Z0))


#opoint1 = FixedPoint(perturbedmap)
#opoint1.find(t=1, guess=[6.0, 0.0], method='scipy.root')
#opointcoords1 = opoint1.coords[0]



opointT = FixedPoint(perturbedmap)
#opointT.find(t=3, guess=[7.6, 0.], method='scipy.root')
opointT.find(t=3, guess=[4.7, 1.3], method='scipy.root')
#opointT.find(t=3, guess=[4.7, -1.2], method='scipy.root')
opointTcoords = opointT.coords[0]


#opointR = FixedPoint(perturbedmap)
#opointR.find(t=3, guess=[7.6, 0.], method='scipy.root')
#opointRcoords = opointR.coords[0]

#opointB = FixedPoint(perturbedmap)
#opointB.find(t=3, guess=[4.7, -1.2], method='scipy.root')
#opointBcoords = opointB.coords[0]



xpointL = FixedPoint(perturbedmap)
xpointL.find(t=3, guess=[4.5, 0.0], method='scipy.root')
#xpointL.find(t=3, guess=[6.4, -1.6], method='scipy.root')
xpointLcoords = xpointL.coords[0]

xpointT = FixedPoint(perturbedmap)
#xpointT.find(t=3, guess=[6.4, 1.6], method='scipy.root')
xpointT.find(t=3, guess=[6.4, -1.6], method='scipy.root')
xpointTcoords = xpointT.coords[0]

#xpointB = FixedPoint(perturbedmap)
#xpointB.find(t=3, guess=[6.4, -1.6], method='scipy.root')
#xpointBcoords = xpointB.coords[0]


opointT.plot(ax=ax, marker='o', color="xkcd:green")
#opointB.plot(ax=ax, marker='o', color="xkcd:green")
#opointR.plot(ax=ax, marker='o', color="xkcd:green")
xpointL.plot(ax=ax, marker='x', color="xkcd:green")
xpointT.plot(ax=ax, marker='x', color="xkcd:green")
#xpointB.plot(ax=ax, marker='x', color="xkcd:green")





#pplot = PoincarePlot.with_linspace(
#    perturbedmap,
#    opointTcoords,
#    xpointTcoords,
#    30
#   )

#pplot = PoincarePlot.with_segments(
#    perturbedmap,
#    [opointTcoords, xpointTcoords, opointTcoords, opointTcoords + np.array([-0.4, 0.25])],
#    neps=[30, 30],
#    connected=False)





#pplot.compute(POINCARE_ITS)
#pplot.plot(ax=ax, color="xkcd:dark grey", s=1.4, linewidths=0)
#np.save('Reso_Hits', pplot._hits)








#Hits=np.load('Reso_Hits.npy')
#ax.scatter(Hits[:,:, 0], Hits[:,:, 1], color="xkcd:dark grey", s=1.4, linewidths=0)



#ax.set_xlabel(r"$R[m]$")
#ax.set_ylabel(r"$Z[m]$")





ResoManifoldDown = Manifold(perturbedmap, xpointL, xpointT)
#ResoManifoldDown.compute(
 # eps_s=9e-6, eps_u=8e-6, nint_s=8, nint_u=8, neps_s=80, neps_u=80
 #  )
#ResoManifoldDown.plot(ax=ax, markersize=0, lw=1.5)

ResoManifoldDown.find_clinics(first_guess_eps_s=9e-6, first_guess_eps_u=8e-6)
ResoManifoldDown.plot_clinics(ax=ax)
traj1=ResoManifoldDown.clinics[0].trajectory
traj1bis=ResoManifoldDown.clinics[1].trajectory


if len(traj1)!=len(traj1bis):
        minlen=min(len(traj1), len(traj1bis))
        traj1=traj1[:minlen]
        traj1bis=traj1bis[:minlen]



ResoManifoldUp = Manifold(perturbedmap, xpointT, xpointL)
#ResoManifoldUp.compute(
 #   eps_s=9e-6, eps_u=8e-6, nint_s=8, nint_u=8, neps_s=80, neps_u=80
#)
ResoManifoldUp.plot(ax=ax, markersize=0, lw=1.5)#

ResoManifoldUp.find_clinics(first_guess_eps_s=9e-6, first_guess_eps_u=8e-6)
ResoManifoldUp.plot_clinics(ax=ax)
traj2=ResoManifoldUp.clinics[0].trajectory
traj2bis=ResoManifoldUp.clinics[1].trajectory

if len(traj2)!=len(traj2bis):
        minlen=min(len(traj2), len(traj2bis))
        traj2=traj2[:minlen]
        traj2bis=traj2bis[:minlen]





int1= np.sum(ResoManifoldDown._compute_lagrangian_in_sections2(traj1))#+np.sum(ResoManifoldDown._compute_lagrangian_in_sections2(traj1bis))
int1bis= np.sum(ResoManifoldDown._compute_lagrangian_in_sections2(traj1bis))
int2=(len(traj1)+len(traj1bis))*ResoManifoldDown._map.lagrangian(xpointLcoords,1)
#int2=len(traj1)*ResoManifoldDown._map.lagrangian(xpointLcoords,1)

int3= np.sum(ResoManifoldUp._compute_lagrangian_in_sections2(traj2))#+np.sum(ResoManifoldUp._compute_lagrangian_in_sections2(traj2bis))
int3bis= np.sum(ResoManifoldUp._compute_lagrangian_in_sections2(traj2bis))
int4=(len(traj2)+len(traj2bis))*ResoManifoldUp._map.lagrangian(xpointTcoords,1)
#int4=(len(traj2)*ResoManifoldUp._map.lagrangian(xpointTcoords,1))





#plt.savefig('Max3mReso.png', bbox_inches='tight', dpi=720)

#print(ResoManifoldUp._map.lagrangian(traj1[0],1),ResoManifoldUp._map.lagrangian(traj1[-1],1),ResoManifoldUp._map.lagrangian(traj2[0],1),ResoManifoldUp._map.lagrangian(traj2[-1],1))
#FluxCheckSL=ComputeShoelace(contour_points)*perturbedfield.B(opointTcoords)
#print(FluxCheckSL)
#FluxCheckADL=ResoManifoldUp._AdL_integral_points(contour_points,dl=None, is_closed=True)
#print(FluxCheckSL, FluxCheckADL,ResoManifoldDown._map.lagrangian(xpointLcoords,1),ResoManifoldUp._map.lagrangian(xpointTcoords,1))

N=2*len(traj1)
points = np.zeros((N+1, 2))
lagrangian = np.zeros(N)
points[0] = xpointLcoords
for i in range(1, N):
    points[i] = perturbedmap.f(1,points[i-1])
for i in range(0, N):    
    lagrangian[i] = ResoManifoldUp._map.lagrangian(points[i],1)

#int2=np.sum(lagrangian)


N=2*len(traj2)
points1 = np.zeros((N+1, 2))
lagrangian1 = np.zeros(N)
points1[0] = xpointTcoords
for i in range(1, N):
    points1[i] = perturbedmap.f(1,points1[i-1])
for i in range(0, N):    
    lagrangian1[i] = ResoManifoldUp._map.lagrangian(points1[i],1)

#int4=np.sum(lagrangian1)

Areso= (int1+int1bis+int3+int3bis-int2-int4)

#print(int1,int1bis,int2,int3,int3bis,int4)
#print(Areso)
#print(len(traj1),len(traj1bis),len(lagrangian), len(traj2),len(traj2bis),len(lagrangian1))


#eps=9e-6
#rEps = x_point_coord + eps * manifold1.vector_s
#feps = manifold1._map.f(-1,rEps)
#vec=(feps - rEps)
#vec= vec/np.linalg.norm(vec)
#vec2=manifold1.vector_s/np.linalg.norm(manifold1.vector_s)

#ax.plot([x_point_coord[0], x_point_coord[0]+vec2[0]],[x_point_coord[1],x_point_coord[1]+vec2[1]], 'b-')
#ax.plot([x_point_coord[0],x_point_coord[0]+vec[0]],[x_point_coord[1],x_point_coord[1]+vec[1]], 'r-')





#print(points)
#print(lagrangian)
#print(points1)
#print(lagrangian1)
#print(np.sum(lagrangian), np.sum(lagrangian1))

print(traj2[6],perturbedmap.f(1,traj2[6]),perturbedmap.f(2,traj2[6]),perturbedmap.f(3,traj2[6]))
