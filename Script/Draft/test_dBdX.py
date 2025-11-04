from pyoculus.fields import AxisymmetricCylindricalGridField
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.fields import AnalyticCylindricalBfield
from pyoculus.solvers import PoincarePlot, FixedPoint
import numpy as np
from matplotlib import pyplot as plt
import logging
logging.basicConfig(level=logging.DEBUG)

plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.size": 10,
    }
)
logging.basicConfig(level=logging.DEBUG)

DoubletField = AxisymmetricCylindricalGridField.from_matlab_file('./Script/Doublet/struct_for_chris.mat', with_perturbation=True)

SHEAR =  0.4
SF = 0.8875

perturbedfield = AnalyticCylindricalBfield.with_new_axis(
    R=6, Z=0, sf=SF, shear=SHEAR, perturbations_args=[]
)
perturbedmap = CylindricalBfieldSection.without_axis(perturbedfield, guess=[6., 0])
axiscoords = np.array((perturbedmap.R0, perturbedmap.Z0))






#R=np.linspace(0.5,1.5,10)
#Z=np.linspace(-1,+1,15) 
#B_R = 0
#B_Z = 0
#B_phi = (1/R*np.ones((15,10))).T
#F_psi =0


#print(DoubletField.pertfield.dBdX([0.8,0.0 , 0.4]))
print(perturbedfield.dBdX([6.3,0.0 , 0.5]))
print((perturbedfield.B([6.3001,0.0 , 0.5])-perturbedfield.B([6.3,0.0 , 0.5]))/0.0001)
print((perturbedfield.B([6.3,0.0 , 0.5001])-perturbedfield.B([6.3,0.0 , 0.5]))/0.0001)


print(DoubletField.dBdX([0.8, 0.0, 0.4]))

print((DoubletField.B([0.80000000001, 0.0, 0.4])-DoubletField.B([0.8, 0.0, 0.4]))/0.00000000001)
print((DoubletField.B([0.8, 0.00000001, 0.4])-DoubletField.B([0.8, 0.0, 0.4]))/0.00000001)
print((DoubletField.B([0.8, 0.0, 0.40000001])-DoubletField.B([0.8, 0.0, 0.4]))/0.00000001)


#print(DoubletField.A([0.8, 0.0, 0.4]))
#print(DoubletField.A_unperturbed([0.8, 0.0, 0.4]))
#print(DoubletField.pertfield.A([0.8, 0.0, 0.4]))
