import os, sys
proj_root = os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..','..','..','..'))  # remonte jusqu'à "Travail master"
if proj_root not in sys.path:
    sys.path.insert(0, proj_root)
from script.function.field_utils import plot_tomographic as tomo
from script.function.field_utils import plot_LIUQE as LIUQE

from pyoculus.fields import AxisymmetricCylindricalGridField
from pyoculus.maps import CylindricalBfieldSection
from pyoculus.solvers import FixedPoint, Manifold, PoincarePlot
import numpy as np
from matplotlib import pyplot as plt
import logging
from scipy.io import loadmat
from matplotlib.tri import Triangulation

fig, ax = plt.subplots(1, 1, figsize=(6, 6.6))

plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "serif",
        "font.size": 9,
    }
)

logging.basicConfig(level=logging.DEBUG)

repository_path = './script/PC_analysis/flux/mat_files/'


file_manifolds = './script/PC_analysis/flux/mat_files/manifold/'


# ...existing code...
# boucle de 1 à 18 (inclus)
for i in range(1, 19):
   for j in range(1,2): 
    print(f"[{i}/18] traitement du fichier mat {j/2}")
    # exemple de pattern : "JF_77021_120_BO78_OS_01.mat", adapte si nécessaire
    pert_mat_file = f"JF_77021_120_BO78_OS_{i:02d}.mat"

    full_path = os.path.join(repository_path, pert_mat_file)

   

    # charge le champ et poursuit le traitement habituel
    JFField = AxisymmetricCylindricalGridField.from_matlab_file(full_path, with_perturbation=True)

    # --- place ici le code que tu veux exécuter pour chaque JFField ---
    # ex: creation de la section, calculs, tracés, sauvegardes...
    section = CylindricalBfieldSection(JFField, phi0=1.9, R0=0.88, Z0=0)
    # ...existing processing code...

    top_o = FixedPoint(section)
    top_o.find(1, [0.9,0.18], method='scipy.root')
    top_o_coord = top_o.coords[0] 
    x_point1 = FixedPoint(section)
    x_point1.find(1, [0.8,-0.27], method='scipy.root')
    x_point1_coord = x_point1.coords[0]


   manifold_1T = Manifold(section, x_point1, x_point1,-x_point1_coord+top_o_coord, -x_point1_coord+top_o_coord)
   manifold_1T.compute(eps_s=9e-6, eps_u=8e-6, nint_s=14, nint_u=14, neps_s=80, neps_u=240) #8 14 240
   manifold_1T.find_clinics(first_guess_eps_s=9e-6, first_guess_eps_u=8e-6)
   manifold_1T.compute_turnstile_areas()
   manifold_1T.save(f"{file_manifolds}coil_{i}_time_{j}.pkl")

   

   

   

# ...existing code...





