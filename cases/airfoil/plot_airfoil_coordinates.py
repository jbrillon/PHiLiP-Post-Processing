import numpy as np


#-----------------------------------------------------
# Import public libraries
import numpy as np # NumPy: contains basic numerical routines
import scipy # SciPy: contains additional numerical routines to numpy
#-----------------------------------------------------
# Import personal libraries
# from finite_difference_library import first_derivative, fd_non_uniform_grid
import os;CURRENT_PATH = os.path.split(os.path.realpath(__file__))[0]+"/";
import sys
sys.path.append(CURRENT_PATH+"../../submodules/quickplotlib/lib"); import quickplotlib as qp

x_store = []
y_store = []

chord=10.0
x,y = np.loadtxt("sd7003-il_alpha_0_deg.csv",skiprows=9,max_rows=61,unpack=True,delimiter=",")
x_store.append(x/chord)
y_store.append(y/chord)

x,y = np.loadtxt("sd7003-il_alpha_8_deg.csv",skiprows=9,max_rows=61,unpack=True,delimiter=",")
x_store.append(x/chord)
y_store.append(y/chord)

qp.plotfxn(x_store,y_store,
    figure_filename="sd7003coordinates",
    figure_size=(9,3),
    legend_labels_tex=["$\\alpha=0^\\circ$","$\\alpha=8^\\circ$"],
    figure_filetype="png",
    title_label="SD7003 Airfoil",
    xlabel="$x$",
    ylabel="$y$")