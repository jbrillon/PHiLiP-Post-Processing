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

def write_coordinates_only(number_of_points,filename,x,y):
    file = open(filename,"w")
    for i in range(0,number_of_points):
        wstr = "%1.6e %1.6e\n" % (x[i],y[i]) # 6 decimals to be consistent with input csv files
        file.write(wstr)
    file.close()
    return

x_store = []
y_store = []

chord=10.0
number_of_points=61
# 0 degree AoA
x,y = np.loadtxt("sd7003-il_alpha_0_deg.csv",skiprows=9,max_rows=number_of_points,unpack=True,delimiter=",")
x_store.append(x/chord)
y_store.append(y/chord)
# write coordinates normalized by chord length
write_coordinates_only(number_of_points,"sd7003-il_alpha_0_deg_surface_only_normalized.txt",x/chord,y/chord)
# 8 degree AoA
x,y = np.loadtxt("sd7003-il_alpha_8_deg.csv",skiprows=9,max_rows=number_of_points,unpack=True,delimiter=",")
x_store.append(x/chord)
y_store.append(y/chord)
# write coordinates normalized by chord length
write_coordinates_only(number_of_points,"sd7003-il_alpha_8_deg_surface_only_normalized.txt",x/chord,y/chord)

qp.plotfxn(x_store,y_store,
    figure_filename="sd7003coordinates",
    figure_size=(9,3),
    legend_labels_tex=["$\\alpha=0^\\circ$","$\\alpha=8^\\circ$"],
    figure_filetype="png",
    title_label="SD7003 Airfoil",
    xlabel="$x$",
    ylabel="$y$")