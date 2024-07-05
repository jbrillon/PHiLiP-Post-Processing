#-----------------------------------------------------
# Import public libraries
import numpy as np # NumPy: contains basic numerical routines
from scipy.optimize import curve_fit
#-----------------------------------------------------
# Import personal libraries
# from finite_difference_library import first_derivative, fd_non_uniform_grid
import os;CURRENT_PATH = os.path.split(os.path.realpath(__file__))[0]+"/";
import sys
sys.path.append(CURRENT_PATH+"../../submodules/quickplotlib/lib"); import quickplotlib as qp
#-----------------------------------------------------
from sys import platform
if platform == "linux" or platform == "linux2":
    # linux
    filesystem="/media/julien/Samsung_T5/"
elif platform == "darwin":
    # OS X
    filesystem="/Volumes/Samsung_T5/"
#-----------------------------------------------------
#=====================================================
# Helper functions
#=====================================================
#-----------------------------------------------------
def u_plus_y_plus(y_plus):
    kappa = 0.38
    C = 4.1
    u_plus = (1.0/kappa)*np.log(1.0+kappa*y_plus) + (C-(1.0/kappa)*np.log(kappa))*(1.0 -np.exp(-1.0*y_plus/11.0) +(-1.0*y_plus/11.0)*np.exp(-y_plus/3.0))
    return u_plus*y_plus
#-----------------------------------------------------
def func(x, a, b):
    # return a * np.exp(-b * x) + c
    return a *x + b
#-----------------------------------------------------
y_plus_store = np.linspace(0,5000,2500)
u_plus_y_plus_store = u_plus_y_plus(y_plus_store)
x_store=[]
y_store=[]
labels_store=[]
x_store.append(y_plus_store)
y_store.append(u_plus_y_plus_store)
labels_store.append("Analytical")

y_plus_sample_points = np.array([0.0,3.0,5.0,8.0,10.0,20.0,35.0,50.0,75.0,100.0,125.0,150.0,200.0,250.0,300.0,350.0,400.0,500.0,575.0,650.0,725.0,800.0,900.0,1000.0,1100.0,1200.0,1300.0,1400.0,1500.0,1600.0,1800.0,2000.0,2500.0,3000.0,3500.0,4000.0,4500.0,5000.0])
u_plus_y_plus_sample = u_plus_y_plus(y_plus_sample_points)
# print(y_plus_sample_points)
# print(u_plus_y_plus_sample)
print("Number of sample points:")
print(np.size(y_plus_sample_points))
print("(y+) sample points:")
for i in range(np.size(y_plus_sample_points)):
    print("%.1f, " % y_plus_sample_points[i],end="")
print(" ")
print("(u+)*(y+) sample points:")
for i in range(np.size(y_plus_sample_points)):
    print("%1.13e, " % u_plus_y_plus_sample[i],end="")
print(" ")

# re_inf = 129245.0
# val_from_philip = np.interp(68156.297044523541,y_plus_sample_points,u_plus_y_plus_sample)
# print(val_from_philip)
# # exit()
# val = u_plus_y_plus(1300)
# print("%1.12e" % val)
# exit()

x_store.append(y_plus_sample_points)
y_store.append(np.interp(y_plus_sample_points,y_plus_sample_points,u_plus_y_plus_sample))
labels_store.append("Linear interpolation")

# xdata=1.0*np.linspace(0,10,2500)
# ydata=1.0*u_plus_y_plus(xdata)
# popt, pcov = curve_fit(func, xdata, ydata)
# x_store.append(xdata)
# y_store.append(func(xdata,*popt))
# labels_store.append("Fit 1")

# xdata=1.0*np.linspace(10,100,2500)
# ydata=1.0*u_plus_y_plus(xdata)
# popt, pcov = curve_fit(func, xdata, ydata)
# x_store.append(xdata)
# y_store.append(func(xdata,*popt))
# labels_store.append("Fit 2")

# xdata=1.0*np.linspace(100,1000,2500)
# ydata=1.0*u_plus_y_plus(xdata)
# popt, pcov = curve_fit(func, xdata, ydata)
# x_store.append(xdata)
# y_store.append(func(xdata,*popt))
# labels_store.append("Fit 3")

# print(np.linalg.norm(func(xdata,*popt)-u_plus_y_plus_store)/np.length)

qp.plotfxn(x_store,y_store,
        figure_filename="wall_model_function",
        figure_size=(8,6),
        legend_labels_tex=labels_store,
        figure_filetype="pdf",
        title_label="Reichardt's Law of the Wall Multiplied by $y^{+}$",
        xlabel="$y^{+}$",
        ylabel="$u^{+}_{||}y^{+}$",
        which_lines_dashed=[1],
        transparent_legend=False,
        legend_border_on=True,
        grid_lines_on=True,
        log_axes="y",
        legend_location="best",
        which_lines_markers=[1])