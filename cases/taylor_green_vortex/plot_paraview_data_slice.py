import sys
from sys import platform
if platform == "linux" or platform == "linux2":
    # linux
    filesystem="/media/julien/Samsung_T5/"
elif platform == "darwin":
    # OS X
    filesystem="/Volumes/Samsung_T5/"

import os;CURRENT_PATH = os.path.split(os.path.realpath(__file__))[0]+"/";
sys.path.append(CURRENT_PATH+"../../submodules/quickplotlib/lib");
import quickplotlib as qp
import numpy as np
import pandas as pd

x_cut=-3.14159265358979

# # path="/Users/Julien/NarvalFiles/2022-11-09_96dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs096_p5_procs512/"
# path="/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64_filter36timeslarger/"


path = filesystem+"NarvalFiles/2023_JCP/verification/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p3_procs1024/solution_files/"
filename = path+"vorticity_plane_from_paraview_test.txt"

df = pd.read_csv(filename,sep=",",header=0)
scalar = df["vorticity_magnitude"].to_numpy()
x = df["x"].to_numpy()
y = df["y"].to_numpy()
z = df["z"].to_numpy()
np.savetxt("dummy_y.txt",y)
np.savetxt("dummy_z.txt",z)

# ----------------------------------------

# plotting code
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import rc as matplotlibrc
matplotlibrc('text.latex', preamble='\\usepackage{color}')
matplotlibrc('text', usetex=True)
matplotlibrc('font', family='serif')

# Figure parameters
fig_directory = "figures"
figure_filetype = "pdf"
# Font sizes
axisTitle_FontSize = 16
axisTickLabel_FontSize = 14
legend_fontSize = 16

y_label = "$z$"
x_label = "$y$"
z_label = "$|\\Omega|$"
fig, ax = plt.subplots(figsize=(6,6))
# fig = plt.figure(figsize=(6,6))
# plt.xlim([np.amin(X),np.amax(X)])
# plt.ylim([np.amin(Y),np.amax(Y)])
# plt.xlim([0.0,1.4])
# plt.ylim([1.6,2.6])
ax.set_xlabel(x_label,fontsize=axisTitle_FontSize)
ax.set_ylabel(y_label,rotation=90,fontsize=axisTitle_FontSize)
plt.setp(ax.get_xticklabels(),fontsize=axisTickLabel_FontSize); plt.setp(ax.get_yticklabels(),fontsize=axisTickLabel_FontSize);

X = y
Y = z
Z = scalar
plt.tricontourf(X, Y, Z, np.linspace(np.amin(Z),np.amax(Z),100),cmap='rainbow')
cbar = fig.colorbar(cs,ticks=np.arange(np.amin(Z),np.amax(Z),50.0))
# cbar = fig.colorbar(cs,ticks=np.linspace(np.amin(Z),np.amax(Z),8))
cbar.set_label(z_label,fontsize=axisTickLabel_FontSize)
# cbar.ax.tick_params(fontsize=axisTickLabel_FontSize)
plt.setp(cbar.ax.get_yticklabels(),fontsize=axisTickLabel_FontSize)

# Fix for the white lines between contour levels
for c in cs.collections:
    c.set_edgecolor("face")


plt.show()

exit()
# # path="/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs096_p5_procs512/"
# filename=path+"flow_field_files/velocity_vorticity-0_reordered.dat"
# # date_for_runs="flow_field_files/velocity_vorticity-0_reordered.dat"
# figure_filename = "vorticity_mag"
# vorticity_magnitude_field = np.loadtxt(filename,skiprows=1,usecols=(0,1,2,6),dtype=np.float64)
# # velocity_field = np.loadtxt("dofs128_p3_velocity/velocity_equidistant_nodes.fld",skiprows=0,usecols=(0,1,2,3,4,5),dtype=np.float64)

# row_indices_of_plane = np.where(np.isclose(vorticity_magnitude_field[:,0],x_cut,atol=1e-10))[0]
# npoints = row_indices_of_plane.size
# y = np.zeros(npoints,dtype=np.float64)
# z = np.zeros(npoints,dtype=np.float64)
# vorticity_mag = np.zeros(npoints,dtype=np.float64)
# for i in range(npoints):
#     index = row_indices_of_plane[i]
#     y[i] = vorticity_magnitude_field[index,1]
#     z[i] = vorticity_magnitude_field[index,2]
#     vorticity_mag[i] = vorticity_magnitude_field[index,3]

# ny = int(np.sqrt(npoints))
# nz = int(np.sqrt(npoints))
# y_mesh = np.zeros((nz,ny),dtype=np.float64)
# z_mesh = np.zeros((nz,ny),dtype=np.float64)
# vorticity_mag_mesh = np.zeros((nz,ny),dtype=np.float64)

# index_global = 0
# for j in range(0,nz):
#     for i in range(0,ny):
#         y_mesh[j,i] = y[index_global]
#         z_mesh[j,i] = z[index_global]
#         vorticity_mag_mesh[j,i] = vorticity_mag[index_global]
#         index_global += 1

# X = y_mesh
# Y = z_mesh
# Z = vorticity_mag_mesh

# # plotting code
# import matplotlib.pyplot as plt
# from matplotlib.lines import Line2D
# from matplotlib import rc as matplotlibrc
# matplotlibrc('text.latex', preamble='\\usepackage{color}')
# matplotlibrc('text', usetex=True)
# matplotlibrc('font', family='serif')

# # Figure parameters
# fig_directory = "figures"
# figure_filetype = "pdf"
# # Font sizes
# axisTitle_FontSize = 16
# axisTickLabel_FontSize = 14
# legend_fontSize = 16

# y_label = "$z$"
# x_label = "$y$"
# z_label = "$|\\Omega|$"

# print('Plotting: ' + figure_filename)
# fig, ax = plt.subplots(figsize=(6,6))
# # fig = plt.figure(figsize=(6,6))
# plt.xlim([np.amin(X),np.amax(X)])
# plt.ylim([np.amin(Y),np.amax(Y)])
# ax.set_xlabel(x_label,fontsize=axisTitle_FontSize)
# ax.set_ylabel(y_label,rotation=90,fontsize=axisTitle_FontSize)
# plt.setp(ax.get_xticklabels(),fontsize=axisTickLabel_FontSize); plt.setp(ax.get_yticklabels(),fontsize=axisTickLabel_FontSize);

# cs = plt.contourf(X,Y,np.transpose(Z), np.linspace(np.amin(Z),np.amax(Z),100),cmap='rainbow')

# cbar = fig.colorbar(cs,ticks=np.arange(np.amin(Z),np.amax(Z),50.0))
# # cbar = fig.colorbar(cs,ticks=np.linspace(np.amin(Z),np.amax(Z),8))
# cbar.set_label(z_label,fontsize=axisTickLabel_FontSize)
# # cbar.ax.tick_params(fontsize=axisTickLabel_FontSize)
# plt.setp(cbar.ax.get_yticklabels(),fontsize=axisTickLabel_FontSize)

# # Fix for the white lines between contour levels
# for c in cs.collections:
#     c.set_edgecolor("face")

# # leg = plt.legend(loc='best', ncol=1, shadow=False, fancybox=True, fontsize=legend_fontSize, framealpha=1.0,edgecolor='inherit')

# plt.tight_layout()
# print('... Saving figure ...')
# plt.savefig(fig_directory+"/"+figure_filename+'.'+figure_filetype,format=figure_filetype,dpi=500)
# plt.close()

# # # qp.plotfield(xdata=x,ydata=y,udata=u,vdata=v,ylabel="$y$",xlabel="$x$",
# # #     # title_label="DHIT Initial Velocity Field at $z=\\pi$\n $24^{3}$ DOF (P5, $N_{el}=4^{3}$)",
# # #     fig_directory=".",figure_filename="velocity_field_48_init",
# # #     xlimits=[0.0,2.0*np.pi],ylimits=[0.0,2.0*np.pi])
