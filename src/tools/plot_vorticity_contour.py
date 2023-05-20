import numpy as np
# import sys; sys.path.append("submodules/quickplotlib/lib"); import quickplotlib as qp
import sys; sys.path.append("/Users/Julien/Python/quickplotlib/lib"); import quickplotlib as qp

# x_cut=-3.14159265358979

# # path="/Users/Julien/NarvalFiles/2022-11-09_96dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs096_p5_procs512/"
# path="/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64_filter36timeslarger/"
# # path="/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs096_p5_procs512/"
# filename=path+"flow_field_files/velocity_vorticity-0_reordered.dat"

nDOFs_per_dim=256
nDOFs_per_plane=nDOFs_per_dim*nDOFs_per_dim

filename="vorticity_x_plane_DNS-1.dat"
data=np.loadtxt(filename,dtype=np.float64)

ny = nDOFs_per_dim
nz = nDOFs_per_dim
y_mesh = np.zeros((ny,nz),dtype=np.float64)
z_mesh = np.zeros((ny,nz),dtype=np.float64)
vorticity_mag_mesh = np.zeros((ny,nz),dtype=np.float64)
index_global = 0
for j in range(0,nz):
    for i in range(0,ny):
        y_mesh[i,j] = data[index_global,0]
        z_mesh[i,j] = data[index_global,1]
        vorticity_mag_mesh[i,j] = data[index_global,2]
        index_global += 1
X = y_mesh
Y = z_mesh
Z = vorticity_mag_mesh

figure_filename="vorticity_mag-1"
# X=data[:,0]
# Y=data[:,1]
# Z=data[:,2]

# plotting code
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import rc as matplotlibrc
matplotlibrc('text.latex', preamble='\\usepackage{color}')
matplotlibrc('text', usetex=True)
matplotlibrc('font', family='serif')

# Figure parameters
fig_directory = "."
figure_filetype = "pdf"
# Font sizes
axisTitle_FontSize = 16
axisTickLabel_FontSize = 14
legend_fontSize = 16

y_label = "$z^{*}$"
x_label = "$y^{*}$"
z_label = "Nondimensional Vorticity Magnitude, $|\\mathbf{\\omega}^{*}|$"

print('Plotting: ' + figure_filename)
fig, ax = plt.subplots(figsize=(9,6))
# fig = plt.figure(figsize=(6,6))
plt.xlim([np.amin(X),np.amax(X)])
plt.ylim([np.amin(Y),np.amax(Y)])
ax.set_xlabel(x_label,fontsize=axisTitle_FontSize)
ax.set_ylabel(y_label,rotation=90,fontsize=axisTitle_FontSize)
plt.setp(ax.get_xticklabels(),fontsize=axisTickLabel_FontSize); plt.setp(ax.get_yticklabels(),fontsize=axisTickLabel_FontSize);
ax.set_title("TGV at Re$_{\\infty}=1600$, P$3$, $256^{3}$ DOFs, CFL$=0.30$, $x^{*}=-\\pi$",fontsize=axisTitle_FontSize)
minZ=1
maxZ=15.0


# levels = [1, 3, 5, 10, 15]
# levels_lines = levels#[1, 5, 10, 15]
levels = [1.5, 3, 4.5, 6, 7.5, 9, 10.5, 12, 13.5]
levels_lines = [1, 3, 6, 9, 12, 15]
# levels_lines = [1, 5, 10, 20, 30]
# levels = [1, 5, 10, 20, 30]
# levels = [3, 5, 7, 10, 15, 20, 25, 30]
# CS3.cmap.set_under('yellow')
# CS3.cmap.set_over('cyan')

# cs = plt.contourf(X,Y,Z, np.linspace(np.amin(Z),np.amax(Z),100),cmap='rainbow')

cs = plt.contourf(X,Y,Z, levels,cmap='rainbow',extend="both")

cs2 = plt.contour(X, Y, Z, levels_lines,colors=('k',),linewidths=(1,))
# cs = plt.contour(X, Y, Z, levels_lines,cmap='rainbow',linewidths=(.3,))
# cs = plt.contourf([X,Y,],Z,cmap='rainbow')

cbar = fig.colorbar(cs,ticks=levels)
# cbar = fig.colorbar(cs,ticks=np.linspace(np.amin(Z),np.amax(Z),8))
cbar.set_label(z_label,fontsize=axisTickLabel_FontSize)
# cbar.ax.tick_params(fontsize=axisTickLabel_FontSize)
plt.setp(cbar.ax.get_yticklabels(),fontsize=axisTickLabel_FontSize)
ax.set_xlim([0,0.5*np.pi])
ax.set_ylim([0.5*np.pi,np.pi])

# Fix for the white lines between contour levels
for c in cs.collections:
    c.set_edgecolor("face")

# leg = plt.legend(loc='best', ncol=1, shadow=False, fancybox=True, fontsize=legend_fontSize, framealpha=1.0,edgecolor='inherit')

plt.tight_layout()
print('... Saving figure ...')
plt.savefig(fig_directory+"/"+figure_filename+'.'+figure_filetype,format=figure_filetype,dpi=500)
plt.close()

# # qp.plotfield(xdata=x,ydata=y,udata=u,vdata=v,ylabel="$y$",xlabel="$x$",
# #     # title_label="DHIT Initial Velocity Field at $z=\\pi$\n $24^{3}$ DOF (P5, $N_{el}=4^{3}$)",
# #     fig_directory=".",figure_filename="velocity_field_48_init",
# #     xlimits=[0.0,2.0*np.pi],ylimits=[0.0,2.0*np.pi])
