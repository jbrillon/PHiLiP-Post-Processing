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
import matplotlib.tri as tri

def plot_vorticity_plane(path,filename_without_extension,file_extension,fig_directory,fig_prepre_fix,
    subdivide=False,title_label=" "):
    # TO DO: move this somewhere else eventually
    figure_filename = fig_prepre_fix+filename_without_extension
    if(subdivide):
        figure_filename += "_with_subdivision"
    else:
        figure_filename += "_no_subdivision"

    # load data
    filename = path+filename_without_extension+"."+file_extension
    y,z,scalar = np.loadtxt(filename,unpack=True,dtype=np.float64)

    X = y
    Y = z
    Z = scalar
    # ----------------------------------------
    # ----------------------------------------
    # Smooth data; ref: https://matplotlib.org/stable/gallery/images_contours_and_fields/tricontour_smooth_user.html
    # ---------------------------------------- 
    if(subdivide):
        min_radius = np.amin(Z)
        # Now create the Triangulation.
        # (Creating a Triangulation without specifying the triangles results in the
        # Delaunay triangulation of the points.)
        triang = tri.Triangulation(X, Y)

        # Mask off unwanted triangles.
        triang.set_mask(np.hypot(X[triang.triangles].mean(axis=1),
                                 Y[triang.triangles].mean(axis=1))
                        < min_radius)

        refiner = tri.UniformTriRefiner(triang)
        tri_refi, z_test_refi = refiner.refine_field(Z, subdiv=3)

    # ----------------------------------------


    contour_levels = [1.2,9.0,11.0]
    colors=["0.25","0.25","0.25"]
    # Figure parameters
    # fig_directory = "figures"
    figure_filetype = "pdf"
    # Font sizes
    axisTitle_FontSize = 16
    axisTickLabel_FontSize = 14
    legend_fontSize = 16

    # plotting code
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib import rc as matplotlibrc
    matplotlibrc('text.latex', preamble='\\usepackage{color}')
    matplotlibrc('text', usetex=True)
    matplotlibrc('font', family='serif')
    print('Plotting: ' + fig_directory + "/" + figure_filename + "." + figure_filetype)
    y_label = "$z$"
    x_label = "$y$"
    z_label = "$|\\Omega|$"
    fig, ax = plt.subplots(figsize=(6,6))
    # fig = plt.figure(figsize=(6,6))
    # plt.xlim([np.amin(X),np.amax(X)]) # if not provided
    # plt.ylim([np.amin(Y),np.amax(Y)]) # if not provided
    plt.xlim([0.0,np.pi])
    plt.ylim([-np.pi,0.0])
    ax.set_xlabel(x_label,fontsize=axisTitle_FontSize)
    ax.set_ylabel(y_label,rotation=90,fontsize=axisTitle_FontSize)
    if(title_label!=" "):
        plt.title(title_label,fontsize=axisTitle_FontSize)
    plt.setp(ax.get_xticklabels(),fontsize=axisTickLabel_FontSize); plt.setp(ax.get_yticklabels(),fontsize=axisTickLabel_FontSize);

    if(subdivide):
        cs = plt.tricontourf(tri_refi, z_test_refi, np.linspace(np.amin(Z),np.amax(Z),100),cmap='rainbow',extend="both")
        # plt.tricontour(tri_refi, z_test_refi)
    else:
        cs = plt.tricontourf(X, Y, Z, np.linspace(np.amin(Z),np.amax(Z),100),cmap='rainbow')
        # plt.tricontour(X, Y, Z,levels=contour_levels,colors=colors)

    cbar = fig.colorbar(cs,ticks=np.arange(np.amin(Z),np.amax(Z),10.0))
    # cbar = fig.colorbar(cs,ticks=np.linspace(np.amin(Z),np.amax(Z),8))
    cbar.set_label(z_label,fontsize=axisTickLabel_FontSize)
    # cbar.ax.tick_params(fontsize=axisTickLabel_FontSize)
    plt.setp(cbar.ax.get_yticklabels(),fontsize=axisTickLabel_FontSize)

    # Fix for the white lines between contour levels
    for c in cs.collections:
        c.set_edgecolor("face")

    plt.tight_layout()
    print('\t ... Saving figure ...')
    plt.savefig(fig_directory+"/"+figure_filename+'.'+figure_filetype,format=figure_filetype,dpi=500)
    plt.close()
    print('\t     Saved.')
    print("---------------------------------------------")
    return


# output_solution_fixed_times_string = [0.0,4.0,5.0,8.0,9.0,10.0,12.0,15.0,16.0,20.0]

# filename="/home/julien/NarvalFiles/2023_JCP/filtered_dns_viscous_tgv/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p7_procs1024/solution_files/vorticity_mag_slice_x0_plane_t_3_quadrant.txt"

# filename="/home/julien/NarvalFiles/2023_JCP/verification/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-4_dofs0256_p3_CFL-0.15_procs1024/solution_files/vorticity_mag_slice_x0_plane_t_3_quadrant.txt"

#-------------------------------------------------------------
paths=(
"/home/julien/NarvalFiles/2023_JCP/filtered_dns_viscous_tgv/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p7_procs1024/solution_files/",\
"/home/julien/NarvalFiles/2023_JCP/verification/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p3_procs1024/solution_files/",\
"/home/julien/NarvalFiles/2023_JCP/verification/viscous_TGV_ILES_NSFR_cDG_IR_2PF-Roe_GL_OI-0_dofs0256_p3_procs1024/solution_files/",\
"/home/julien/NarvalFiles/2023_JCP/verification/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-4_dofs0256_p3_CFL-0.15_procs1024/solution_files/",\
)
files_per_path=[10,10,7,7]
prefix="vorticity_mag_slice_x0_plane_t_"
file_extension="txt"
scalar_field_name = "vorticity_magnitude"
n_paths=len(paths)
output_solution_fixed_times_string = [0.0,4.0,5.0,8.0,9.0,10.0,12.0,15.0,16.0,20.0]
fig_directory = "./figures/2023_JCP"
labels_for_plot=[\
    "$256^3$ P$7$ $c_{DG}$ NSFR.IR-GL",\
    "$256^3$ P$3$ $c_{DG}$ NSFR.IR-GL",\
    "$256^3$ P$3$ $c_{DG}$ NSFR.IR-GL-Roe",\
    "$256^3$ P$3$ Strong DG-Roe-GL-OI",\
    ]
fig_prepre_fix = [\
    "256_p7_NSFR_cDG_IR_GL_",\
    "256_p3_NSFR_cDG_IR_GL_",\
    "256_p3_NSFR_cDG_IR_Roe_GL_",\
    "256_p3_strong_DG_Roe_GL_OI-4_",\
]


for i in range(0,n_paths):
    for j in range(0,files_per_path[i]):
        filename_without_extension = prefix+str(j)+"_quadrant"
        plot_title = labels_for_plot[i]+", $t^{*}=%i$" % output_solution_fixed_times_string[j]
        plot_vorticity_plane(paths[i],filename_without_extension,file_extension,fig_directory,fig_prepre_fix[i],subdivide=False,title_label=plot_title)

