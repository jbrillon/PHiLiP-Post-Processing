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

from struct import unpack # for loading the reference data

def get_reference_result():
    ''' 
    Reference: W. M. van Rees W. M., A. Leonard, D. I. Pullin and P. Koumoutsakos, A comparison of vortex and pseudo-spectral methods for the simulation of periodic vortical flows at high Reynolds numbers, J. Comput. Phys., 230(2011), 2794-2805.
    Source: (1) http://www.as.dlr.de/hiocfd/case_c3.5.pdf 
            (2) http://www.as.dlr.de/hiocfd/
    '''
    Pi = np.arccos(-1)
    x = np.arange(-Pi,Pi-2.*Pi/512.,2.*Pi/512.)
    y = np.arange(-Pi,Pi-2.*Pi/512.,2.*Pi/512.)

    size = (len(x),len(y))

    w = np.zeros(size)

    file = open('./data/van_rees/wn_slice_x0_08000.out', 'rb')
    data = file.read()

    n = 0
    for i in range(len(x)):
      for j in range(len(y)):
        a = unpack('d', data[n:n+8])
        w[j][i] = a[0]
        n = n+8
    return [x,y,w]

def plot_vorticity_plane(path,filename_without_extension,file_extension,fig_directory,fig_prepre_fix,
    subdivide=False,title_label=" ",fill_contour=False,plot_reference_result=False,
    shift_domain_by_minus_pi=False,plot_single_element_domain=False,nElements_per_direction=16):
    # TO DO: move this somewhere else eventually
    figure_filename = fig_prepre_fix+filename_without_extension
    if(subdivide):
        figure_filename += "_with_subdivision"
    else:
        figure_filename += "_no_subdivision"
    if(fill_contour):
        figure_filename += "_filled"
    else:
        figure_filename += "_no_fill"

    if(plot_reference_result):
        figure_filename += "_comparison_to_reference"
    # load data
    filename = path+filename_without_extension+"."+file_extension
    y,z,scalar = np.loadtxt(filename,unpack=True,dtype=np.float64)

    if(shift_domain_by_minus_pi):
        y = y-np.pi
        z = z-np.pi
        for i in range(0,len(z)):
            if(z[i]>=0.0000000000000000e+00):
                z[i] -= 3.1415926535897931e+00
            else:
                z[i] += 3.1415926535897931e+00

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

    # contour_levels = [1, 3, 5, 10, 15, 20, 30]
    if(fill_contour):
        contour_levels = [1, 5, 10, 20, 30]
    elif(plot_reference_result):
        contour_levels = [1, 5, 10, 20, 30]
    else:
        contour_levels = np.linspace(1,30,20)

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
    if(plot_single_element_domain):
        element_edges = np.linspace(-np.pi,np.pi,nElements_per_direction)
        x_edge_index = int(0.33*nElements_per_direction)
        y_edge_index = int(0.75*nElements_per_direction)
        plt.xlim([element_edges[x_edge_index],element_edges[x_edge_index+1]])
        plt.ylim([element_edges[y_edge_index],element_edges[y_edge_index+1]])
    else:
        plt.xlim([0.0,np.pi])
        plt.ylim([-np.pi,0.0])

    ax.set_xlabel(x_label,fontsize=axisTitle_FontSize)
    ax.set_ylabel(y_label,rotation=90,fontsize=axisTitle_FontSize)
    if(title_label!=" "):
        plt.title(title_label,fontsize=axisTitle_FontSize)
    plt.setp(ax.get_xticklabels(),fontsize=axisTickLabel_FontSize); plt.setp(ax.get_yticklabels(),fontsize=axisTickLabel_FontSize);

    if(subdivide):
        if(fill_contour):
            cs = plt.tricontourf(tri_refi, z_test_refi, np.linspace(np.amin(Z),np.amax(Z),100),cmap='rainbow',extend="both")
        if(plot_reference_result):
            contour_line_clr_for_result = 'r'
        else:
            contour_line_clr_for_result = 'k'
        cs_lines = plt.tricontour(tri_refi,z_test_refi,levels=contour_levels,colors=(contour_line_clr_for_result,),linewidths=(1,))
    else:
        if(fill_contour):
            # cs = plt.tricontourf(X, Y, Z, np.linspace(np.amin(Z),np.amax(Z),100),cmap='rainbow')
            # cs = plt.tricontourf(X, Y, Z, np.linspace(np.amin(contour_levels),np.amax(contour_levels),100),cmap='rainbow')
            cs = plt.tricontourf(X, Y, Z, np.linspace(0.0,32,100),cmap='rainbow')
        cs_lines = plt.tricontour(X,Y,Z,levels=contour_levels,colors=('k',),linewidths=(0.5,))

    for level in cs_lines.collections:
        for kp,lpath in reversed(list(enumerate(level.get_paths()))):
            # go in reversed order due to deletions!

            # include test for "smallness" of your choice here:
            # I'm using a simple estimation for the diameter based on the
            #    x and y diameter...
            verts = lpath.vertices # (N,2)-shape array of contour line coordinates
            diameter = np.max(verts.max(axis=0) - verts.min(axis=0))

            if diameter<0.15: # threshold to be refined for your actual dimensions!
                del(level.get_paths()[kp])  # no remove() for Path objects:(

    if(fill_contour):
        cbar = fig.colorbar(cs,ticks=contour_levels)
        # cbar = fig.colorbar(cs,ticks=np.linspace(np.amin(Z),np.amax(Z),8))
        cbar.set_label(z_label,fontsize=axisTickLabel_FontSize)
        # cbar.ax.tick_params(fontsize=axisTickLabel_FontSize)
        plt.setp(cbar.ax.get_yticklabels(),fontsize=axisTickLabel_FontSize)

        # Fix for the white lines between contour levels
        for c in cs.collections:
            c.set_edgecolor("face")

    if(plot_reference_result):
        [x,y,w] = get_reference_result()
        plt.contour(x,y,w, levels = contour_levels, colors=('k',),linewidths=(1,))

    cbar.remove()

    plt.tight_layout()
    print('\t ... Saving figure ...')
    plt.savefig(fig_directory+"/"+figure_filename+'.'+figure_filetype,format=figure_filetype,dpi=500)#,rasterized=True
    plt.close()
    print('\t     Saved.')
    print("---------------------------------------------")
    return

# lets check the 96^3 DOFs P5 results
paths=(
filesystem+"NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/flow_field_files/",\
filesystem+"NarvalFiles/2023_JCP/spectra_fix/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/flow_field_files/",\
)
files_per_path=[2,2]

file_extension="dat"
n_paths=len(paths)

fig_directory = "./figures/2023_JCP"
labels_for_plot=[\
    "$96^3$ P$5$ $c_{DG}$ NSFR.IR-GL",\
    "Oversampled $96^3$ P$5$ $c_{DG}$ NSFR.IR-GL",\
]
fig_prepre_fix = [\
    "new_96_p5_NSFR_cDG_IR_GL_",\
    "new_96_p5_NSFR_cDG_IR_GL_oversampled_",\
]
shift_domain_by_minus_pi_store=[False,True]

for i in range(0,n_paths):
    for j in range(1,files_per_path[i]):
            filename_without_extension = "velocity_vorticity-%i_reordered_slice_interface_averaged" % j
            plot_title = labels_for_plot[i]
            plot_vorticity_plane(paths[i],filename_without_extension,file_extension,fig_directory,fig_prepre_fix[i],
                subdivide=False,title_label=plot_title,fill_contour=True,plot_reference_result=False,
                shift_domain_by_minus_pi=shift_domain_by_minus_pi_store[i],
                plot_single_element_domain=False)
exit()
# DNS STUFF FROM PAPER BELOW
# output_solution_fixed_times_string = [0.0,4.0,5.0,8.0,9.0,10.0,12.0,15.0,16.0,20.0]

# filename="/home/julien/NarvalFiles/2023_JCP/filtered_dns_viscous_tgv/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p7_procs1024/solution_files/vorticity_mag_slice_x0_plane_t_3_quadrant.txt"

# filename="/home/julien/NarvalFiles/2023_JCP/verification/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-4_dofs0256_p3_CFL-0.15_procs1024/solution_files/vorticity_mag_slice_x0_plane_t_3_quadrant.txt"

#-------------------------------------------------------------
# paths=(
# "/home/julien/NarvalFiles/2023_JCP/filtered_dns_viscous_tgv/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p7_procs1024/solution_files/",\
# "/home/julien/NarvalFiles/2023_JCP/verification/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p3_procs1024/solution_files/",\
# "/home/julien/NarvalFiles/2023_JCP/verification/viscous_TGV_ILES_NSFR_cDG_IR_2PF-Roe_GL_OI-0_dofs0256_p3_procs1024/solution_files/",\
# "/home/julien/NarvalFiles/2023_JCP/verification/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-4_dofs0256_p3_CFL-0.15_procs1024/solution_files/",\
# )
# files_per_path=[10,10,7,7]

paths=(
# filesystem+"NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/solution_files/",\
filesystem+"NarvalFiles/2023_JCP/verification/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0256_p7_procs1024/solution_files/",\
)
files_per_path=[10]

prefix="vorticity_mag_slice_x0_plane_t_"
file_extension="txt"
scalar_field_name = "vorticity_magnitude"
n_paths=len(paths)
output_solution_fixed_times_string = [0.0,4.0,5.0,8.0,9.0,10.0,12.0,15.0,16.0,20.0]
fig_directory = "./figures/2023_JCP"
labels_for_plot=[\
    # "$256^3$ P$7$ $c_{DG}$ NSFR.IR-GL",\
    # "$256^3$ P$3$ $c_{DG}$ NSFR.IR-GL",\
    # "$256^3$ P$3$ $c_{DG}$ NSFR.IR-GL-Roe",\
    # "$256^3$ P$3$ Strong DG-Roe-GL-OI",\
    # "$96^3$ P$5$ $c_{DG}$ NSFR.IR-GL",\
    "$256^3$ P$7$ $c_{DG}$ NSFR.IR-GLL",\
]
fig_prepre_fix = [\
    # "256_p7_NSFR_cDG_IR_GL_",\
    # "256_p3_NSFR_cDG_IR_GL_",\
    # "256_p3_NSFR_cDG_IR_Roe_GL_",\
    # "256_p3_strong_DG_Roe_GL_OI-4_",\
    # "96_p5_NSFR_cDG_IR_GL_",\
    "256_p7_NSFR_cDG_IR_GLL_",\
]


for i in range(0,n_paths):
# for i in range(0,1):
    for j in range(0,files_per_path[i]):
    # for j in range(3,4):
        '''filename_without_extension = prefix+str(j)+"_quadrant"
        plot_title = labels_for_plot[i]+", $t^{*}=%i$" % output_solution_fixed_times_string[j]
        plot_vorticity_plane(paths[i],filename_without_extension,file_extension,fig_directory,fig_prepre_fix[i],
            subdivide=True,title_label=plot_title,fill_contour=False)'''
        if(output_solution_fixed_times_string[j]==8.0):
            filename_without_extension = prefix+str(j)+"_quadrant"
            plot_title = labels_for_plot[i]+", $t^{*}=%i$" % output_solution_fixed_times_string[j]
            plot_vorticity_plane(paths[i],filename_without_extension,file_extension,fig_directory,fig_prepre_fix[i],
                subdivide=True,title_label=plot_title,fill_contour=False,plot_reference_result=True)

