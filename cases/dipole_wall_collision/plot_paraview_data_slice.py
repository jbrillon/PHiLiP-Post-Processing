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
#-------------------------------------------------------------
def get_reference_result():
    path = filesystem+"NarvalFiles/2024_CSME/dipole_wall_collision/viscous_DWC_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs1024_p7/solution_files/"
    filename = path+"vorticity"+"."+"csv"
    return get_paraview_file_data(filename)
#-------------------------------------------------------------
def get_paraview_file_data(
    input_filename,
    scalar_field_name="vorticity_magnitude"):
    # Load values
    df = pd.read_csv(input_filename,sep=",",header=0)
    scalar = df[scalar_field_name].to_numpy()
    x = df["Points:0"].to_numpy()
    y = df["Points:1"].to_numpy()
    return [x,y,scalar]
#-------------------------------------------------------------
def plot_vorticity_plane(path,filename_without_extension,file_extension,fig_directory,fig_prepre_fix,
    subdivide=False,title_label=" ",fill_contour=False,plot_reference_result=False,
    shift_domain_by_minus_pi=False,plot_single_element_domain=False,nElements_per_direction=16,
    plot_at_wall=False,color_bar_off=False,plot_wall_element_edge=False,wall_element_size_dx=0.0):
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
    x,y,scalar = get_paraview_file_data(filename)

    scalar_max = 317.0 # np.amax(scalar)
    scalar_min = 0.0 #np.amin(scalar)
    print("Scalar max value: %1.6e" % scalar_max)
    print("Scalar min value: %1.6e" % scalar_min)

    if(shift_domain_by_minus_pi):
        y = y-np.pi
        z = z-np.pi
        for i in range(0,len(z)):
            if(z[i]>=0.0000000000000000e+00):
                z[i] -= 3.1415926535897931e+00
            else:
                z[i] += 3.1415926535897931e+00

    X = x
    Y = y
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
    # if(fill_contour):
    #     contour_levels = [1, 5, 10, 20, 30]
    # elif(plot_reference_result):
    #     contour_levels = [1, 5, 10, 20, 30]
    # else:
    #     contour_levels = np.linspace(1,30,20)
    if(fill_contour):
        contour_levels = np.linspace(scalar_min,scalar_max,6)
    else:
        contour_levels = np.arange(10.0,271.0,20.0)
        # contour_levels = np.linspace(scalar_min,scalar_max,6)

    # Figure parameters
    # fig_directory = "figures"
    figure_filetype = "pdf"
    # Font sizes
    axisTitle_FontSize = 16#16
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
    y_label = "$y$"
    x_label = "$x$"
    z_label = "$|\\omega|$"
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
    # else:
    #     plt.xlim([0.0,np.pi])
    #     plt.ylim([-np.pi,0.0])

    if(plot_at_wall):
        plt.ylim([0.0,0.5])
        plt.xlim([0.5,1.0])

    ax.set_xlabel(x_label,fontsize=axisTitle_FontSize)
    ax.set_ylabel(y_label,rotation=90,fontsize=axisTitle_FontSize)
    if(title_label!=" "):
        plt.title(title_label,fontsize=2*axisTitle_FontSize)
    plt.setp(ax.get_xticklabels(),fontsize=axisTickLabel_FontSize); plt.setp(ax.get_yticklabels(),fontsize=axisTickLabel_FontSize);

    if(plot_reference_result):
        [x,y,w] = get_reference_result()
        plt.tricontour(x,y,w,levels=contour_levels,colors=('k',),linewidths=(1.0,),linestyles=("dashed",))
        # plt.contour(x,y,w, levels = contour_levels, colors=('k',),linewidths=(1,))

    if(plot_reference_result):
        contour_line_clr_for_result = 'r'
    else:
        contour_line_clr_for_result = 'k'

    if(subdivide):
        if(fill_contour):
            cs = plt.tricontourf(tri_refi, z_test_refi, np.linspace(scalar_min,scalar_max,100),cmap='rainbow',extend="both")
        cs_lines = plt.tricontour(tri_refi,z_test_refi,levels=contour_levels,colors=(contour_line_clr_for_result,),linewidths=(0.5,)) # change to 1
    else:
        if(fill_contour):
            # cs = plt.tricontourf(X, Y, Z, np.linspace(np.amin(Z),np.amax(Z),100),cmap='rainbow')
            cs = plt.tricontourf(X, Y, Z, np.linspace(scalar_min,scalar_max,100),cmap='rainbow',extend="both")
            # cs = plt.tricontourf(X, Y, Z, np.linspace(0.0,32,100),cmap='rainbow')
        cs_lines = plt.tricontour(X,Y,Z,levels=contour_levels,colors=(contour_line_clr_for_result,),linewidths=(0.5,))
    '''
    for level in cs_lines.collections:
        for kp,lpath in reversed(list(enumerate(level.get_paths()))):
            # go in reversed order due to deletions!

            # include test for "smallness" of your choice here:
            # I'm using a simple estimation for the diameter based on the
            #    x and y diameter...
            verts = lpath.vertices # (N,2)-shape array of contour line coordinates
            diameter = np.max(verts.max(axis=0) - verts.min(axis=0))

            if diameter<0.035: # threshold to be refined for your actual dimensions!
                del(level.get_paths()[kp])  # no remove() for Path objects:(
    '''
    if(fill_contour):
        cbar = fig.colorbar(cs,ticks=contour_levels)
        # cbar = fig.colorbar(cs,ticks=np.linspace(np.amin(Z),np.amax(Z),8))
        cbar.set_label(z_label,fontsize=axisTickLabel_FontSize)
        # cbar.ax.tick_params(fontsize=axisTickLabel_FontSize)
        plt.setp(cbar.ax.get_yticklabels(),fontsize=axisTickLabel_FontSize)

        # Fix for the white lines between contour levels
        for c in cs.collections:
            c.set_edgecolor("face")

    if((plot_at_wall or color_bar_off) and (fill_contour==True)):
        cbar.remove()
    if(plot_wall_element_edge):
        ax.axvline(x=(1.0-wall_element_size_dx),linestyle="dashed",color="k",alpha=0.5)
    plt.tight_layout()
    print('\t ... Saving figure ...')
    plt.savefig(fig_directory+"/"+figure_filename+'.'+figure_filetype,format=figure_filetype,dpi=500)#,rasterized=True
    plt.close()
    print('\t     Saved.')
    print("---------------------------------------------")
    return

if(True):
    # compare solution at wall
    paths=(
    filesystem+"NarvalFiles/2024_CSME/dipole_wall_collision/viscous_DWC_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0192_p2/solution_files/",\
    filesystem+"NarvalFiles/2024_CSME/dipole_wall_collision/viscous_DWC_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0384_p5/solution_files/",\
    filesystem+"NarvalFiles/2024_CSME/dipole_wall_collision/viscous_DWC_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0512_p7/solution_files/",\
    filesystem+"NarvalFiles/2024_CSME/dipole_wall_collision/viscous_DWC_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs1024_p7/solution_files/",\
    filesystem+"NarvalFiles/2024_CSME/dipole_wall_collision/viscous_DWC_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0192_p2_stretched_mesh/solution_files/",\
    )
    file_extension="csv"
    n_paths=len(paths)
    fig_directory = "./figures/2024_CSME"
    labels_for_plot=[\
        "$192^2$ P$2$ Uniform Mesh",\
        "$384^2$ P$5$ Uniform Mesh",\
        "$512^2$ P$7$ Uniform Mesh",\
        "$1024^2$ P$7$ Uniform Mesh",\
        "$192^2$ P$2$ Stretched Mesh",\
    ]
    fig_prepre_fix = [\
        "t4_192_p2_NSFR_cDG_IR_GL_",\
        "t4_384_p5_NSFR_cDG_IR_GL_",\
        "t4_512_p7_NSFR_cDG_IR_GL_",\
        "t4_1024_p7_NSFR_cDG_IR_GL_",\
        "t4_192_p2_stretched_mesh_NSFR_cDG_IR_GL_",\
    ]
    subdivide_input=[True,False,False,False,True]
    # subdivide_input=[False,False,False,True]
    # subdivide_input=[False,False,False,False,False]
    # subdivide_input=[True,True,True,True,False]
    time = [0.0,0.1,0.2,0.3,0.35,0.5,0.6,0.8,1.0]
    wall_element_size_dx_store=[\
        (2.0/64.0),\
        (2.0/64.0),\
        (2.0/64.0),\
        (2.0/128.0),\
        (1.0-9.9879544973373413e-01),\
    ]
    for i in range(0,n_paths):
        filename_without_extension = "vorticity"
        plot_title = labels_for_plot[i] #"$t=%1.2f$" % time[4]
        plot_vorticity_plane(paths[i],filename_without_extension,file_extension,fig_directory,fig_prepre_fix[i],
            subdivide=subdivide_input[i],
            #title_label=plot_title,
            fill_contour=False,
            plot_reference_result=True,
            shift_domain_by_minus_pi=False,
            plot_at_wall=True,
            plot_single_element_domain=False,
            plot_wall_element_edge=False,wall_element_size_dx=wall_element_size_dx_store[i])

if(True):
    # simulation overview
    paths=(
    # filesystem+"NarvalFiles/2024_CSME/dipole_wall_collision/viscous_DWC_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0384_p5/solution_files/",\
    filesystem+"NarvalFiles/2024_CSME/dipole_wall_collision/viscous_DWC_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0512_p7/solution_files/",\
    )
    files_per_path=[9,9]#[9,9]

    file_extension="csv"
    n_paths=len(paths)

    fig_directory = "./figures/2024_CSME"
    labels_for_plot=[\
        # "$384^2$ P$5$ $c_{DG}$ NSFR.IR-GL",\
        "$512^2$ P$7$ $c_{DG}$ NSFR.IR-GL",\
    ]
    fig_prepre_fix = [\
        # "384_p5_NSFR_cDG_IR_GL_",\
        "512_p7_NSFR_cDG_IR_GL_",\
    ]
    shift_domain_by_minus_pi_store=[False,False]
    time = [0.0,0.1,0.2,0.3,0.35,0.5,0.6,0.8,1.0]

    for i in range(0,n_paths):
        for j in range(0,files_per_path[i]):
                filename_without_extension = "vorticity_%i" % j
                plot_title = "$t=%1.2f$" % time[j]
                plot_vorticity_plane(paths[i],filename_without_extension,file_extension,fig_directory,fig_prepre_fix[i],
                    subdivide=False,
                    #title_label=plot_title,
                    fill_contour=False,plot_reference_result=False,
                    shift_domain_by_minus_pi=shift_domain_by_minus_pi_store[i],
                    plot_single_element_domain=False,
                    color_bar_off=True)
