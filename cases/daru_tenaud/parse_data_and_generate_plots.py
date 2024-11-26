import os;CURRENT_PATH = os.path.split(os.path.realpath(__file__))[0]+"/";
import sys
sys.path.append(CURRENT_PATH+"../../src");
import numpy as np
sys.path.append(CURRENT_PATH+"../../submodules/quickplotlib/lib"); import quickplotlib as qp
from sys import platform
if platform == "linux" or platform == "linux2":
    # linux
    filesystem="/media/julien/Samsung_T5/"
elif platform == "darwin":
    # OS X
    filesystem="/Volumes/Samsung_T5/"
#-----------------------------------------------------
#=====================================================
# Input variables for plotting
#=====================================================
global subdirectories, filenames, labels, black_line_flag, \
dashed_line_flag, figure_filename_postfix, figure_title, \
ylimits_kinetic_energy_input, ylimits_dissipation_input, \
log_axes_input, legend_on_input, legend_inside_input, \
plot_reference_result, nlegendcols_input, \
figure_subdirectory, data_directory_base, figure_directory_base
#=====================================================
def plot_for_presentation(
    subdirectories_for_plot,
    labels_for_plot,
    black_line_flag_for_plot,
    dashed_line_flag_for_plot):
    
    global subdirectories, filenames, labels, black_line_flag, \
    dashed_line_flag, figure_filename_postfix, figure_title, \
    ylimits_kinetic_energy_input, ylimits_dissipation_input, \
    log_axes_input, legend_on_input, legend_inside_input, \
    plot_reference_result, nlegendcols_input, \
    figure_subdirectory, data_directory_base, figure_directory_base
    #-----------------------------------------------------
    # data store
    triple_point_x_location_value_store = []
    triple_point_x_location_time_store = []
    triple_point_y_location_value_store = []
    triple_point_y_location_time_store = []
    vorticity_near_wall_value_store = []
    vorticity_near_wall_x_store = []
    density_along_wall_value_store = []
    density_along_wall_x_store = []
    density_near_wall_value_store = []
    density_near_wall_x_store = []
    skin_friction_value_store = []
    skin_friction_x_store = []
    #-----------------------------------------------------
    # load reference data
    x, density = np.loadtxt("./data/reference/daru2009numerical/density_Re1000_along_wall.csv",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
    density_along_wall_value_store.append(density)
    density_along_wall_x_store.append(x)
    x, density = np.loadtxt("./data/reference/daru2009numerical/density_Re1000_y0.05.csv",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
    density_near_wall_value_store.append(density)
    density_near_wall_x_store.append(x)
    x, vorticity = np.loadtxt("./data/reference/daru2009numerical/vorticity_Re1000_y0.05.csv",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
    vorticity_near_wall_value_store.append(vorticity)
    vorticity_near_wall_x_store.append(x)
    x, skin_friction = np.loadtxt("./data/reference/daru2009numerical/skinFriction_Re1000.csv",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
    skin_friction_value_store.append(skin_friction)
    skin_friction_x_store.append(x)
    time, location = np.loadtxt("./data/reference/daru2009numerical/triplePointLocation_x.csv",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
    triple_point_x_location_value_store.append(location)
    triple_point_x_location_time_store.append(time)
    time, location = np.loadtxt("./data/reference/daru2009numerical/triplePointLocation_y.csv",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
    triple_point_y_location_value_store.append(location)
    triple_point_y_location_time_store.append(time)
    labels.append("OSMP7 [Daru \\& Tenaud]\n $4000\\times2000$ DOF")
    black_line_flag.append(True)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    clr_input_store = ['k','tab:blue','tab:red','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    mrkr_input_store = ['None','None','None','None','None','None','None','None']
    lnstl_input_store = ['solid','solid','solid','solid','solid','dashed','solid','dashed','solid']
    #-----------------------------------------------------
    number_of_result_curves=len(subdirectories_for_plot)
    for i in range(0,number_of_result_curves):
        figure_filename_postfix_input=figure_filename_postfix
        #-----------------------------------------------------
        subdirectories.append(subdirectories_for_plot[i])
        labels.append(labels_for_plot[i])
        black_line_flag.append(black_line_flag_for_plot[i])
        dashed_line_flag.append(dashed_line_flag_for_plot[i])
        filenames.append("turbulent_quantities.txt")
        #-----------------------------------------------------
        # load file
        directory = data_directory_base+"/"+subdirectories[i]+"/"
        x, density = np.loadtxt(directory+"density_along_wall.csv",skiprows=0,dtype=np.float64,unpack=True,delimiter=",")
        density_along_wall_value_store.append(density)
        density_along_wall_x_store.append(x)
        x, density = np.loadtxt(directory+"density_along_y0point05.csv",skiprows=0,dtype=np.float64,unpack=True,delimiter=",")
        density_near_wall_value_store.append(density)
        density_near_wall_x_store.append(x)
        x, vorticity = np.loadtxt(directory+"vorticity_along_y0point05.csv",skiprows=0,dtype=np.float64,unpack=True,delimiter=",")
        vorticity_near_wall_value_store.append(vorticity/1000.0) # WARNING SCALED RESULTS MANUALLY
        vorticity_near_wall_x_store.append(x)
        # x, skin_friction = np.loadtxt(directory+"skinFriction_along_wall.csv",skiprows=0,dtype=np.float64,unpack=True,delimiter=",")
        # skin_friction_value_store.append(skin_friction)
        # skin_friction_x_store.append(x)
        time, x, y = np.loadtxt(directory+"triple_point_vs_time.txt",skiprows=0,dtype=np.float64,unpack=True)
        triple_point_x_location_value_store.append(x)
        triple_point_x_location_time_store.append(time)
        triple_point_y_location_value_store.append(y)
        triple_point_y_location_time_store.append(time)

    qp.plotfxn(xdata=density_along_wall_x_store,
            ydata=density_along_wall_value_store,
            ylabel='Density $\\rho^{*}$ along wall',
            xlabel='Nondimensional $x$-coordinate, $x^{*}$',
            figure_filename=figure_subdirectory+'density_along_wall'+figure_filename_postfix,
            title_label=figure_title,
            markers=False,
            legend_labels_tex=labels,
            black_lines=False,
            xlimits=[0.3,1.0],
            ylimits=[20,120],
            log_axes=log_axes_input,
            which_lines_black=black_line_flag,
            which_lines_dashed=dashed_line_flag,
            which_lines_only_markers=[0],
            legend_on=legend_on_input,
            legend_inside=legend_inside_input,
            nlegendcols=nlegendcols_input,
            figure_size=(6,6),
            transparent_legend=True,#transparent_legend_input,
            legend_border_on=False,
            grid_lines_on=False,
            clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store,
            legend_fontSize=12,#14
            legend_location="best")

    qp.plotfxn(xdata=density_near_wall_x_store,
            ydata=density_near_wall_value_store,
            ylabel='Density $\\rho^{*}$ near wall at $y^{*}=0.05$',
            xlabel='Nondimensional $x$-coordinate, $x^{*}$',
            figure_filename=figure_subdirectory+'density_near_wall'+figure_filename_postfix,
            title_label=figure_title,
            markers=False,
            legend_labels_tex=labels,
            black_lines=False,
            xlimits=[0.3,1.0],
            ylimits=[10,110],
            log_axes=log_axes_input,
            which_lines_black=black_line_flag,
            which_lines_dashed=dashed_line_flag,
            which_lines_only_markers=[0],
            legend_on=legend_on_input,
            legend_inside=legend_inside_input,
            nlegendcols=nlegendcols_input,
            figure_size=(6,6),
            transparent_legend=True,#transparent_legend_input,
            legend_border_on=False,
            grid_lines_on=False,
            clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store,
            legend_fontSize=12,#14
            legend_location="best")

    qp.plotfxn(xdata=vorticity_near_wall_x_store,
            ydata=vorticity_near_wall_value_store,
            ylabel='Vorticity near wall at $y^{*}=0.05$',
            xlabel='Nondimensional $x$-coordinate, $x^{*}$',
            figure_filename=figure_subdirectory+'vorticity_near_wall'+figure_filename_postfix,
            title_label=figure_title,
            markers=False,
            legend_labels_tex=labels,
            black_lines=False,
            xlimits=[0.5,1.0],
            # ylimits=[20,120],
            log_axes=log_axes_input,
            which_lines_black=black_line_flag,
            which_lines_dashed=dashed_line_flag,
            which_lines_only_markers=[0],
            legend_on=legend_on_input,
            legend_inside=legend_inside_input,
            nlegendcols=nlegendcols_input,
            figure_size=(6,6),
            transparent_legend=True,#transparent_legend_input,
            legend_border_on=False,
            grid_lines_on=False,
            clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store,
            legend_fontSize=12,#14
            legend_location="best")

    # qp.plotfxn(xdata=skin_friction_x_store,
    #         ydata=skin_friction_value_store,
    #         ylabel='Skin friction',
    #         xlabel='Nondimensional $x$-coordinate, $x^{*}$',
    #         figure_filename=figure_subdirectory+'skin_friction_wall'+figure_filename_postfix,
    #         title_label=figure_title,
    #         markers=False,
    #         legend_labels_tex=labels,
    #         black_lines=False,
    #         xlimits=[0.3,1.0],
    #         ylimits=[20,120],
    #         log_axes=log_axes_input,
    #         which_lines_black=black_line_flag,
    #         which_lines_dashed=dashed_line_flag,
    #         which_lines_only_markers=[0],
    #         legend_on=legend_on_input,
    #         legend_inside=legend_inside_input,
    #         nlegendcols=nlegendcols_input,
    #         figure_size=(6,6),
    #         transparent_legend=True,#transparent_legend_input,
    #         legend_border_on=False,
    #         grid_lines_on=False,
    #         clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store,
    #         legend_fontSize=12,#14
    #         legend_location="best")

    mrkr_input_store = ['o','^','d','s','+','d','v','>','<']
    qp.plotfxn(xdata=triple_point_x_location_time_store,
            ydata=triple_point_x_location_value_store,
            ylabel='Nondimensional $x$-coordinate, $x^{*}$',
            xlabel='Nondimensional Time, $t^{*}$',
            figure_filename=figure_subdirectory+'triple_point_x_location'+figure_filename_postfix,
            title_label=figure_title,
            markers=False,
            legend_labels_tex=labels,
            black_lines=False,
            xlimits=[0.4,1.0],
            ylimits=[0.55,0.95],
            log_axes=log_axes_input,
            which_lines_black=black_line_flag,
            which_lines_dashed=dashed_line_flag,
            which_lines_only_markers=[0],
            legend_on=legend_on_input,
            legend_inside=legend_inside_input,
            nlegendcols=nlegendcols_input,
            figure_size=(6,6),
            transparent_legend=True,#transparent_legend_input,
            legend_border_on=False,
            grid_lines_on=False,
            clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store,
            legend_fontSize=12,#14
            legend_location="best")

    qp.plotfxn(xdata=triple_point_y_location_time_store,
            ydata=triple_point_y_location_value_store,
            ylabel='Nondimensional $y$-coordinate, $y^{*}$',
            xlabel='Nondimensional Time, $t^{*}$',
            figure_filename=figure_subdirectory+'triple_point_y_location'+figure_filename_postfix,
            title_label=figure_title,
            markers=False,
            legend_labels_tex=labels,
            black_lines=False,
            xlimits=[0.4,1.0],
            ylimits=[0.07,0.16],
            log_axes=log_axes_input,
            which_lines_black=black_line_flag,
            which_lines_dashed=dashed_line_flag,
            which_lines_only_markers=[0],
            legend_on=legend_on_input,
            legend_inside=legend_inside_input,
            nlegendcols=nlegendcols_input,
            figure_size=(6,6),
            transparent_legend=True,#transparent_legend_input,
            legend_border_on=False,
            grid_lines_on=False,
            clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store,
            legend_fontSize=12,#14
            legend_location="best")
#-----------------------------------------------------
#=====================================================
def reinit_inputs():
    global subdirectories, filenames, labels, black_line_flag, \
    dashed_line_flag, figure_filename_postfix, figure_title, \
    ylimits_kinetic_energy_input, ylimits_dissipation_input, \
    log_axes_input, legend_on_input, legend_inside_input, \
    plot_reference_result, nlegendcols_input, \
    figure_subdirectory, data_directory_base, figure_directory_base, \
    smoothing_input, plot_PHiLiP_DNS_result_as_reference_input

    subdirectories = []
    filenames = []
    labels = []
    black_line_flag = []
    dashed_line_flag = []
    figure_filename_postfix = "" # default
    figure_title = "" # default
    ylimits_kinetic_energy_input = [] # default
    ylimits_dissipation_input = [] # default
    log_axes_input=None # default
    legend_on_input=True # default
    legend_inside_input=False # default
    plot_reference_result=False # default
    nlegendcols_input=1
    figure_subdirectory="" # default
    data_directory_base = "/Users/Julien/julien_phd/post_processing/data/taylor_green_vortex"
    # figure_directory_base = "/Users/Julien/julien_phd/post_processing/figures/taylor_green_vortex"
    figure_directory_base = "figures"
    plot_PHiLiP_DNS_result_as_reference_input=True # default
#=====================================================
#-----------------------------------------------------

#=====================================================
# All results
#-----------------------------------------------------
if(True):
    #-----------------------------------------------------
    # clr_input = ['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2024_JCP/daru_tenaud"
    date_for_runs="."
    figure_subdirectory="./"
    # figure_title = "DT at Re$_{\\infty}=1000$" # comment to turn off
    figure_filename_postfix = "_all"
    legend_inside_input=True
    plot_reference_result=True
    plot_PHiLiP_DNS_result_as_reference_input=False
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "viscous_DT_Re1000_p3_cDG_Ra_2PF_GLL_OI-0_p3_250x125_elements",\
    # "viscous_DT_Re1000_p3_cDG_Ra_2PF_GLL_OI-0_p3_500x250_elements",\
    ]
    # labels
    labels_for_plot=[\
    "p$3$ $c_{DG}$ NSFR.CH$_{\\mathrm{RA}}$-GLL-Roe-PPL\n $1000\\times500$ DOF",\
    # "p$3$ $c_{DG}$ NSFR.CH$_{\\mathrm{RA}}$",\
    ]
    black_line_flag_for_plot=[False,False,False,False,False,False,False,False]
    dashed_line_flag_for_plot=[False,False,False,False,False,True,True]
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot)
