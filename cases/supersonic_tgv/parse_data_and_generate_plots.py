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
    time_store = []
    kinetic_energy_store = []
    solenoidal_dissipation_store = []
    dilatational_dissipation_store = []
    #-----------------------------------------------------
    time, kinetic_energy = np.loadtxt("./data/chapelier2024/kinetic_energy.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
    time_store.append(time)
    kinetic_energy_store.append(kinetic_energy)
    labels.append("Chapelier et al. ($2048^3$ DOFs)")
    black_line_flag.append(True)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    clr_input_store = ['k','tab:blue','tab:blue','tab:red','tab:red','tab:green','tab:green','tab:orange','tab:orange']#,'tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    mrkr_input_store = ['o','None','None','None','None','None','None','None']
    lnstl_input_store = ['None','solid','solid','dashed','solid','dashed','solid','dashed','solid']
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
        filename = data_directory_base+"/"+subdirectories[i]+"/"+filenames[i]
        time, kinetic_energy, enstrophy, vorticity_based_dissipation, pressure_dilatation_based_dissipation, strain_rate_based_dissipation, deviatoric_strain_rate_based_dissipation, solenoidal_dissipation, dilatational_dissipation = np.loadtxt(filename,skiprows=1,dtype=np.float64,unpack=True)
        time_store.append(time)
        kinetic_energy_store.append(kinetic_energy)
        solenoidal_dissipation_store.append(solenoidal_dissipation)
        dilatational_dissipation_store.append(dilatational_dissipation)

    qp.plotfxn(xdata=time_store,#[time,time],
            ydata=kinetic_energy_store,#[kinetic_energy,kolmogorov_slope],
            ylabel='Nondimensional Kinetic Energy, $K^{*}$',#=\\frac{1}{\\rho_{\\infty}V_{\\infty}^{2}|\\Omega|}\\int_{\\Omega}\\rho(u\\cdot\\u)d\\Omega$',
            xlabel='Nondimensional Time, $t^{*}$',
            figure_filename=figure_subdirectory+'kinetic_energy_vs_time'+figure_filename_postfix,
            title_label=figure_title,
            markers=False,
            legend_labels_tex=labels,
            black_lines=False,
            xlimits=[0,20.0],
            ylimits=[0.0,0.13],
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
            legend_fontSize=14,
            legend_location="lower left")#,

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
# DOFs: 256^3 | All results
#-----------------------------------------------------
if(True):
    #-----------------------------------------------------
    # clr_input = ['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2024_JCP/"
    date_for_runs="."
    figure_subdirectory="2024_JCP"
    # figure_title = "TGV at Re$_{\\infty}=1600$, $256^{3}$ DOFs, CFL=$0.10$" # comment to turn off
    figure_filename_postfix = "16p7"
    legend_inside_input=True
    plot_reference_result=True
    plot_PHiLiP_DNS_result_as_reference_input=False
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0128_p7_procs512",\
    ]
    # labels
    labels_for_plot=[\
    "$c_{DG}$ NSFR 16p$7$ ($128^2$ DOF)",\
    ]
    black_line_flag_for_plot=[False,False,False,False,False,False,False,False]
    dashed_line_flag_for_plot=[False,False,False,False,False,True,True]
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot)
