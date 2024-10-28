import os;CURRENT_PATH = os.path.split(os.path.realpath(__file__))[0]+"/";
import sys
sys.path.append(CURRENT_PATH+"../../src");
import numpy as np
import pandas as pd
sys.path.append(CURRENT_PATH+"../../submodules/quickplotlib/lib"); import quickplotlib as qp
from sys import platform
if platform == "linux" or platform == "linux2":
    # linux
    filesystem="/media/julien/Samsung_T5/"
elif platform == "darwin":
    # OS X
    filesystem="/Volumes/Samsung_T5/"
#-----------------------------------------------------
#-----------------------------------------------------
# define functions
#-----------------------------------------------------
from scipy.interpolate import splrep, splev

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
def plot_mach_number_profile(files,labels_,DOF):
    global subdirectories, filenames, labels, black_line_flag, \
    dashed_line_flag, figure_filename_postfix, figure_title, \
    ylimits_kinetic_energy_input, ylimits_dissipation_input, \
    log_axes_input, legend_on_input, legend_inside_input, \
    plot_reference_result, nlegendcols_input, \
    figure_subdirectory, data_directory_base, figure_directory_base
    #-----------------------------------------------------
    # data store
    y_store = []
    mach_store = []

    # reference data
    y, scalar = np.loadtxt("./data/chapelier2024/reference_teno6_512_mach_number_profile.txt",skiprows=1,dtype=np.float64,unpack=True,delimiter=",")
    y_store.append(y)
    mach_store.append(scalar)
    labels.append("Ref. TENO6 $512^3$ DOF\n[Chapelier et al.]")
    # labels.append("$512^3$ DOF\n[Chapelier et al.]")

    # load results
    scalar_field_name = "mach_number"
    for i,input_filename in enumerate(files):
        # Load values
        df = pd.read_csv(data_directory_base+input_filename,sep=",",header=0)
        scalar = df[scalar_field_name].to_numpy()
        y = df["Points:1"].to_numpy()
        # y = y-np.pi() # shift data for comparison
        # store the data
        y_store.append(y)
        mach_store.append(scalar)
        labels.append(labels_[i])

    #-----------------------------------------------------
    clr_input_store = ['k','tab:blue','tab:red','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    mrkr_input_store = ['o','None','None','None','None','None','None','None','None']
    lnstl_input_store = ['None','solid','solid','solid','solid','dashed','solid','dashed','solid']
    # if(plotting_subsonic_result):
    #     lnstl_input_store = ['None','solid','dashed','solid','solid','dashed','solid','dashed','solid']
    # if(DOF==256):
    #     lnstl_input_store = ['None','solid','dashed','solid','solid','dashed','solid','dashed','solid']
    # if(DOF==128):
    #     lnstl_input_store = ['None','solid','solid','solid','dashed','dashed','solid','dashed','solid']

    qp.plotfxn(xdata=y_store,
            ydata=mach_store,
            xlabel='$y^{*}$',#=\\frac{1}{\\rho_{\\infty}V_{\\infty}^{2}|\\Omega|}\\int_{\\Omega}\\rho(u\\cdot\\u)d\\Omega$',
            ylabel='Local Mach Number, $M$',
            figure_filename=figure_subdirectory+'mach_number_profile'+figure_filename_postfix,
            title_label=figure_title,
            markers=False,
            legend_labels_tex=labels,
            black_lines=False,
            xlimits=[0,0.5*np.pi],
            ylimits=[0.2,1.8],
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
    return
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
# DOFs: 64^3 | All results
#-----------------------------------------------------
if(True):
    #-----------------------------------------------------
    # clr_input = ['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2024_JCP/"
    date_for_runs="."
    figure_subdirectory="./"
    # figure_title = "TGV at Re$_{\\infty}=1600$, $256^{3}$ DOFs, CFL=$0.10$" # comment to turn off
    figure_filename_postfix = "_64"
    legend_inside_input=True
    plot_reference_result=True
    plot_PHiLiP_DNS_result_as_reference_input=False
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "supersonic_viscous_TGV_ILES_NSFR_cHU_Ra_2PF_GLL_OI-0_dofs0064_p7_procs128/paraview_mach_number_vs_y_t2point5.txt",\
    "supersonic_viscous_TGV_ILES_NSFR_cDG_Ra_2PF_GLL_OI-0_dofs0064_p7_procs128/paraview_mach_number_vs_y_t2point5.txt",\
    "supersonic_viscous_TGV_ILES_NSFR_cPlus_Ra_2PF_GLL_OI-0_dofs0064_p3_procs128/paraview_mach_number_vs_y_t2point5.txt",\
    "supersonic_viscous_TGV_ILES_NSFR_cDG_Ra_2PF_GLL_OI-0_dofs0064_p3_procs128/paraview_mach_number_vs_y_t2point5.txt",\
    "supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0064_p3_procs128/paraview_mach_number_vs_y_t2point5.txt",\
    "supersonic_viscous_TGV_ILES_NSFR_cDG_KG_2PF_GLL_OI-0_dofs0064_p3_procs128/paraview_mach_number_vs_y_t2point5.txt",\
    ]
    # labels
    labels_for_plot=[\
    "p$7$ $c_{HU}$ NSFR.CH$_{\\mathrm{RA}}$",\
    "p$7$ $c_{DG}$ NSFR.CH$_{\\mathrm{RA}}$",\
    "p$3$ $c_{+}$ NSFR.CH$_{\\mathrm{RA}}$",\
    "p$3$ $c_{DG}$ NSFR.CH$_{\\mathrm{RA}}$",\
    "p$3$ $c_{DG}$ NSFR.IR",\
    "p$3$ $c_{DG}$ NSFR.KG",\
    ]
    plot_mach_number_profile(subdirectories_for_plot,labels_for_plot,64)

#=====================================================
# DOFs: 128^3 | All results
#-----------------------------------------------------
if(True):
    #-----------------------------------------------------
    # clr_input = ['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2024_JCP/"
    date_for_runs="."
    figure_subdirectory="./"
    # figure_title = "TGV at Re$_{\\infty}=1600$, $256^{3}$ DOFs, CFL=$0.10$" # comment to turn off
    figure_filename_postfix = "_128"
    legend_inside_input=True
    plot_reference_result=True
    plot_PHiLiP_DNS_result_as_reference_input=False
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "supersonic_viscous_TGV_ILES_NSFR_cDG_Ra_2PF_GLL_OI-0_dofs0128_p3_procs512/paraview_mach_number_vs_y_t2point5.txt",\
    "supersonic_viscous_TGV_ILES_NSFR_cPlus_Ra_2PF_GLL_OI-0_dofs0128_p3_procs512/paraview_mach_number_vs_y_t2point5.txt",\
    "supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0128_p7_procs512/paraview_mach_number_vs_y_t2point5.txt",\
    ]
    # labels
    labels_for_plot=[\
    "p$3$ $c_{DG}$ NSFR.CH$_{\\mathrm{RA}}$",\
    "p$3$ $c_{+}$ NSFR.CH$_{\\mathrm{RA}}$",\
    "p$7$ $c_{DG}$ NSFR.CH$_{\\mathrm{RA}}$",\
    ]
    plot_mach_number_profile(subdirectories_for_plot,labels_for_plot,128)
#=====================================================
# DOFs: 256^3 | All results
#-----------------------------------------------------
if(True):
    #-----------------------------------------------------
    # clr_input = ['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2024_JCP/"
    date_for_runs="."
    figure_subdirectory="./"
    # figure_title = "TGV at Re$_{\\infty}=1600$, $256^{3}$ DOFs, CFL=$0.10$" # comment to turn off
    figure_filename_postfix = "_256"
    legend_inside_input=True
    plot_reference_result=True
    plot_PHiLiP_DNS_result_as_reference_input=False
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "supersonic_viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0256_p7_procs512/paraview_mach_number_vs_y_t2point5.txt",\
    ]
    # labels
    labels_for_plot=[\
    "p$7$ $c_{DG}$ NSFR.CH$_{\\mathrm{RA}}$",\
    ]
    plot_mach_number_profile(subdirectories_for_plot,labels_for_plot,256)