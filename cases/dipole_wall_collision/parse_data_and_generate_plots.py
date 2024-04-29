import os;CURRENT_PATH = os.path.split(os.path.realpath(__file__))[0]+"/";
import sys
sys.path.append(CURRENT_PATH+"../../src");
from plot_dipole_wall_collision_quantities import plot_periodic_turbulence

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
figure_subdirectory, data_directory_base, figure_directory_base, \
smoothing_input
#=====================================================
def get_smoothing_parameters_from_subdirectories(subdirectories_for_plot):
    smoothing_parameters_store=[]
    number_of_result_curves=len(subdirectories_for_plot)
    for i in range(0,number_of_result_curves):
        name=subdirectories_for_plot[i]
        if(name=="viscous_DWC_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0192_p2"):
            smoothing_parameter=200
        elif(name=="viscous_DWC_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0192_p2_small_dt"):
            smoothing_parameter=200 # copied
        elif(name=="viscous_DWC_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0384_p5"):
            smoothing_parameter=500
        elif(name=="viscous_DWC_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0512_p7"):
            smoothing_parameter=700
        elif(name=="viscous_DWC_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0512_p7_small_dt"):
            smoothing_parameter=700 # copied
        elif(name=="viscous_DWC_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs1024_p7"):
            smoothing_parameter=1500
        elif(name=="viscous_DWC_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0192_p2_stretched_mesh"):
            smoothing_parameter=1500
        else:
            print("ERROR: Invalid subdirectories_for_plot entry for smoothing parameter. Aborting...")
            exit()
        smoothing_parameters_store.append(smoothing_parameter)
    return smoothing_parameters_store
#=====================================================
def plot_for_presentation(
    subdirectories_for_plot,
    labels_for_plot,
    black_line_flag_for_plot,
    dashed_line_flag_for_plot,
    final_time_for_plot=1.0,
    legend_fontSize_input=14,
    plot_filtered_dns_input=False,
    plot_zoomed_section_dissipation_rate=False,
    plot_zoomed_section_numerical_dissipation_components=False,
    plot_zoomed_section_enstrophy=False,
    clr_input=[],mrkr_input=[],lnstl_input=[],
    solid_and_dashed_lines=False,
    dashed_and_solid_lines=False,
    dofs_for_zoomed_section=256,
    smoothing_parameters_input=[]):
    
    # flag to generate only the final plot with all the curves to save plotting time
    generate_only_final_plot_with_all_curves=True

    global subdirectories, filenames, labels, black_line_flag, \
    dashed_line_flag, figure_filename_postfix, figure_title, \
    ylimits_kinetic_energy_input, ylimits_dissipation_input, \
    log_axes_input, legend_on_input, legend_inside_input, \
    plot_reference_result, nlegendcols_input, \
    figure_subdirectory, data_directory_base, figure_directory_base, \
    smoothing_input, plot_PHiLiP_DNS_result_as_reference_input
    #-----------------------------------------------------
    number_of_result_curves=len(subdirectories_for_plot)
    for i in range(0,number_of_result_curves):
        figure_filename_postfix_input=figure_filename_postfix
        if(generate_only_final_plot_with_all_curves==False):
            figure_filename_postfix_input+="_%i" % i
        #-----------------------------------------------------
        subdirectories.append(subdirectories_for_plot[i])
        labels.append(labels_for_plot[i])
        black_line_flag.append(black_line_flag_for_plot[i])
        dashed_line_flag.append(dashed_line_flag_for_plot[i])
        filenames.append("turbulent_quantities.txt")
        #-----------------------------------------------------
        if(generate_only_final_plot_with_all_curves==False or i==(number_of_result_curves-1)):
            plot_periodic_turbulence(
                figure_subdirectory,
                subdirectories,
                filenames,
                labels,
                black_line_flag,
                dashed_line_flag,
                figure_directory_base,
                data_directory_base,
                True, # plot_reference_result,
                figure_filename_postfix_input,
                figure_title,
                log_axes_input,
                legend_on_input,
                legend_inside_input,
                nlegendcols_input,
                solid_and_dashed_lines=solid_and_dashed_lines,
                dashed_and_solid_lines=dashed_and_solid_lines,
                clr_input=clr_input,mrkr_input=mrkr_input,lnstl_input=lnstl_input,
                transparent_legend_input=True,
                tmax=final_time_for_plot,
                legend_fontSize_input=legend_fontSize_input,
                plot_kinetic_energy=True,
                plot_enstrophy=True,
                plot_palinstrophy=True,
                plot_PHiLiP_DNS_result_as_reference=plot_PHiLiP_DNS_result_as_reference_input,
                palinstrophy_smoothing=smoothing_input,
                plot_filtered_dns=plot_filtered_dns_input,
                plot_zoomed_section_dissipation_rate=plot_zoomed_section_dissipation_rate,
                plot_zoomed_section_numerical_dissipation_components=plot_zoomed_section_numerical_dissipation_components,
                plot_zoomed_section_enstrophy=plot_zoomed_section_enstrophy,
                dofs_for_zoomed_section=dofs_for_zoomed_section,
                check_smoothing_parameters=True,
                smoothing_parameters=smoothing_parameters_input)
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
    smoothing_input = True # set as [] for no smoothing
    plot_PHiLiP_DNS_result_as_reference_input=True # default
#=====================================================
#-----------------------------------------------------

#=====================================================
# DOFs: 256^3 | All results
#-----------------------------------------------------
if(False):
    #-----------------------------------------------------
    # clr_input = ['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2024_CSME/dipole_wall_collision/"
    date_for_runs="."
    figure_subdirectory="2024_CSME"
    # figure_title = "TGV at Re$_{\\infty}=1600$, $256^{3}$ DOFs, CFL=$0.10$" # comment to turn off
    figure_filename_postfix = "first_result"
    legend_inside_input=True
    plot_reference_result=True
    plot_PHiLiP_DNS_result_as_reference_input=False
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "viscous_DWC_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0192_p2",\
    "viscous_DWC_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0384_p5",\
    "viscous_DWC_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0512_p7",\
    "viscous_DWC_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0512_p7_small_dt",\
    "viscous_DWC_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs1024_p7",\
    "viscous_DWC_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0192_p2_small_dt",\
    "viscous_DWC_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0192_p2_stretched_mesh",\
    ]
    # labels
    labels_for_plot=[\
    # "$c_{DG}$ NSFR.IR-GL 64p$2$ ($192^2$ DOF)",\
    # "$c_{DG}$ NSFR.IR-GL 64p$5$ ($384^2$ DOF)",\
    # "$c_{DG}$ NSFR.IR-GL 64p$7$ ($512^2$ DOF)",\
    # "$c_{DG}$ NSFR.IR-GL 64p$2$ ($192^2$ DOF) $\\Delta t/2$",\
    "$c_{DG}$ NSFR 64p$2$ ($192^2$ DOF)\n $\\Delta t=10e^{-6}$",\
    "$c_{DG}$ NSFR 64p$5$ ($384^2$ DOF)\n $\\Delta t=10e^{-6}$",\
    "$c_{DG}$ NSFR 64p$7$ ($512^2$ DOF)\n $\\Delta t=5e^{-6}$",\
    "$c_{DG}$ NSFR 64p$7$ ($512^2$ DOF)\n $\\Delta t=2.5e^{-6}$",\
    "$c_{DG}$ NSFR 128p$7$ ($1024^2$ DOF)\n $\\Delta t=2.5e^{-6}$",\
    "$c_{DG}$ NSFR 64p$2$ ($192^2$ DOF)\n $\\Delta t=5e^{-6}$",\
    "$c_{DG}$ NSFR 64p$2$ ($192^2$ DOF)\n $\\Delta t=2.5e^{-6}$, Stretched mesh",\

    ]
    black_line_flag_for_plot=[False,False,False,False,False,False,False,False]
    dashed_line_flag_for_plot=[False,False,False,True,False,True,True]
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot,
        final_time_for_plot=1.0,legend_fontSize_input=12)

#=====================================================
# DOFs: 256^3 | All results
#-----------------------------------------------------
if(True):
    #-----------------------------------------------------
    # clr_input = ['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2024_CSME/dipole_wall_collision/"
    date_for_runs="."
    figure_subdirectory="2024_CSME"
    # figure_title = "TGV at Re$_{\\infty}=1600$, $256^{3}$ DOFs, CFL=$0.10$" # comment to turn off
    if(smoothing_input==True):
        figure_filename_postfix = "first_result_smoothed"
    else:
        figure_filename_postfix = "first_result"
    legend_inside_input=True
    plot_reference_result=True
    plot_PHiLiP_DNS_result_as_reference_input=False
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "viscous_DWC_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0192_p2",\
    "viscous_DWC_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0384_p5",\
    "viscous_DWC_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0512_p7",\
    "viscous_DWC_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0192_p2_stretched_mesh",\
    "viscous_DWC_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs1024_p7",\
    ]
    # labels
    labels_for_plot=[\
    "$c_{DG}$ NSFR 64p$2$ ($192^2$ DOF)\n $\\Delta t=10e^{-6}$",\
    "$c_{DG}$ NSFR 64p$5$ ($384^2$ DOF)\n $\\Delta t=10e^{-6}$",\
    "$c_{DG}$ NSFR 64p$7$ ($512^2$ DOF)\n $\\Delta t=5e^{-6}$",\
    "$c_{DG}$ NSFR 64p$2$ ($192^2$ DOF)\n $\\Delta t=2.5e^{-6}$, Stretched mesh",\
    "$c_{DG}$ NSFR 128p$7$ ($1024^2$ DOF)\n $\\Delta t=2.5e^{-6}$",\
    ]
    black_line_flag_for_plot=[False,False,False,False,False,False,False,False]
    dashed_line_flag_for_plot=[False,False,False,False,False,True,True]
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot,
        final_time_for_plot=1.0,legend_fontSize_input=12,
        smoothing_parameters_input=get_smoothing_parameters_from_subdirectories(subdirectories_for_plot))#smoothing_parameters_input=[400,750,1000]

#=====================================================
# DOFs: 256^3 | All results
#-----------------------------------------------------
if(False):
    #-----------------------------------------------------
    # clr_input = ['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2024_CSME/dipole_wall_collision/"
    date_for_runs="."
    figure_subdirectory="2024_CSME"
    # figure_title = "TGV at Re$_{\\infty}=1600$, $256^{3}$ DOFs, CFL=$0.10$" # comment to turn off
    figure_filename_postfix = "time_steps"
    legend_inside_input=True
    plot_reference_result=True
    plot_PHiLiP_DNS_result_as_reference_input=False
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "viscous_DWC_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0192_p2",\
    "viscous_DWC_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0192_p2_small_dt",\
    "viscous_DWC_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0512_p7",\
    "viscous_DWC_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0512_p7_small_dt",\
    ]
    # labels
    labels_for_plot=[\
    # "$c_{DG}$ NSFR.IR-GL 64p$2$ ($192^2$ DOF)",\
    # "$c_{DG}$ NSFR.IR-GL 64p$5$ ($384^2$ DOF)",\
    # "$c_{DG}$ NSFR.IR-GL 64p$7$ ($512^2$ DOF)",\
    # "$c_{DG}$ NSFR.IR-GL 64p$2$ ($192^2$ DOF) $\\Delta t/2$",\
    "$c_{DG}$ NSFR 64p$2$ ($192^2$ DOF)\n $\\Delta t=10e^{-6}$",\
    "$c_{DG}$ NSFR 64p$2$ ($192^2$ DOF)\n $\\Delta t=5e^{-6}$",\
    "$c_{DG}$ NSFR 64p$7$ ($512^2$ DOF)\n $\\Delta t=5e^{-6}$",\
    "$c_{DG}$ NSFR 64p$7$ ($512^2$ DOF)\n $\\Delta t=2.5e^{-6}$",\
    ]
    black_line_flag_for_plot=[False,False,False,False,False,False,False,False]
    dashed_line_flag_for_plot=[False,False,False,False,False,True,True]
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot,
        final_time_for_plot=1.0,legend_fontSize_input=12)
