import os;CURRENT_PATH = os.path.split(os.path.realpath(__file__))[0]+"/";
import sys
sys.path.append(CURRENT_PATH+"../../src");
from plot_unsteady_integrated_turbulent_flow_quantities import plot_periodic_turbulence

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
def plot_for_presentation(
    subdirectories_for_plot,
    labels_for_plot,
    black_line_flag_for_plot,
    dashed_line_flag_for_plot,
    final_time_for_plot=10.0,
    legend_fontSize_input=14,
    plot_filtered_dns_input=False,
    plot_zoomed_section_dissipation_rate=False,
    plot_zoomed_section_numerical_dissipation_components=False,
    plot_zoomed_section_enstrophy=False):
    
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
                plot_reference_result,
                figure_filename_postfix_input,
                figure_title,
                log_axes_input,
                legend_on_input,
                legend_inside_input,
                nlegendcols_input,
                # clr_input=clr_input,
                plot_numerical_viscosity=True,
                plot_dissipation_components=True,
                transparent_legend_input=True,
                tmax=final_time_for_plot,
                legend_fontSize_input=legend_fontSize_input,
                solid_and_dashed_lines=False,
                plot_kinetic_energy=True,
                plot_enstrophy=True,
                plot_numerical_dissipation=True,
                plot_PHiLiP_DNS_result_as_reference=plot_PHiLiP_DNS_result_as_reference_input,
                dissipation_rate_smoothing=smoothing_input,
                plot_filtered_dns=plot_filtered_dns_input,
                plot_zoomed_section_dissipation_rate=plot_zoomed_section_dissipation_rate,
                plot_zoomed_section_numerical_dissipation_components=plot_zoomed_section_numerical_dissipation_components,
                plot_zoomed_section_enstrophy=plot_zoomed_section_enstrophy)
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
    smoothing_input = []
    plot_PHiLiP_DNS_result_as_reference_input=True # default
#=====================================================
#-----------------------------------------------------
#=====================================================
# DOFs: 96^3 | Correction Parameter Accuracy
#-----------------------------------------------------
if(True):
    #-----------------------------------------------------
    # clr_input = ['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    # figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL=$0.10$" # comment to turn off
    figure_filename_postfix = "96_p5_correction_parameter"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512",\
    "correction_parameter/viscous_TGV_ILES_NSFR_cSD_IR_2PF_GL_OI-0_dofs096_p5_procs512",\
    "correction_parameter/viscous_TGV_ILES_NSFR_cHU_IR_2PF_GL_OI-0_dofs096_p5_procs512",\
    "correction_parameter/viscous_TGV_ILES_NSFR_cPlus_IR_2PF_GL_OI-0_dofs096_p5_procs512",\
    ]
    # labels
    labels_for_plot=[\
    "$c_{DG}$",\
    "$c_{SD}$",\
    "$c_{HU}$",\
    "$c_{+}$",\
    ]
    black_line_flag_for_plot=[False,False,False,False]
    dashed_line_flag_for_plot=[False,False,False,False]
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot,
        final_time_for_plot=20.0,plot_filtered_dns_input=True,legend_fontSize_input=12,
        plot_zoomed_section_dissipation_rate=True,
        plot_zoomed_section_numerical_dissipation_components=False,
        plot_zoomed_section_enstrophy=False)
#=====================================================
# DOFs: 96^3 | Correction Parameter Time-Step
#-----------------------------------------------------
if(True):
    #-----------------------------------------------------
    # clr_input = ['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    # figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs" # comment to turn off
    figure_filename_postfix = "96_p5_correction_parameter_cfl_advantage"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512",\
    "time_step_advantage_with_physical_check/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.26_procs512",\
    "correction_parameter/viscous_TGV_ILES_NSFR_cPlus_IR_2PF_GL_OI-0_dofs096_p5_procs512",\
    "time_step_advantage_with_physical_check/viscous_TGV_ILES_NSFR_cPlus_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.36_procs512",\
    "flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512",\
    "time_step_advantage_strong_DG/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_CFL-0.14_procs512",\
    ] 
    # labels
    labels_for_plot=[\
    "$c_{DG}$: CFL=$0.10$",\
    "$c_{DG}$: CFL=$0.26$",\
    "$c_{+}$: CFL=$0.10$",\
    "$c_{+}$: CFL=$0.36$",\
    "sDG: CFL=$0.10$",\
    "sDG: CFL=$0.14$",\
    ]
    black_line_flag_for_plot=[False,False,False,False,False,False]
    dashed_line_flag_for_plot=[False,True,False,True,False,True]
    smoothing_input = [True,True,True,True,True,True]
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot,
        final_time_for_plot=20.0,plot_filtered_dns_input=True,legend_fontSize_input=12,
        plot_zoomed_section_dissipation_rate=True,
        plot_zoomed_section_numerical_dissipation_components=False,
        plot_zoomed_section_enstrophy=False)

#=====================================================
# DOFs: 256^3 | All results
#-----------------------------------------------------
if(True):
    #-----------------------------------------------------
    # clr_input = ['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    # figure_title = "TGV at Re$_{\\infty}=1600$, $256^{3}$ DOFs, CFL=$0.10$" # comment to turn off
    figure_filename_postfix = "256_verification"
    legend_inside_input=True
    plot_reference_result=True
    plot_PHiLiP_DNS_result_as_reference_input=False
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "filtered_dns_viscous_tgv/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p7_procs1024",\
    "verification/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p3_procs1024",\
    "verification/viscous_TGV_ILES_NSFR_cDG_IR_2PF-Roe_GL_OI-0_dofs0256_p3_procs1024",\
    "verification/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-4_dofs0256_p3_CFL-0.15_procs1024",\
    ]
    # labels
    labels_for_plot=[\
    "$256^3$ P$7$ $c_{DG}$ NSFR.IR-GL",\
    "$256^3$ P$3$ $c_{DG}$ NSFR.IR-GL",\
    "$256^3$ P$3$ $c_{DG}$ NSFR.IR-GL-Roe",\
    "$256^3$ P$3$ Strong DG-Roe-GL-OI",\
    ]
    black_line_flag_for_plot=[False,False,False,False]
    dashed_line_flag_for_plot=[False,False,False,True]
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot,
        final_time_for_plot=20.0,plot_filtered_dns_input=False,legend_fontSize_input=12,
        plot_zoomed_section_dissipation_rate=True,
        plot_zoomed_section_numerical_dissipation_components=True,
        plot_zoomed_section_enstrophy=True)

#=====================================================
# DOFs: 96^3 | NSFR vs Strong DG on GL and GLL flux nodes
#-----------------------------------------------------
if(True):
    # clr_input = ['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    # figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL$=0.10$" # comment to turn off
    figure_filename_postfix = "96_p5_flux_nodes_with_sDG"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories_for_plot = [\
    "flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512", \
    "flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs096_p5_procs512", \
    "flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512", \
    "flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GLL_OI-6_dofs096_p5_procs512", \
    ]
    labels_for_plot=[\
    "$c_{DG}$ NSFR.IR-GL", \
    "$c_{DG}$ NSFR.IR-GLL", \
    "Strong DG-Roe-GL-OI", \
    "Strong DG-Roe-GLL-OI", \
    ]
    black_line_flag_for_plot=[False,False,False,False,False,False]
    dashed_line_flag_for_plot=[False,False,False,True,False,True]
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot,
        final_time_for_plot=20.0,plot_filtered_dns_input=True,legend_fontSize_input=12,
        plot_zoomed_section_dissipation_rate=True,
        plot_zoomed_section_numerical_dissipation_components=False,
        plot_zoomed_section_enstrophy=False)
exit()

#=====================================================
# DOFs: 96^3 | Strong DG GLL without OI, SGS model stabilization
#-----------------------------------------------------
if(True):
    # clr_input = ['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    # figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL$=0.10$" # comment to turn off
    figure_filename_postfix = "96_sDG_gll_sgs_model_stabilization"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories_for_plot = [\
    "flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GLL_OI-6_dofs096_p5_procs512", \
    "filter_width_stabilization/viscous_TGV_ILES_std_strong_DG_Roe_GLL_OI-0_dofs096_p5_procs512", \
    "filter_width_stabilization/viscous_TGV_LES_SMAG_MC-0.18_std_strong_DG_Roe_GLL_OI-0_dofs096_p5_procs512", \
    ]
    labels_for_plot=[\
    "Strong DG-Roe-GLL-OI", \
    "Strong DG-Roe-GLL", \
    "Strong DG-Roe-GLL-Smag. $C_{S}=0.18$", \
    ]
    black_line_flag_for_plot=[False,False,False,False,False,False]
    dashed_line_flag_for_plot=[False,True,False,False,False]
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot,final_time_for_plot=20.0,plot_filtered_dns_input=True,legend_fontSize_input=12)

#=====================================================
# DOFs: 64^3 | Over-integration stabilization
#-----------------------------------------------------
if(True):
    # clr_input = ['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    # figure_title = "TGV at Re$_{\\infty}=1600$, P$7$, $64^{3}$ DOFs, CFL$=0.10$" # comment to turn off
    figure_filename_postfix = "64_p7_overintegration_stability"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories_for_plot = [\
    "high_poly_degree_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs064_p7_procs512",\
    "high_poly_degree_GL_flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-8_dofs064_p7_procs512",\
    "high_poly_degree_GL_flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-4_dofs064_p7_CFL-0.10_procs512",\
    "high_poly_degree_GL_flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-2_dofs064_p7_CFL-0.10_procs512",\
    "high_poly_degree_GL_flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-1_dofs064_p7_CFL-0.10_procs512",\
    "high_poly_degree_GL_flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-0_dofs064_p7_CFL-0.10_procs512",\
    ]
    labels_for_plot=[\
    "$c_{DG}$ NSFR.IR-GL",\
    "sDG-OI.8",\
    "sDG-OI.4",\
    "sDG-OI.2",\
    "sDG-OI.1",\
    "sDG-OI.0",\
    ]
    black_line_flag_for_plot=[False,False,False,False,False,False]
    dashed_line_flag_for_plot=[False,False,True,False,False,True]
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot,final_time_for_plot=20.0,plot_filtered_dns_input=False,legend_fontSize_input=12)

#=====================================================
# DOFs: 48,64,96 | Strong DG vs NSFR CONVERGENCE
#-----------------------------------------------------
if(True):
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    # figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL$=0.10$" # comment to turn off
    figure_filename_postfix = "convergence_comparison"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512",\
    "flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512",\
    "high_poly_degree_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs064_p7_procs512",\
    "high_poly_degree_GL_flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-8_dofs064_p7_procs512",\
    "robustness/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs048_p5_procs64",\
    "robustness/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs048_p5_procs64",\
    # "robustness/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs024_p5_procs16",\
    # "robustness/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs024_p5_procs16",\
    # "robustness/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs012_p5_procs16",\
    # "robustness/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs012_p5_procs16",\
    ]
    # labels
    labels_for_plot=[\
    "$96^{3}$, P$5$, $c_{DG}$", \
    "$96^{3}$, P$5$, sDG", \
    "$64^{3}$, P$7$, $c_{DG}$", \
    "$64^{3}$, P$7$, sDG", \
    "$48^{3}$, P$5$, $c_{DG}$",\
    "$48^{3}$, P$5$, sDG",\
    # "$24^{3}$, P$5$",\
    # "$24^{3}$, P$5$",\
    # "$12^{3}$, P$5$",\
    # "$12^{3}$, P$5$",\
    ]
    black_line_flag_for_plot=[False,False,False,False,False,False,False]
    dashed_line_flag_for_plot=[False,False,False,False,False,False,False]
    smoothing_input = [False,False,False,False,False,False,False]
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot,final_time_for_plot=20.0,plot_filtered_dns_input=True,legend_fontSize_input=12)

#=====================================================
# DOFs: 96^3 | Basic SGS Models on GLL flux nodes (no filter width modifications)
#-----------------------------------------------------
if(False):
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    # figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL$=0.10$" # comment to turn off
    figure_filename_postfix = "96_sgs_models_gll"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs096_p5_procs512",\
    "sgs_model_GLL_flux_nodes/viscous_TGV_LES_SMAG_MC-0.10_NSFR_cDG_IR_2PF_GLL_OI-0_dofs096_p5_procs512",\
    "sgs_model_GLL_flux_nodes/viscous_TGV_LES_SMAG_MC-0.18_NSFR_cDG_IR_2PF_GLL_OI-0_dofs096_p5_procs512",\
    "sgs_model_GLL_flux_nodes/viscous_TGV_LES_WALE_MC-0.50_NSFR_cDG_IR_2PF_GLL_OI-0_dofs096_p5_procs512",\
    "sgs_model_GLL_flux_nodes/viscous_TGV_LES_VRMN_MC-0.081_NSFR_cDG_IR_2PF_GLL_OI-0_dofs096_p5_procs512",\
    ]
    # labels
    labels_for_plot=[\
    "$c_{DG}$ NSFR.IR-GLL", \
    "$c_{DG}$ NSFR.IR-GLL-Smag. $C_{S}=0.10$", \
    "$c_{DG}$ NSFR.IR-GLL-Smag. $C_{S}=0.18$",\
    "$c_{DG}$ NSFR.IR-GLL-WALE $C_{W}=0.50$",\
    "$c_{DG}$ NSFR.IR-GLL-VRMN $C_{V}=0.081$",\
    ]
    black_line_flag_for_plot=[False,False,False,False,False,False,False]
    dashed_line_flag_for_plot=[False,False,False,False,False,False,False]
    smoothing_input = [False,False,False,False,False,False,False]
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot,final_time_for_plot=20.0,plot_filtered_dns_input=True,legend_fontSize_input=12)

#=====================================================
# DOFs: 96^3 | sDG with and without OI vs NSFR
#-----------------------------------------------------
if(True):
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    # figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL$=0.10$" # comment to turn off
    figure_filename_postfix = "96_sDG_with_and_without_OI_vs_NSFR"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512",\
    "filter_width_stabilization/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-0_dofs096_p5_procs512",\
    "flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512",\
    "upwind_dissipation_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF-Roe_GL_OI-0_dofs096_p5_procs512",\
    ]
    # labels
    labels_for_plot=[\
    "Strong DG-Roe-GL-OI", \
    "Strong DG-Roe-GL", \
    "$c_{DG}$ NSFR.IR-GL",\
    "$c_{DG}$ NSFR.IR-GL-Roe",\
    ]
    black_line_flag_for_plot=[False,False,False,False,False,False,False]
    dashed_line_flag_for_plot=[False,False,False,False,False,False,False]
    smoothing_input = [False,False,False,False,False,False,False]
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot,final_time_for_plot=20.0,plot_filtered_dns_input=True,legend_fontSize_input=12)

#=====================================================
# DOFs: 96^3 | OI Stability on GL flux nodes
#-----------------------------------------------------
if(True):
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    # figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL$=0.10$" # comment to turn off
    figure_filename_postfix = "96_OI_stability_GL"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    # "upwind_dissipation_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF-Roe_GL_OI-0_dofs096_p5_procs512",\
    "flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512",\
    "flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512",\
    "over_integration_accuracy_strong_DG/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-4_dofs096_p5_CFL-0.10_procs512",\
    "over_integration_accuracy_strong_DG/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-2_dofs096_p5_CFL-0.10_procs512",\
    "filter_width_stabilization/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-0_dofs096_p5_procs512",\
    ]
    # labels
    labels_for_plot=[\
    # "$c_{DG}$ NSFR.IR-GL-Roe",\
    "$c_{DG}$ NSFR.IR-GL",\
    "sDG-OI-6", \
    "sDG-OI-4", \
    "sDG-OI-2",\
    "sDG-OI-0",\
    ]
    black_line_flag_for_plot=[False,False,False,False,False,False,False]
    dashed_line_flag_for_plot=[False,False,False,False,False,False,False]
    smoothing_input = [False,False,False,False,False,False,False]
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot,final_time_for_plot=20.0,plot_filtered_dns_input=True,legend_fontSize_input=12)
exit()
#=====================================================
# DOFs: ALL | Strong DG CONVERGENCE
#-----------------------------------------------------
if(True):
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    # figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL$=0.10$" # comment to turn off
    figure_filename_postfix = "strong_DG_convergence"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512",\
    "high_poly_degree_GL_flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-8_dofs064_p7_procs512",\
    "robustness/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs048_p5_procs64",\
    "robustness/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs024_p5_procs16",\
    # "robustness/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs012_p5_procs16",\
    ]
    # labels
    labels_for_plot=[\
    "$96^{3}$ DOFs, P$5$", \
    "$64^{3}$ DOFs, P$7$", \
    "$48^{3}$ DOFs, P$5$",\
    "$24^{3}$ DOFs, P$5$",\
    # "$12^{3}$ DOFs, P$5$",\
    ]
    black_line_flag_for_plot=[False,False,False,False,False,False,False]
    dashed_line_flag_for_plot=[False,False,False,False,False,False,False]
    smoothing_input = [False,False,False,False,False,False,False]
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot,final_time_for_plot=20.0,plot_filtered_dns_input=True,legend_fontSize_input=12)

#=====================================================
# DOFs: 96^3 | Basic SGS Models on GL flux nodes (no filter width modifications)
#-----------------------------------------------------
if(True):
    # batch_labels = [ \
    # "$c_{DG}$ NSFR", \
    # "$c_{DG}$ NSFR.IR-Smag. $C_{S}=0.10$", \
    # "$c_{DG}$ NSFR.IR-Smag. $C_{S}=0.18$", \
    # "$c_{DG}$ NSFR.IR-WALE", \
    # "$c_{DG}$ NSFR.IR-VRMN", \
    # "$c_{DG}$ NSFR.IR-Smag. $C_{S}=0.1$, $\\Delta_{min}$", \
    # ]
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    # figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL$=0.10$" # comment to turn off
    figure_filename_postfix = "96_p5_sgs_models_gl"
    legend_inside_input=True
    # #-----------------------------------------------------
    # subdirectories.append("sgs_model_GL_flux_nodes/viscous_TGV_LES_SMAG_MC-0.18_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512_filter_width_flad_and_gassner")
    # filenames.append("turbulent_quantities.txt")
    # labels.append("$c_{DG}$ NSFR.IR-GL-Smag. $C_{S}=0.18$ $\\Delta=\\frac{V}{(P+1)^{3}}$")
    # black_line_flag.append(False)
    # dashed_line_flag.append(True)
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512",\
    "sgs_model_GL_flux_nodes/viscous_TGV_LES_SMAG_MC-0.18_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512",\
    "sgs_model_GL_flux_nodes/viscous_TGV_LES_SMAG_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512",\
    "sgs_model_GL_flux_nodes/viscous_TGV_LES_WALE_MC-0.50_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512",\
    "sgs_model_GL_flux_nodes/viscous_TGV_LES_VRMN_MC-0.081_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512",\
    "sgs_model_GL_flux_nodes/viscous_TGV_LES_SI.SMAG_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs512",\
    ]
    # labels
    labels_for_plot=[\
    "$c_{DG}$ NSFR.IR-GL (No Model)", \
    "Smag. $C_{S}=0.18$", \
    "Smag. $C_{S}=0.10$", \
    "WALE $C_{W}=0.50$", \
    "VRMN $C_{V}=0.081$", \
    "SI.Smag. $C_{S}=0.10$", \
    ]
    black_line_flag_for_plot=[False,False,False,False,False,False]
    dashed_line_flag_for_plot=[False,False,False,False,False,True]
    smoothing_input=[False,False,False,False,False,False]
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot,final_time_for_plot=20.0,plot_filtered_dns_input=True,legend_fontSize_input=12)

#=====================================================
# DOFs: 96^3 | LRNC Advanced SGS Models on GL flux nodes (no filter width modifications)
#-----------------------------------------------------
if(True):
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    # figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL$=0.10$" # comment to turn off
    figure_filename_postfix = "96_p5_lrnc_sgs_models_gl"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512",\
    "sgs_model_GL_flux_nodes/viscous_TGV_LES_SMAG_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512",\
    "sgs_model_GL_flux_nodes/viscous_TGV_LES_WALE_MC-0.50_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512",\
    "sgs_model_GL_flux_nodes/viscous_TGV_LES_VRMN_MC-0.081_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512",\
    "sgs_model_GL_flux_nodes/viscous_TGV_LES_SMAG.LRNC_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs512",\
    "sgs_model_GL_flux_nodes/viscous_TGV_LES_WALE.LRNC_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs16",\
    "sgs_model_GL_flux_nodes/viscous_TGV_LES_VRMN.LRNC_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs16",\
    ]
    # labels
    labels_for_plot=[\
    "$c_{DG}$ NSFR.IR-GL (No Model)", \
    "Smag. $C_{S}=0.10$", \
    "WALE $C_{W}=0.50$", \
    "VRMN $C_{V}=0.081$", \
    "Smag.LRNC $C_{S}=0.10$", \
    "WALE.LRNC $C_{W}=0.10$", \
    "VRMN.LRNC $C_{V}=0.10$", \
    ]
    black_line_flag_for_plot=[False,False,False,False,False,False,False]
    dashed_line_flag_for_plot=[False,False,False,False,True,True,True]
    smoothing_input=[False,False,False,False,True,True,True]
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot,final_time_for_plot=20.0,plot_filtered_dns_input=True,legend_fontSize_input=12)

#=====================================================
# DOFs: 96^3 | LRNC Advanced SGS Models on GL flux nodes (no filter width modifications)
#-----------------------------------------------------
if(True):
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    # figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL$=0.10$" # comment to turn off
    figure_filename_postfix = "96_p5_lrnc_advanced_sgs_models_gl"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512",\
    "sgs_model_GL_flux_nodes/viscous_TGV_LES_SMAG_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512",\
    "sgs_model_GL_flux_nodes/viscous_TGV_LES_SMAG.LRNC_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs512",\
    "sgs_model_GL_flux_nodes/viscous_TGV_LES_SI.SMAG.LRNC_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs16",\
    "sgs_model_GL_flux_nodes/viscous_TGV_LES_filtered_pL3_SMAG.LRNC_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs16",\
    "sgs_model_GL_flux_nodes/viscous_TGV_LES_filtered_pL3_SI.SMAG.LRNC_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs16",\
    "sgs_model_GL_flux_nodes/viscous_TGV_LES_DYNAMIC.SMAG.LRNC_CLIPMC-0.01-pL3_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs16",\
    ]
    # labels
    labels_for_plot=[\
    "$c_{DG}$ NSFR.IR-GL (No Model)", \
    "Smag. $C_{S}=0.10$", \
    "Smag.LRNC $C_{S}=0.10$", \
    "SI.Smag.LRNC $C_{S}=0.10$", \
    "HPF.Smag.LRNC $C_{S}=0.10$ $P_{L}=3$", \
    "HPF.SI.Smag.LRNC $C_{S}=0.10$ $P_{L}=3$", \
    "Dyn.Smag.LRNC ($P_{TF}=3$, $C_{max}=0.1$)", \
    ]
    black_line_flag_for_plot=[False,False,False,False,False,False,False,False,False,False]
    dashed_line_flag_for_plot=[False,False,False,False,False,False,False,False,False,False]
    smoothing_input=[False,False,False,False,False,False,False,False,False,False]
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot,final_time_for_plot=20.0,plot_filtered_dns_input=True,legend_fontSize_input=12)

#=====================================================
# DOFs: 96^3 | Advanced SGS Models on GL flux nodes (no filter width modifications)
#-----------------------------------------------------
if(True):
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    # figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL$=0.10$" # comment to turn off
    figure_filename_postfix = "96_p5_advanced_sgs_models_gl"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512",\
    "sgs_model_GL_flux_nodes/viscous_TGV_LES_SMAG_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512",\
    # "sgs_model_GL_flux_nodes/viscous_TGV_LES_filtered_pL2_SMAG_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs512",\
    "sgs_model_GL_flux_nodes/viscous_TGV_LES_filtered_pL3_SMAG_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs512",\
    # "sgs_model_GL_flux_nodes/viscous_TGV_LES_filtered_pL4_SMAG_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs512",\
    "sgs_model_GL_flux_nodes/viscous_TGV_LES_SI.SMAG_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs512",\
    # "sgs_model_GL_flux_nodes/viscous_TGV_LES_filtered_pL3_SS.VMS_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs512",\
    "sgs_model_GL_flux_nodes/viscous_TGV_LES_filtered_pL3_SI.SMAG_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs512",\
    "sgs_model_GL_flux_nodes/viscous_TGV_LES_DYNAMIC.SMAG-pL3_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs512",\
    "sgs_model_GL_flux_nodes/viscous_TGV_LES_SMAG.LRNC_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.1_procs512",\
    ]
    # labels
    labels_for_plot=[\
    "$c_{DG}$ NSFR.IR-GL (No Model)", \
    "Smag. $C_{S}=0.10$", \
    #HPF.Smag. $C_{S}=0.10$ $P_{L}=2$", \
    "HPF.Smag. $C_{S}=0.10$ $P_{L}=3$", \
    #HPF.Smag. $C_{S}=0.10$ $P_{L}=4$", \
    "SI.Smag. $C_{S}=0.10$", \
    # "$SS.VMS $P_{L}=3$", \
    "HPF.SI.Smag. $C_{S}=0.10$ $P_{L}=3$", \
    "Dyn.Smag. ($P_{TF}=3$, $C_{max}=0.1$)", \
    "Smag.LRNC $C_{S}=0.10$", \
    ]
    black_line_flag_for_plot=[False,False,False,False,False,False,False]
    dashed_line_flag_for_plot=[False,False,False,False,False,False,True]
    smoothing_input=[False,False,False,False,False,False,False]
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot,final_time_for_plot=20.0,plot_filtered_dns_input=True,legend_fontSize_input=12)

#=====================================================
# DOFs: ALL | NSFR CONVERGENCE
#-----------------------------------------------------
if(True):
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    # figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL$=0.10$" # comment to turn off
    figure_filename_postfix = "cDG_NSFR_convergence"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512",\
    "high_poly_degree_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs064_p7_procs512",\
    "robustness/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs048_p5_procs64",\
    "robustness/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs024_p5_procs16",\
    # "robustness/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs012_p5_procs16",\
    ]
    # labels
    labels_for_plot=[\
    "$96^{3}$ DOFs, P$5$", \
    "$64^{3}$ DOFs, P$7$", \
    "$48^{3}$ DOFs, P$5$",\
    "$24^{3}$ DOFs, P$5$",\
    # "$12^{3}$ DOFs, P$5$",\
    ]
    black_line_flag_for_plot=[False,False,False,False,False,False,False]
    dashed_line_flag_for_plot=[False,False,False,False,False,False,False]
    smoothing_input = [False,False,False,False,False,False,False]
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot,final_time_for_plot=20.0,plot_filtered_dns_input=True,legend_fontSize_input=12)

#=====================================================
# DOFs: ALL | Strong DG CONVERGENCE
#-----------------------------------------------------
if(False):
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, CFL$=0.10$" # comment to turn off
    figure_filename_postfix = "strong_DG_convergence"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories.append("flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$96^{3}$ DOFs, Strong DG-Roe-GL-OI")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("robustness/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs048_p5_procs64")
    filenames.append("turbulent_quantities.txt")
    labels.append("$48^{3}$ DOFs, Strong DG-Roe-GL-OI")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("robustness/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-0_dofs048_p5_procs64")
    filenames.append("turbulent_quantities.txt")
    labels.append("$48^{3}$ DOFs, Strong DG-Roe-GL")
    black_line_flag.append(False)
    dashed_line_flag.append(True)
    #-----------------------------------------------------
    subdirectories.append("robustness/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs024_p5_procs16")
    filenames.append("turbulent_quantities.txt")
    labels.append("$24^{3}$ DOFs, Strong DG-Roe-GL-OI")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("robustness/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-0_dofs024_p5_procs16")
    filenames.append("turbulent_quantities.txt")
    labels.append("$24^{3}$ DOFs, Strong DG-Roe-GL")
    black_line_flag.append(False)
    dashed_line_flag.append(True)
    #-----------------------------------------------------
    subdirectories.append("robustness/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs012_p5_procs16")
    filenames.append("turbulent_quantities.txt")
    labels.append("$12^{3}$ DOFs, Strong DG-Roe-GL-OI")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("robustness/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-0_dofs012_p5_procs16")
    filenames.append("turbulent_quantities.txt")
    labels.append("$12^{3}$ DOFs, Strong DG-Roe-GL")
    black_line_flag.append(False)
    dashed_line_flag.append(True)
    #-----------------------------------------------------
    plot_periodic_turbulence(
        figure_subdirectory,
        subdirectories,
        filenames,
        labels,
        black_line_flag,
        dashed_line_flag,
        figure_directory_base,
        data_directory_base,
        plot_reference_result,
        figure_filename_postfix,
        figure_title,
        log_axes_input,
        legend_on_input,
        legend_inside_input,
        nlegendcols_input,
        transparent_legend_input=True,
        tmax=20.0,
        legend_fontSize_input=14,
        solid_and_dashed_lines=False)
    #-----------------------------------------------------
#=====================================================
# DOFs: 96^3 | Collocated vs Uncollocated
#-----------------------------------------------------
if(True):
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    # figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL$=0.10$" # comment to turn off
    figure_filename_postfix = "96_p5_flux_nodes"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512",\
    "flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs096_p5_procs512",\
    ]
    # labels
    labels_for_plot=[\
    "$c_{DG}$ NSFR.IR-GL", \
    "$c_{DG}$ NSFR.IR-GLL",\
    ]
    black_line_flag_for_plot=[False,False,False,False]
    dashed_line_flag_for_plot=[False,False,False,False]
    smoothing_input = [False,False,False,False]
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot,final_time_for_plot=20.0,plot_filtered_dns_input=True,legend_fontSize_input=12)

#=====================================================
# DOFs: 96^3 | Upwind dissipation on GL flux nodes
#-----------------------------------------------------
if(True):
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    # figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL$=0.10$" # comment to turn off
    figure_filename_postfix = "96_p5_upwind_gl"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512",\
    "upwind_dissipation_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF-LxF_GL_OI-0_dofs096_p5_procs512",\
    "upwind_dissipation_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF-Roe_GL_OI-0_dofs096_p5_procs512",\
    "upwind_dissipation_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF-L2R_GL_OI-0_dofs096_p5_procs512",\
    ]
    # labels
    labels_for_plot=[\
    "$c_{DG}$ NSFR.IR-GL", \
    "$c_{DG}$ NSFR.IR-GL-LxF", \
    "$c_{DG}$ NSFR.IR-GL-Roe", \
    "$c_{DG}$ NSFR.IR-GL-L2R", \
    ]
    black_line_flag_for_plot=[False,False,False,False]
    dashed_line_flag_for_plot=[False,False,False,False]
    smoothing_input = [False,False,False,False]
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot,final_time_for_plot=20.0,plot_filtered_dns_input=True,legend_fontSize_input=12)

#=====================================================
# DOFs: 96^3 | NSFR Two-Point-Flux
#-----------------------------------------------------
if(True):
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    # figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL$=0.10$" # comment to turn off
    figure_filename_postfix = "96_p5_two_point_flux"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512",\
    "two_point_flux/viscous_TGV_ILES_NSFR_cDG_KG_2PF_GL_OI-0_dofs096_p5_procs512",\
    "two_point_flux/viscous_TGV_ILES_NSFR_cDG_CH_2PF_GL_OI-0_dofs096_p5_procs512",\
    "two_point_flux/viscous_TGV_ILES_NSFR_cDG_Ra_2PF_GL_OI-0_dofs096_p5_procs512",\
    ]
    # labels
    labels_for_plot=[\
    "$c_{DG}$ NSFR.IR-GL", \
    "$c_{DG}$ NSFR.KG-GL", \
    "$c_{DG}$ NSFR.CH-GL", \
    "$c_{DG}$ NSFR.CH$_{\\mathrm{RA}}$-GL", \
    ]
    black_line_flag_for_plot=[False,False,False,False]
    dashed_line_flag_for_plot=[False,True,True,True]
    smoothing_input = [False,False,False,False]
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot,final_time_for_plot=20.0,plot_filtered_dns_input=True,legend_fontSize_input=12)

exit()
#////////// everything after has not been updated
#=====================================================
# DOFs: 96^3 | cDG vs cPlus
#-----------------------------------------------------
if(True):
    #-----------------------------------------------------
    # clr_input = ['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL=$0.10$" # comment to turn off
    figure_filename_postfix = "96_p5_cDG_vs_cPlus"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512",\
    # "correction_parameter/viscous_TGV_ILES_NSFR_cSD_IR_2PF_GL_OI-0_dofs096_p5_procs512",\
    # "correction_parameter/viscous_TGV_ILES_NSFR_cHU_IR_2PF_GL_OI-0_dofs096_p5_procs512",\
    "correction_parameter/viscous_TGV_ILES_NSFR_cPlus_IR_2PF_GL_OI-0_dofs096_p5_procs512",\
    ]
    # labels
    labels_for_plot=[\
    "$c_{DG}$ NSFR.IR-GL",\
    # "$c_{SD}$ NSFR.IR-GL",\
    # "$c_{HU}$ NSFR.IR-GL",\
    "$c_{+}$ NSFR.IR-GL",\
    ]
    black_line_flag_for_plot=[False,False,False,False]
    dashed_line_flag_for_plot=[False,False,False,False]
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot,final_time_for_plot=20.0)

#=====================================================
# DOFs: 96^3 | Over-integration vs Split Form De-aliasing strategy
#-----------------------------------------------------
if(False):
    #-----------------------------------------------------
    # clr_input = ['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    reinit_inputs()
    # /Volumes/Samsung_T5/NarvalFiles
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL=$0.10$" # comment to turn off
    figure_filename_postfix = "96_p5_OI_vs_SF"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512",\
    "flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512",\
    ] 
    # labels
    labels_for_plot=[\
    "Strong DG-Roe-GL-OI",\
    "$c_{DG}$ NSFR.IR-GL",\
    ]
    black_line_flag_for_plot=[False,False]
    dashed_line_flag_for_plot=[False,False]
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot)

#=====================================================
# DOFs: 96^3 | Time step advantage with physical solution check
#-----------------------------------------------------
if(False):
    #-----------------------------------------------------
    # clr_input = ['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    reinit_inputs()
    # /Volumes/Samsung_T5/NarvalFiles
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs" # comment to turn off
    figure_filename_postfix = "96_p5_cPlus_cfl_physically_consistent"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories.append("time_step_advantage_with_physical_check/viscous_TGV_ILES_NSFR_cPlus_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.36_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{+}$ NSFR.IR-GL CFL=$0.36$")
    black_line_flag.append(False)
    dashed_line_flag.append(True)
    # #-----------------------------------------------------
    # subdirectories.append("time_step_advantage_with_physical_check/viscous_TGV_ILES_NSFR_cPlus_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.34_procs512")
    # filenames.append("turbulent_quantities.txt")
    # labels.append("$c_{+}$ NSFR.IR-GL CFL=$0.34$")
    # black_line_flag.append(False)
    # dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("time_step_advantage_with_physical_check/viscous_TGV_ILES_NSFR_cPlus_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.32_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{+}$ NSFR.IR-GL CFL=$0.32$")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("time_step_advantage/viscous_TGV_ILES_NSFR_cPlus_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.3_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{+}$ NSFR.IR-GL CFL=$0.30$")
    black_line_flag.append(False)
    dashed_line_flag.append(True)
    #-----------------------------------------------------
    subdirectories.append("time_step_advantage/viscous_TGV_ILES_NSFR_cPlus_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.2_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{+}$ NSFR.IR-GL CFL=$0.20$")
    black_line_flag.append(False)
    dashed_line_flag.append(True)
    #-----------------------------------------------------
    subdirectories.append("correction_parameter/viscous_TGV_ILES_NSFR_cPlus_IR_2PF_GL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{+}$ NSFR.IR-GL CFL=$0.10$")
    black_line_flag.append(False)
    dashed_line_flag.append(True)
    #-----------------------------------------------------
    plot_periodic_turbulence(
        figure_subdirectory,
        subdirectories,
        filenames,
        labels,
        black_line_flag,
        dashed_line_flag,
        figure_directory_base,
        data_directory_base,
        plot_reference_result,
        figure_filename_postfix,
        figure_title,
        log_axes_input,
        legend_on_input,
        legend_inside_input,
        nlegendcols_input,
        # clr_input=clr_input,
        transparent_legend_input=True,
        tmax=12.5,#14
        legend_fontSize_input=14,
        solid_and_dashed_lines=False,
        plot_numerical_dissipation=True)
    #-----------------------------------------------------
    # clr_input = ['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs" # comment to turn off
    figure_filename_postfix = "96_p5_cDG_cfl_physically_consistent"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories.append("time_step_advantage_with_physical_check/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.26_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GL CFL=$0.26$")
    black_line_flag.append(False)
    dashed_line_flag.append(True)
    #-----------------------------------------------------
    subdirectories.append("time_step_advantage_with_physical_check/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.24_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GL CFL=$0.24$")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("time_step_advantage_with_physical_check/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.22_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GL CFL=$0.22$")
    black_line_flag.append(False)
    dashed_line_flag.append(True)
    #-----------------------------------------------------
    subdirectories.append("time_step_advantage/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.2_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GL CFL=$0.20$")
    black_line_flag.append(False)
    dashed_line_flag.append(True)
    #-----------------------------------------------------
    subdirectories.append("flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GL CFL=$0.10$")
    black_line_flag.append(False)
    dashed_line_flag.append(True)
    #-----------------------------------------------------
    plot_periodic_turbulence(
        figure_subdirectory,
        subdirectories,
        filenames,
        labels,
        black_line_flag,
        dashed_line_flag,
        figure_directory_base,
        data_directory_base,
        plot_reference_result,
        figure_filename_postfix,
        figure_title,
        log_axes_input,
        legend_on_input,
        legend_inside_input,
        nlegendcols_input,
        # clr_input=clr_input,
        transparent_legend_input=True,
        tmax=12.5,#14
        legend_fontSize_input=14,
        solid_and_dashed_lines=False,
        plot_numerical_dissipation=True)
    # #-----------------------------------------------------
    # reinit_inputs()
    # data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    # date_for_runs="."
    # figure_subdirectory="2023_JCP"
    # figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL$=0.10$" # comment to turn off
    # figure_filename_postfix = "96_p5_strong_DG_cfl"
    # legend_inside_input=True
    # #-----------------------------------------------------
    # subdirectories.append("flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512")
    # filenames.append("turbulent_quantities.txt")
    # # labels.append("Strong DG, Roe, GL flux nodes,\n $n_{quad}=2(P+1)$")
    # labels.append("Strong DG-Roe-GL-OI CFL=$0.10$")
    # black_line_flag.append(False)
    # dashed_line_flag.append(False)
    # #-----------------------------------------------------
    # subdirectories.append("time_step_advantage_strong_DG/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_CFL-0.13_procs512")
    # filenames.append("turbulent_quantities.txt")
    # # labels.append("Strong DG, Roe, GL flux nodes,\n $n_{quad}=2(P+1)$")
    # labels.append("Strong DG-Roe-GL-OI CFL=$0.13$")
    # black_line_flag.append(False)
    # dashed_line_flag.append(True)
    # #-----------------------------------------------------
    # subdirectories.append("time_step_advantage_strong_DG/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_CFL-0.14_procs512")
    # filenames.append("turbulent_quantities.txt")
    # # labels.append("Strong DG, Roe, GL flux nodes,\n $n_{quad}=2(P+1)$")
    # labels.append("Strong DG-Roe-GL-OI CFL=$0.14$")
    # black_line_flag.append(False)
    # dashed_line_flag.append(True)
    # #-----------------------------------------------------
    # subdirectories.append("time_step_advantage_strong_DG/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_CFL-0.15_procs512")
    # filenames.append("turbulent_quantities.txt")
    # # labels.append("Strong DG, Roe, GL flux nodes,\n $n_{quad}=2(P+1)$")
    # labels.append("Strong DG-Roe-GL-OI CFL=$0.15$")
    # black_line_flag.append(False)
    # dashed_line_flag.append(True)
    # #-----------------------------------------------------
    # plot_periodic_turbulence(
    #     figure_subdirectory,
    #     subdirectories,
    #     filenames,
    #     labels,
    #     black_line_flag,
    #     dashed_line_flag,
    #     figure_directory_base,
    #     data_directory_base,
    #     plot_reference_result,
    #     figure_filename_postfix,
    #     figure_title,
    #     log_axes_input,
    #     legend_on_input,
    #     legend_inside_input,
    #     nlegendcols_input,
    #     # clr_input=clr_input,
    #     transparent_legend_input=True,
    #     tmax=12.5,#14
    #     legend_fontSize_input=14,
    #     solid_and_dashed_lines=False,
    #     plot_numerical_dissipation=True)
    # #-----------------------------------------------------
    # # clr_input = ['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    # reinit_inputs()
    # data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    # date_for_runs="."
    # figure_subdirectory="2023_JCP"
    # figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL$=0.10$" # comment to turn off
    # figure_filename_postfix = "96_p5_cfl_advantage"
    # legend_inside_input=True
    # #-----------------------------------------------------
    # subdirectories.append("flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512")
    # filenames.append("turbulent_quantities.txt")
    # # labels.append("Strong DG, Roe, GL flux nodes,\n $n_{quad}=2(P+1)$")
    # labels.append("Strong DG-Roe-GL-OI CFL=$0.10$")
    # black_line_flag.append(False)
    # dashed_line_flag.append(False)
    # #-----------------------------------------------------
    # subdirectories.append("flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512")
    # filenames.append("turbulent_quantities.txt")
    # labels.append("$c_{DG}$ NSFR.IR-GL CFL=$0.10$")
    # black_line_flag.append(False)
    # dashed_line_flag.append(False)
    # #-----------------------------------------------------
    # subdirectories.append("correction_parameter/viscous_TGV_ILES_NSFR_cPlus_IR_2PF_GL_OI-0_dofs096_p5_procs512")
    # filenames.append("turbulent_quantities.txt")
    # labels.append("$c_{+}$ NSFR.IR-GL CFL=$0.10$")
    # black_line_flag.append(False)
    # dashed_line_flag.append(False)
    # #-----------------------------------------------------
    # subdirectories.append("time_step_advantage_strong_DG/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_CFL-0.14_procs512")
    # filenames.append("turbulent_quantities.txt")
    # # labels.append("Strong DG, Roe, GL flux nodes,\n $n_{quad}=2(P+1)$")
    # labels.append("Strong DG-Roe-GL-OI CFL=$0.14$")
    # black_line_flag.append(False)
    # dashed_line_flag.append(True)
    # #-----------------------------------------------------
    # subdirectories.append("time_step_advantage/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.34_procs512")
    # filenames.append("turbulent_quantities.txt")
    # labels.append("$c_{DG}$ NSFR.IR-GL CFL=$0.34$")
    # black_line_flag.append(False)
    # dashed_line_flag.append(True)
    # #-----------------------------------------------------
    # subdirectories.append("time_step_advantage/viscous_TGV_ILES_NSFR_cPlus_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.44_procs512")
    # filenames.append("turbulent_quantities.txt")
    # labels.append("$c_{+}$ NSFR.IR-GL CFL=$0.44$")
    # black_line_flag.append(False)
    # dashed_line_flag.append(True)
    # #-----------------------------------------------------
    # plot_periodic_turbulence(
    #     figure_subdirectory,
    #     subdirectories,
    #     filenames,
    #     labels,
    #     black_line_flag,
    #     dashed_line_flag,
    #     figure_directory_base,
    #     data_directory_base,
    #     plot_reference_result,
    #     figure_filename_postfix,
    #     figure_title,
    #     log_axes_input,
    #     legend_on_input,
    #     legend_inside_input,
    #     nlegendcols_input,
    #     # clr_input=clr_input,
    #     transparent_legend_input=True,
    #     tmax=12.5,#14
    #     legend_fontSize_input=14,
    #     solid_and_dashed_lines=False,
    #     plot_numerical_dissipation=True)
    # #-----------------------------------------------------

#=====================================================
# DOFs: 96^3 | Time step advantage (CFL advantage)
#-----------------------------------------------------
if(False):
    #-----------------------------------------------------
    # clr_input = ['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    reinit_inputs()
    # /Volumes/Samsung_T5/NarvalFiles
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL$=0.10$" # comment to turn off
    figure_filename_postfix = "96_p5_cPlus_cfl"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories.append("correction_parameter/viscous_TGV_ILES_NSFR_cPlus_IR_2PF_GL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{+}$ NSFR.IR-GL CFL=$0.10$")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("time_step_advantage/viscous_TGV_ILES_NSFR_cPlus_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.2_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{+}$ NSFR.IR-GL CFL=$0.20$")
    black_line_flag.append(False)
    dashed_line_flag.append(True)
    #-----------------------------------------------------
    subdirectories.append("time_step_advantage/viscous_TGV_ILES_NSFR_cPlus_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.3_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{+}$ NSFR.IR-GL CFL=$0.30$")
    black_line_flag.append(False)
    dashed_line_flag.append(True)
    #-----------------------------------------------------
    subdirectories.append("time_step_advantage/viscous_TGV_ILES_NSFR_cPlus_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.42_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{+}$ NSFR.IR-GL CFL=$0.42$")
    black_line_flag.append(False)
    dashed_line_flag.append(True)
    #-----------------------------------------------------
    subdirectories.append("time_step_advantage/viscous_TGV_ILES_NSFR_cPlus_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.44_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{+}$ NSFR.IR-GL CFL=$0.44$")
    black_line_flag.append(False)
    dashed_line_flag.append(True)
    #-----------------------------------------------------
    plot_periodic_turbulence(
        figure_subdirectory,
        subdirectories,
        filenames,
        labels,
        black_line_flag,
        dashed_line_flag,
        figure_directory_base,
        data_directory_base,
        plot_reference_result,
        figure_filename_postfix,
        figure_title,
        log_axes_input,
        legend_on_input,
        legend_inside_input,
        nlegendcols_input,
        # clr_input=clr_input,
        transparent_legend_input=True,
        tmax=12.5,#14
        legend_fontSize_input=14,
        solid_and_dashed_lines=False,
        plot_numerical_dissipation=True)
    #-----------------------------------------------------
    # clr_input = ['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL$=0.10$" # comment to turn off
    figure_filename_postfix = "96_p5_cDG_cfl"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories.append("flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GL CFL=$0.10$")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("time_step_advantage/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.2_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GL CFL=$0.20$")
    black_line_flag.append(False)
    dashed_line_flag.append(True)
    #-----------------------------------------------------
    subdirectories.append("time_step_advantage/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.32_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GL CFL=$0.32$")
    black_line_flag.append(False)
    dashed_line_flag.append(True)
    #-----------------------------------------------------
    subdirectories.append("time_step_advantage/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.34_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GL CFL=$0.34$")
    black_line_flag.append(False)
    dashed_line_flag.append(True)
    #-----------------------------------------------------
    plot_periodic_turbulence(
        figure_subdirectory,
        subdirectories,
        filenames,
        labels,
        black_line_flag,
        dashed_line_flag,
        figure_directory_base,
        data_directory_base,
        plot_reference_result,
        figure_filename_postfix,
        figure_title,
        log_axes_input,
        legend_on_input,
        legend_inside_input,
        nlegendcols_input,
        # clr_input=clr_input,
        transparent_legend_input=True,
        tmax=12.5,#14
        legend_fontSize_input=14,
        solid_and_dashed_lines=False,
        plot_numerical_dissipation=True)
    #-----------------------------------------------------
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL$=0.10$" # comment to turn off
    figure_filename_postfix = "96_p5_strong_DG_cfl"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories.append("flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    # labels.append("Strong DG, Roe, GL flux nodes,\n $n_{quad}=2(P+1)$")
    labels.append("Strong DG-Roe-GL-OI CFL=$0.10$")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    # #-----------------------------------------------------
    # subdirectories.append("time_step_advantage_strong_DG/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_CFL-0.12_procs512")
    # filenames.append("turbulent_quantities.txt")
    # # labels.append("Strong DG, Roe, GL flux nodes,\n $n_{quad}=2(P+1)$")
    # labels.append("Strong DG-Roe-GL-OI CFL=$0.12$")
    # black_line_flag.append(False)
    # dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("time_step_advantage_strong_DG/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_CFL-0.13_procs512")
    filenames.append("turbulent_quantities.txt")
    # labels.append("Strong DG, Roe, GL flux nodes,\n $n_{quad}=2(P+1)$")
    labels.append("Strong DG-Roe-GL-OI CFL=$0.13$")
    black_line_flag.append(False)
    dashed_line_flag.append(True)
    #-----------------------------------------------------
    subdirectories.append("time_step_advantage_strong_DG/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_CFL-0.14_procs512")
    filenames.append("turbulent_quantities.txt")
    # labels.append("Strong DG, Roe, GL flux nodes,\n $n_{quad}=2(P+1)$")
    labels.append("Strong DG-Roe-GL-OI CFL=$0.14$")
    black_line_flag.append(False)
    dashed_line_flag.append(True)
    #-----------------------------------------------------
    subdirectories.append("time_step_advantage_strong_DG/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_CFL-0.15_procs512")
    filenames.append("turbulent_quantities.txt")
    # labels.append("Strong DG, Roe, GL flux nodes,\n $n_{quad}=2(P+1)$")
    labels.append("Strong DG-Roe-GL-OI CFL=$0.15$")
    black_line_flag.append(False)
    dashed_line_flag.append(True)
    #-----------------------------------------------------
    plot_periodic_turbulence(
        figure_subdirectory,
        subdirectories,
        filenames,
        labels,
        black_line_flag,
        dashed_line_flag,
        figure_directory_base,
        data_directory_base,
        plot_reference_result,
        figure_filename_postfix,
        figure_title,
        log_axes_input,
        legend_on_input,
        legend_inside_input,
        nlegendcols_input,
        # clr_input=clr_input,
        transparent_legend_input=True,
        tmax=12.5,#14
        legend_fontSize_input=14,
        solid_and_dashed_lines=False,
        plot_numerical_dissipation=True)
    #-----------------------------------------------------
    # clr_input = ['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL$=0.10$" # comment to turn off
    figure_filename_postfix = "96_p5_cfl_advantage"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories.append("flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    # labels.append("Strong DG, Roe, GL flux nodes,\n $n_{quad}=2(P+1)$")
    labels.append("Strong DG-Roe-GL-OI CFL=$0.10$")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GL CFL=$0.10$")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("correction_parameter/viscous_TGV_ILES_NSFR_cPlus_IR_2PF_GL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{+}$ NSFR.IR-GL CFL=$0.10$")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("time_step_advantage_strong_DG/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_CFL-0.14_procs512")
    filenames.append("turbulent_quantities.txt")
    # labels.append("Strong DG, Roe, GL flux nodes,\n $n_{quad}=2(P+1)$")
    labels.append("Strong DG-Roe-GL-OI CFL=$0.14$")
    black_line_flag.append(False)
    dashed_line_flag.append(True)
    #-----------------------------------------------------
    subdirectories.append("time_step_advantage/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.34_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GL CFL=$0.34$")
    black_line_flag.append(False)
    dashed_line_flag.append(True)
    #-----------------------------------------------------
    subdirectories.append("time_step_advantage/viscous_TGV_ILES_NSFR_cPlus_IR_2PF_GL_OI-0_dofs096_p5_CFL-0.44_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{+}$ NSFR.IR-GL CFL=$0.44$")
    black_line_flag.append(False)
    dashed_line_flag.append(True)
    #-----------------------------------------------------
    plot_periodic_turbulence(
        figure_subdirectory,
        subdirectories,
        filenames,
        labels,
        black_line_flag,
        dashed_line_flag,
        figure_directory_base,
        data_directory_base,
        plot_reference_result,
        figure_filename_postfix,
        figure_title,
        log_axes_input,
        legend_on_input,
        legend_inside_input,
        nlegendcols_input,
        # clr_input=clr_input,
        transparent_legend_input=True,
        tmax=12.5,#14
        legend_fontSize_input=14,
        solid_and_dashed_lines=False,
        plot_numerical_dissipation=True)
    #-----------------------------------------------------

#=====================================================
# DOFs: 96^3 | NSFR.IR-cDG-GL and NSFR.IR-cDG-GLL (BASELINE)
#-----------------------------------------------------
if(False):
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL$=0.10$" # comment to turn off
    figure_filename_postfix = "96_cDG_roe_vs_sDG"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories.append("upwind_dissipation_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF-Roe_GL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GL-Roe")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    # labels.append("Strong DG, Roe, GL flux nodes,\n $n_{quad}=2(P+1)$")
    labels.append("Strong DG-Roe-GL-OI")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    plot_periodic_turbulence(
        figure_subdirectory,
        subdirectories,
        filenames,
        labels,
        black_line_flag,
        dashed_line_flag,
        figure_directory_base,
        data_directory_base,
        plot_reference_result,
        figure_filename_postfix,
        figure_title,
        log_axes_input,
        legend_on_input,
        legend_inside_input,
        nlegendcols_input,
        transparent_legend_input=True,
        tmax=20.0,
        legend_fontSize_input=12,
        solid_and_dashed_lines=False,
        plot_numerical_dissipation=True)
    #-----------------------------------------------------

#=====================================================
# DOFs: 96^3 | NSFR.IR-cDG-GL and NSFR.IR-cDG-GLL (BASELINE)
#-----------------------------------------------------
if(True):
    #-----------------------------------------------------
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL$=0.10$" # comment to turn off
    figure_filename_postfix = "96_baseline_scheme"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512",\
    "flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs096_p5_procs512",\
    ]
    # labels
    labels_for_plot=[\
    "$c_{DG}$ NSFR.IR-GL",\
    "$c_{DG}$ NSFR.IR-GLL",\
    ]
    black_line_flag_for_plot=[\
    False,\
    False,\
    ]
    dashed_line_flag_for_plot=[\
    False,\
    False,\
    ]
    smoothing_input = [\
    False,\
    False,\
    ]
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot)

#=====================================================
# DOFs: 256^3 | NSFR.IR-cDG-GL (BASELINE)
#-----------------------------------------------------
if(False):
    clr_input = ['tab:blue','tab:red','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    # figure_title = "TGV at Re$_{\\infty}=1600$, P$3$, $256^{3}$ DOFs, CFL$=0.30$" # comment to turn off
    figure_title = "TGV at Re$_{\\infty}=1600$, P$3$, $256^{3}$ DOFs" # comment to turn off
    figure_filename_postfix = "256_verification"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories.append("verification/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p3_procs1024")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GL, CFL$=0.30$")
    # labels.append("$c_{DG}$ NSFR.IR-GL")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("verification/viscous_TGV_ILES_NSFR_cDG_IR_2PF-Roe_GL_OI-0_dofs0256_p3_procs1024")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GL-Roe, CFL$=0.15$")
    black_line_flag.append(False)
    dashed_line_flag.append(True)
    #-----------------------------------------------------
    plot_periodic_turbulence(
        figure_subdirectory,
        subdirectories,
        filenames,
        labels,
        black_line_flag,
        dashed_line_flag,
        figure_directory_base,
        data_directory_base,
        plot_reference_result,
        figure_filename_postfix,
        figure_title,
        log_axes_input,
        legend_on_input,
        legend_inside_input,
        nlegendcols_input,
        clr_input=clr_input,
        transparent_legend_input=True,
        tmax=20.0,
        legend_fontSize_input=12,
        solid_and_dashed_lines=False,
        plot_numerical_dissipation=True)
    #-----------------------------------------------------
if(True):
    clr_input = ['tab:blue','tab:red','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    # figure_title = "TGV at Re$_{\\infty}=1600$, P$3$, $256^{3}$ DOFs, CFL$=0.30$" # comment to turn off
    figure_title = "TGV at Re$_{\\infty}=1600$, P$3$, $256^{3}$ DOFs" # comment to turn off
    figure_filename_postfix = "256_p3_OI_vs_SF"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories_for_plot=[\
    "verification/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-4_dofs0256_p3_CFL-0.15_procs1024",\
    "verification/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p3_procs1024",\
    ]
    # labels
    labels_for_plot=[\
    "Strong DG-Roe-GL-OI",\
    "$c_{DG}$ NSFR.IR-GL",\
    ]
    black_line_flag_for_plot=[\
    False,\
    False,\
    ]
    dashed_line_flag_for_plot=[\
    False,\
    False,\
    ]
    smoothing_input = [\
    False,\
    False,\
    ]
    plot_for_presentation(subdirectories_for_plot,labels_for_plot,black_line_flag_for_plot,dashed_line_flag_for_plot,final_time_for_plot=10.0)

#=====================================================
# DOFs: 96^3 | High order poly instabilities
#-----------------------------------------------------
if(False):
    # clr_input = ['tab:red','tab:blue','tab:green','tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    figure_title = "TGV at Re$_{\\infty}=1600$, P$7$, $64^{3}$ DOFs, CFL$=0.10$" # comment to turn off
    figure_filename_postfix = "64_high_poly_degree"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories.append("high_poly_degree_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs064_p7_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GL")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("high_poly_degree_GL_flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-0_dofs064_p7_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("Strong DG-Roe-GL")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("high_poly_degree_GL_flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-8_dofs064_p7_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("Strong DG-Roe-GL-OI")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    plot_periodic_turbulence(
        figure_subdirectory,
        subdirectories,
        filenames,
        labels,
        black_line_flag,
        dashed_line_flag,
        figure_directory_base,
        data_directory_base,
        plot_reference_result,
        figure_filename_postfix,
        figure_title,
        log_axes_input,
        legend_on_input,
        legend_inside_input,
        nlegendcols_input,
        # clr_input=clr_input,
        transparent_legend_input=True,
        tmax=14.0,
        legend_fontSize_input=14,
        solid_and_dashed_lines=False)
    #-----------------------------------------------------

#=====================================================
# DOFs: 96^3 | NSFR.IR-cDG (GL and GLL) WITH OVERINTEGRATION
#-----------------------------------------------------
if(False):
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL$=0.10$" # comment to turn off
    figure_filename_postfix = "96_overintegration"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories.append("flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GL")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GLL")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("over_integration/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-6_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GL-OI")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("over_integration/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-6_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GLL-OI")
    black_line_flag.append(False)
    dashed_line_flag.append(True)
    #-----------------------------------------------------
    plot_periodic_turbulence(
        figure_subdirectory,
        subdirectories,
        filenames,
        labels,
        black_line_flag,
        dashed_line_flag,
        figure_directory_base,
        data_directory_base,
        plot_reference_result,
        figure_filename_postfix,
        figure_title,
        log_axes_input,
        legend_on_input,
        legend_inside_input,
        nlegendcols_input,
        transparent_legend_input=True,
        tmax=20.0,
        legend_fontSize_input=14,
        solid_and_dashed_lines=False)
    #-----------------------------------------------------

#=====================================================
# DOFs: 96^3 | Flux Nodes (cDG NSFR vs Strong DG)
#-----------------------------------------------------
if(False):
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL$=0.10$" # comment to turn off
    figure_filename_postfix = "96_flux_nodes"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories.append("flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GL")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GLL")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    # labels.append("Strong DG, Roe, GL flux nodes,\n $n_{quad}=2(P+1)$")
    labels.append("Strong DG-Roe-GL-OI")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GLL_OI-6_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    # labels.append("Strong DG, Roe, GLL flux nodes,\n $n_{quad}=2(P+1)$")
    labels.append("Strong DG-Roe-GLL-OI")
    black_line_flag.append(False)
    dashed_line_flag.append(True)
    #-----------------------------------------------------
    plot_periodic_turbulence(
        figure_subdirectory,
        subdirectories,
        filenames,
        labels,
        black_line_flag,
        dashed_line_flag,
        figure_directory_base,
        data_directory_base,
        plot_reference_result,
        figure_filename_postfix,
        figure_title,
        log_axes_input,
        legend_on_input,
        legend_inside_input,
        nlegendcols_input,
        transparent_legend_input=True,
        tmax=20.0,
        legend_fontSize_input=14,
        solid_and_dashed_lines=False)
    #-----------------------------------------------------

#=====================================================
# DOFs: 96^3 | Upwind dissipation on GLL flux nodes
#-----------------------------------------------------
if(False):
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL$=0.10$" # comment to turn off
    figure_filename_postfix = "96_upwind_gll"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories.append("flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GLL")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("upwind_dissipation_GLL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF-LxF_GLL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GLL-LxF")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("upwind_dissipation_GLL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF-Roe_GLL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GLL-Roe")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("upwind_dissipation_GLL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF-L2R_GLL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GLL-L2R")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    plot_periodic_turbulence(
        figure_subdirectory,
        subdirectories,
        filenames,
        labels,
        black_line_flag,
        dashed_line_flag,
        figure_directory_base,
        data_directory_base,
        plot_reference_result,
        figure_filename_postfix,
        figure_title,
        log_axes_input,
        legend_on_input,
        legend_inside_input,
        nlegendcols_input,
        transparent_legend_input=True,
        tmax=20.0,
        legend_fontSize_input=14,
        solid_and_dashed_lines=False)
    #-----------------------------------------------------

#=====================================================
# DOFs: 96^3 | Filter width stabilization
#-----------------------------------------------------
if(False):
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL$=0.10$" # comment to turn off
    figure_filename_postfix = "96_filter_width_stabilization"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories.append("flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    # labels.append("Strong DG, Roe, GL flux nodes,\n $n_{quad}=2(P+1)$")
    labels.append("Strong DG-Roe-GL-OI")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("filter_width_stabilization/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("Strong DG-Roe-GL")
    black_line_flag.append(False)
    dashed_line_flag.append(True)
    #-----------------------------------------------------
    subdirectories.append("filter_width_stabilization/viscous_TGV_LES_SMAG_MC-0.18_std_strong_DG_Roe_GL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("Strong DG-Roe-GL-Smag. $C_{S}=0.18$")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("filter_width_stabilization/viscous_TGV_LES_SMAG_MC-0.18_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("Strong DG-Roe-GL-OI-Smag. $C_{S}=0.18$")
    black_line_flag.append(False)
    dashed_line_flag.append(True)
    #-----------------------------------------------------
    # subdirectories.append("flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GLL_OI-6_dofs096_p5_procs512")
    # filenames.append("turbulent_quantities.txt")
    # # labels.append("Strong DG, Roe, GL flux nodes,\n $n_{quad}=2(P+1)$")
    # labels.append("Strong DG-Roe-GLL-OI")
    # black_line_flag.append(False)
    # dashed_line_flag.append(False)
    # #-----------------------------------------------------
    # subdirectories.append("filter_width_stabilization/viscous_TGV_ILES_std_strong_DG_Roe_GLL_OI-0_dofs096_p5_procs512")
    # filenames.append("turbulent_quantities.txt")
    # labels.append("Strong DG-Roe-GLL")
    # black_line_flag.append(False)
    # dashed_line_flag.append(True)
    # #-----------------------------------------------------
    # subdirectories.append("filter_width_stabilization/viscous_TGV_LES_SMAG_MC-0.18_std_strong_DG_Roe_GLL_OI-0_dofs096_p5_procs512")
    # filenames.append("turbulent_quantities.txt")
    # labels.append("Strong DG-Roe-GLL-Smag. $C_{S}=0.18$")
    # black_line_flag.append(False)
    # dashed_line_flag.append(False)
    # #-----------------------------------------------------
    # subdirectories.append("filter_width_stabilization/viscous_TGV_LES_SMAG_MC-0.18_std_strong_DG_Roe_GLL_OI-6_dofs096_p5_procs512")
    # filenames.append("turbulent_quantities.txt")
    # labels.append("Strong DG-Roe-GLL-OI-Smag. $C_{S}=0.18$")
    # black_line_flag.append(False)
    # dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("filter_width_stabilization/viscous_TGV_LES_SMAG_MC-0.18_std_strong_DG_Roe_GL_OI-0_dofs096_p5_procs512_filter_width_flad_and_gassner")
    filenames.append("turbulent_quantities.txt")
    labels.append("Strong DG-Roe-GL-Smag. $C_{S}=0.18$ $\\Delta=\\frac{V}{(P+1)^{3}}$")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("filter_width_stabilization/viscous_TGV_LES_SMAG_MC-0.18_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512_filter_width_flad_and_gassner")
    filenames.append("turbulent_quantities.txt")
    labels.append("Strong DG-Roe-GL-OI-Smag. $C_{S}=0.18$ $\\Delta=\\frac{V}{(P+1)^{3}}$")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    plot_periodic_turbulence(
        figure_subdirectory,
        subdirectories,
        filenames,
        labels,
        black_line_flag,
        dashed_line_flag,
        figure_directory_base,
        data_directory_base,
        plot_reference_result,
        figure_filename_postfix,
        figure_title,
        log_axes_input,
        legend_on_input,
        legend_inside_input,
        nlegendcols_input,
        transparent_legend_input=True,
        tmax=20.0,
        legend_fontSize_input=14,
        solid_and_dashed_lines=False)
    #-----------------------------------------------------

#=====================================================
# DOFs: 96^3 | Filter width stabilization
#-----------------------------------------------------
if(False):
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL$=0.10$" # comment to turn off
    figure_filename_postfix = "96_filter_width_stabilization_gll"
    legend_inside_input=True
    # #-----------------------------------------------------
    # subdirectories.append("flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512")
    # filenames.append("turbulent_quantities.txt")
    # # labels.append("Strong DG, Roe, GL flux nodes,\n $n_{quad}=2(P+1)$")
    # labels.append("Strong DG-Roe-GL-OI")
    # black_line_flag.append(False)
    # dashed_line_flag.append(False)
    # #-----------------------------------------------------
    # subdirectories.append("filter_width_stabilization/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-0_dofs096_p5_procs512")
    # filenames.append("turbulent_quantities.txt")
    # labels.append("Strong DG-Roe-GL")
    # black_line_flag.append(False)
    # dashed_line_flag.append(True)
    # #-----------------------------------------------------
    # subdirectories.append("filter_width_stabilization/viscous_TGV_LES_SMAG_MC-0.18_std_strong_DG_Roe_GL_OI-0_dofs096_p5_procs512")
    # filenames.append("turbulent_quantities.txt")
    # labels.append("Strong DG-Roe-GL-Smag. $C_{S}=0.18$")
    # black_line_flag.append(False)
    # dashed_line_flag.append(False)
    # #-----------------------------------------------------
    # subdirectories.append("filter_width_stabilization/viscous_TGV_LES_SMAG_MC-0.18_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512")
    # filenames.append("turbulent_quantities.txt")
    # labels.append("Strong DG-Roe-GL-OI-Smag. $C_{S}=0.18$")
    # black_line_flag.append(False)
    # dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GLL_OI-6_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    # labels.append("Strong DG, Roe, GL flux nodes,\n $n_{quad}=2(P+1)$")
    labels.append("Strong DG-Roe-GLL-OI")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("filter_width_stabilization/viscous_TGV_ILES_std_strong_DG_Roe_GLL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("Strong DG-Roe-GLL")
    black_line_flag.append(False)
    dashed_line_flag.append(True)
    #-----------------------------------------------------
    subdirectories.append("filter_width_stabilization/viscous_TGV_LES_SMAG_MC-0.18_std_strong_DG_Roe_GLL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("Strong DG-Roe-GLL-Smag. $C_{S}=0.18$")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    # #-----------------------------------------------------
    # subdirectories.append("filter_width_stabilization/viscous_TGV_LES_SMAG_MC-0.18_std_strong_DG_Roe_GLL_OI-6_dofs096_p5_procs512")
    # filenames.append("turbulent_quantities.txt")
    # labels.append("Strong DG-Roe-GLL-OI-Smag. $C_{S}=0.18$")
    # black_line_flag.append(False)
    # dashed_line_flag.append(False)
    # #-----------------------------------------------------
    # subdirectories.append("filter_width_stabilization/viscous_TGV_LES_SMAG_MC-0.18_std_strong_DG_Roe_GL_OI-0_dofs096_p5_procs512_filter_width_flad_and_gassner")
    # filenames.append("turbulent_quantities.txt")
    # labels.append("Strong DG-Roe-GL-Smag. $C_{S}=0.18$ $\\Delta=\\frac{V}{(P+1)^{3}}$")
    # black_line_flag.append(False)
    # dashed_line_flag.append(False)
    # #-----------------------------------------------------
    # subdirectories.append("filter_width_stabilization/viscous_TGV_LES_SMAG_MC-0.18_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512_filter_width_flad_and_gassner")
    # filenames.append("turbulent_quantities.txt")
    # labels.append("Strong DG-Roe-GL-OI-Smag. $C_{S}=0.18$ $\\Delta=\\frac{V}{(P+1)^{3}}$")
    # black_line_flag.append(False)
    # dashed_line_flag.append(False)
    #-----------------------------------------------------
    plot_periodic_turbulence(
        figure_subdirectory,
        subdirectories,
        filenames,
        labels,
        black_line_flag,
        dashed_line_flag,
        figure_directory_base,
        data_directory_base,
        plot_reference_result,
        figure_filename_postfix,
        figure_title,
        log_axes_input,
        legend_on_input,
        legend_inside_input,
        nlegendcols_input,
        transparent_legend_input=True,
        tmax=20.0,
        legend_fontSize_input=14,
        solid_and_dashed_lines=False)
    #-----------------------------------------------------

# #=====================================================
# # DOFs: 96^3 | Filter width stabilization
# #-----------------------------------------------------
# if(False):
#     reinit_inputs()
#     data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
#     date_for_runs="."
#     figure_subdirectory="2023_JCP"
#     figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL$=0.10$" # comment to turn off
#     figure_filename_postfix = "96_filter_width_stabilization"
#     legend_inside_input=True
#     #-----------------------------------------------------
#     subdirectories.append("flux_nodes/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512")
#     filenames.append("turbulent_quantities.txt")
#     # labels.append("Strong DG, Roe, GL flux nodes,\n $n_{quad}=2(P+1)$")
#     labels.append("Strong DG-Roe-GL-OI")
#     black_line_flag.append(False)
#     dashed_line_flag.append(False)
#     #-----------------------------------------------------
#     subdirectories.append("filter_width_stabilization/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-0_dofs096_p5_procs512")
#     filenames.append("turbulent_quantities.txt")
#     labels.append("Strong DG-Roe-GL")
#     black_line_flag.append(False)
#     dashed_line_flag.append(True)
#     #-----------------------------------------------------
#     subdirectories.append("filter_width_stabilization/viscous_TGV_LES_SMAG_MC-0.18_std_strong_DG_Roe_GL_OI-0_dofs096_p5_procs512")
#     filenames.append("turbulent_quantities.txt")
#     labels.append("Strong DG-Roe-GL-Smag. $C_{S}=0.18$")
#     black_line_flag.append(False)
#     dashed_line_flag.append(False)
#     #-----------------------------------------------------
#     subdirectories.append("filter_width_stabilization/viscous_TGV_ILES_std_strong_DG_Roe_GLL_OI-0_dofs096_p5_procs512")
#     filenames.append("turbulent_quantities.txt")
#     labels.append("Strong DG-Roe-GLL")
#     black_line_flag.append(False)
#     dashed_line_flag.append(True)
#     #-----------------------------------------------------
#     subdirectories.append("filter_width_stabilization/viscous_TGV_LES_SMAG_MC-0.18_std_strong_DG_Roe_GLL_OI-0_dofs096_p5_procs512")
#     filenames.append("turbulent_quantities.txt")
#     labels.append("Strong DG-Roe-GLL-Smag. $C_{S}=0.18$")
#     black_line_flag.append(False)
#     dashed_line_flag.append(False)
#     #-----------------------------------------------------
#     subdirectories.append("filter_width_stabilization/viscous_TGV_LES_SMAG_MC-0.18_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512")
#     filenames.append("turbulent_quantities.txt")
#     labels.append("Strong DG-Roe-GL-OI-Smag. $C_{S}=0.18$")
#     black_line_flag.append(False)
#     dashed_line_flag.append(False)
#     #-----------------------------------------------------
#     subdirectories.append("filter_width_stabilization/viscous_TGV_LES_SMAG_MC-0.18_std_strong_DG_Roe_GLL_OI-6_dofs096_p5_procs512")
#     filenames.append("turbulent_quantities.txt")
#     labels.append("Strong DG-Roe-GLL-OI-Smag. $C_{S}=0.18$")
#     black_line_flag.append(False)
#     dashed_line_flag.append(False)
#     #-----------------------------------------------------
#     subdirectories.append("filter_width_stabilization/viscous_TGV_LES_SMAG_MC-0.18_std_strong_DG_Roe_GL_OI-0_dofs096_p5_procs512_filter_width_flad_and_gassner")
#     filenames.append("turbulent_quantities.txt")
#     labels.append("Strong DG-Roe-GL-Smag. $C_{S}=0.18$ $\\Delta=\\frac{V}{(P+1)^{3}}$")
#     black_line_flag.append(False)
#     dashed_line_flag.append(False)
#     #-----------------------------------------------------
#     subdirectories.append("filter_width_stabilization/viscous_TGV_LES_SMAG_MC-0.18_std_strong_DG_Roe_GL_OI-6_dofs096_p5_procs512_filter_width_flad_and_gassner")
#     filenames.append("turbulent_quantities.txt")
#     labels.append("Strong DG-Roe-GL-OI-Smag. $C_{S}=0.18$ $\\Delta=\\frac{V}{(P+1)^{3}}$")
#     black_line_flag.append(False)
#     dashed_line_flag.append(False)
#     #-----------------------------------------------------
#     plot_periodic_turbulence(
#         figure_subdirectory,
#         subdirectories,
#         filenames,
#         labels,
#         black_line_flag,
#         dashed_line_flag,
#         figure_directory_base,
#         data_directory_base,
#         plot_reference_result,
#         figure_filename_postfix,
#         figure_title,
#         log_axes_input,
#         legend_on_input,
#         legend_inside_input,
#         nlegendcols_input,
#         transparent_legend_input=True,
#         tmax=20.0,
#         legend_fontSize_input=14,
#         solid_and_dashed_lines=False)
#     #-----------------------------------------------------

#=====================================================
# DOFs: 96^3 | NSFR Correction Parameter
#-----------------------------------------------------
if(False):
    reinit_inputs()
    data_directory_base=filesystem+"NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL$=0.10$" # comment to turn off
    figure_filename_postfix = "96_correction_parameter"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories.append("flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GL")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("correction_parameter/viscous_TGV_ILES_NSFR_cPlus_IR_2PF_GL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{+}$ NSFR.IR-GL")
    black_line_flag.append(False)
    dashed_line_flag.append(True)
    #-----------------------------------------------------
    subdirectories.append("correction_parameter/viscous_TGV_ILES_NSFR_cHU_IR_2PF_GL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{HU}$ NSFR.IR-GL")
    black_line_flag.append(False)
    dashed_line_flag.append(True)
    #-----------------------------------------------------
    subdirectories.append("correction_parameter/viscous_TGV_ILES_NSFR_cSD_IR_2PF_GL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{SD}$ NSFR.IR-GL")
    black_line_flag.append(False)
    dashed_line_flag.append(True)
    #-----------------------------------------------------
    plot_periodic_turbulence(
        figure_subdirectory,
        subdirectories,
        filenames,
        labels,
        black_line_flag,
        dashed_line_flag,
        figure_directory_base,
        data_directory_base,
        plot_reference_result,
        figure_filename_postfix,
        figure_title,
        log_axes_input,
        legend_on_input,
        legend_inside_input,
        nlegendcols_input,
        transparent_legend_input=True,
        tmax=20.0,
        legend_fontSize_input=14,
        solid_and_dashed_lines=False)
    #-----------------------------------------------------



