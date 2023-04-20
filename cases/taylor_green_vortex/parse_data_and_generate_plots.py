import os;CURRENT_PATH = os.path.split(os.path.realpath(__file__))[0]+"/";
import sys
sys.path.append(CURRENT_PATH+"../../src");
from plot_unsteady_integrated_turbulent_flow_quantities import plot_periodic_turbulence
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
path_to_reference_result
#=====================================================
def reinit_inputs():
    global subdirectories, filenames, labels, black_line_flag, \
    dashed_line_flag, figure_filename_postfix, figure_title, \
    ylimits_kinetic_energy_input, ylimits_dissipation_input, \
    log_axes_input, legend_on_input, legend_inside_input, \
    plot_reference_result, nlegendcols_input, \
    figure_subdirectory, data_directory_base, figure_directory_base, \
    path_to_reference_result

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
    plot_reference_result=True # default
    nlegendcols_input=1
    figure_subdirectory="" # default
    data_directory_base = "/Users/Julien/julien_phd/post_processing/data/taylor_green_vortex"
    path_to_reference_result = "/Users/Julien/julien_phd/post_processing/data/taylor_green_vortex/dns"
    # figure_directory_base = "/Users/Julien/julien_phd/post_processing/figures/taylor_green_vortex"
    figure_directory_base = "figures"
#=====================================================
#-----------------------------------------------------

#=====================================================
# DOFs: 96^3 | NSFR.IR-cDG-GL and NSFR.IR-cDG-GLL (BASELINE)
#-----------------------------------------------------
if(True):
    reinit_inputs()
    data_directory_base="/Users/Julien/NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL$=0.10$" # comment to turn off
    figure_filename_postfix = "96_baseline_scheme"
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
        path_to_reference_result,
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
# DOFs: 96^3 | NSFR.IR-cDG (GL and GLL) WITH OVERINTEGRATION
#-----------------------------------------------------
if(True):
    reinit_inputs()
    data_directory_base="/Users/Julien/NarvalFiles/2023_JCP/"
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
        path_to_reference_result,
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
if(True):
    reinit_inputs()
    data_directory_base="/Users/Julien/NarvalFiles/2023_JCP/"
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
        path_to_reference_result,
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
# DOFs: 96^3 | Basic SGS Models on GLL flux nodes (no filter width modifications)
#-----------------------------------------------------
if(True):
    # batch_labels = [ \
    # "$c_{DG}$ NSFR", \
    # "$c_{DG}$ NSFR.IR-Smag. $C_{S}=0.10$", \
    # "$c_{DG}$ NSFR.IR-Smag. $C_{S}=0.18$", \
    # "$c_{DG}$ NSFR.IR-WALE", \
    # "$c_{DG}$ NSFR.IR-Vreman", \
    # "$c_{DG}$ NSFR.IR-Smag. $C_{S}=0.1$, $\\Delta_{min}$", \
    # ]
    reinit_inputs()
    data_directory_base="/Users/Julien/NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL$=0.10$" # comment to turn off
    figure_filename_postfix = "96_sgs_models_gll"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories.append("flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GLL")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("sgs_model_GLL_flux_nodes/viscous_TGV_LES_SMAG_MC-0.10_NSFR_cDG_IR_2PF_GLL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GLL-Smag. $C_{S}=0.10$")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("sgs_model_GLL_flux_nodes/viscous_TGV_LES_SMAG_MC-0.18_NSFR_cDG_IR_2PF_GLL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GLL-Smag. $C_{S}=0.18$")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("sgs_model_GLL_flux_nodes/viscous_TGV_LES_WALE_MC-0.50_NSFR_cDG_IR_2PF_GLL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GLL-WALE $C_{W}=0.50$")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("sgs_model_GLL_flux_nodes/viscous_TGV_LES_VRMN_MC-0.081_NSFR_cDG_IR_2PF_GLL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GLL-Vreman $C_{V}=0.081$")
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
        path_to_reference_result,
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
# DOFs: 96^3 | Basic SGS Models on GL flux nodes (no filter width modifications)
#-----------------------------------------------------
if(True):
    # batch_labels = [ \
    # "$c_{DG}$ NSFR", \
    # "$c_{DG}$ NSFR.IR-Smag. $C_{S}=0.10$", \
    # "$c_{DG}$ NSFR.IR-Smag. $C_{S}=0.18$", \
    # "$c_{DG}$ NSFR.IR-WALE", \
    # "$c_{DG}$ NSFR.IR-Vreman", \
    # "$c_{DG}$ NSFR.IR-Smag. $C_{S}=0.1$, $\\Delta_{min}$", \
    # ]
    reinit_inputs()
    data_directory_base="/Users/Julien/NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL$=0.10$" # comment to turn off
    figure_filename_postfix = "96_sgs_models_gl"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories.append("flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GL")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("sgs_model_GL_flux_nodes/viscous_TGV_LES_SMAG_MC-0.10_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GL-Smag. $C_{S}=0.10$")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("sgs_model_GL_flux_nodes/viscous_TGV_LES_SMAG_MC-0.18_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GL-Smag. $C_{S}=0.18$")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("sgs_model_GL_flux_nodes/viscous_TGV_LES_WALE_MC-0.50_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GL-WALE. $C_{W}=0.50$")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("sgs_model_GL_flux_nodes/viscous_TGV_LES_VRMN_MC-0.081_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GL-Vreman. $C_{V}=0.081$")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("sgs_model_GL_flux_nodes/viscous_TGV_LES_SMAG_MC-0.18_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512_filter_width_flad_and_gassner")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GL-Smag. $C_{S}=0.18$ $\\Delta=\\frac{V}{(P+1)^{3}}$")
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
        path_to_reference_result,
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
# DOFs: 96^3 | Upwind dissipation on GL flux nodes
#-----------------------------------------------------
if(True):
    reinit_inputs()
    data_directory_base="/Users/Julien/NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL$=0.10$" # comment to turn off
    figure_filename_postfix = "96_upwind_gl"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories.append("flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GL")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("upwind_dissipation_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF-LxF_GL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GL-LxF")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("upwind_dissipation_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF-Roe_GL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GL-Roe")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("upwind_dissipation_GL_flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF-L2R_GL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GL-L2R")
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
        path_to_reference_result,
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
if(True):
    reinit_inputs()
    data_directory_base="/Users/Julien/NarvalFiles/2023_JCP/"
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
        path_to_reference_result,
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
if(True):
    reinit_inputs()
    data_directory_base="/Users/Julien/NarvalFiles/2023_JCP/"
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
        path_to_reference_result,
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
if(True):
    reinit_inputs()
    data_directory_base="/Users/Julien/NarvalFiles/2023_JCP/"
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
        path_to_reference_result,
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
# if(True):
#     reinit_inputs()
#     data_directory_base="/Users/Julien/NarvalFiles/2023_JCP/"
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
#         path_to_reference_result,
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
if(True):
    reinit_inputs()
    data_directory_base="/Users/Julien/NarvalFiles/2023_JCP/"
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
        path_to_reference_result,
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
# DOFs: 96^3 | NSFR Two-Point-Flux
#-----------------------------------------------------
if(True):
    reinit_inputs()
    data_directory_base="/Users/Julien/NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, $96^{3}$ DOFs, CFL$=0.10$" # comment to turn off
    figure_filename_postfix = "96_two_point_flux"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories.append("flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.IR-GL")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("two_point_flux/viscous_TGV_ILES_NSFR_cDG_KG_2PF_GL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.KG-GL")
    black_line_flag.append(False)
    dashed_line_flag.append(True)
    #-----------------------------------------------------
    subdirectories.append("two_point_flux/viscous_TGV_ILES_NSFR_cDG_CH_2PF_GL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.CH-GL")
    black_line_flag.append(False)
    dashed_line_flag.append(True)
    #-----------------------------------------------------
    subdirectories.append("two_point_flux/viscous_TGV_ILES_NSFR_cDG_Ra_2PF_GL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$c_{DG}$ NSFR.CH$_{\\mathrm{RA}}$-GL")
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
        path_to_reference_result,
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
# DOFs: ALL | NSFR CONVERGENCE
#-----------------------------------------------------
if(True):
    reinit_inputs()
    data_directory_base="/Users/Julien/NarvalFiles/2023_JCP/"
    date_for_runs="."
    figure_subdirectory="2023_JCP"
    figure_title = "TGV at Re$_{\\infty}=1600$, P$5$, CFL$=0.10$" # comment to turn off
    figure_filename_postfix = "cDG_NSFR_convergence"
    legend_inside_input=True
    #-----------------------------------------------------
    subdirectories.append("flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512")
    filenames.append("turbulent_quantities.txt")
    labels.append("$96^{3}$ DOFs, $c_{DG}$ NSFR.IR-GL")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("robustness/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs048_p5_procs64")
    filenames.append("turbulent_quantities.txt")
    labels.append("$48^{3}$ DOFs, $c_{DG}$ NSFR.IR-GL")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("robustness/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs024_p5_procs16")
    filenames.append("turbulent_quantities.txt")
    labels.append("$24^{3}$ DOFs, $c_{DG}$ NSFR.IR-GL")
    black_line_flag.append(False)
    dashed_line_flag.append(False)
    #-----------------------------------------------------
    subdirectories.append("robustness/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs012_p5_procs16")
    filenames.append("turbulent_quantities.txt")
    labels.append("$12^{3}$ DOFs, $c_{DG}$ NSFR.IR-GL")
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
        path_to_reference_result,
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
# DOFs: ALL | Strong DG CONVERGENCE
#-----------------------------------------------------
if(True):
    reinit_inputs()
    data_directory_base="/Users/Julien/NarvalFiles/2023_JCP/"
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
        path_to_reference_result,
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