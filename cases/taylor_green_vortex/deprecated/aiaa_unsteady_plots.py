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
    figure_subdirectory=""
    data_directory_base = "/Users/Julien/julien_phd/post_processing/data/taylor_green_vortex"
    path_to_reference_result = "/Users/Julien/julien_phd/post_processing/data/taylor_green_vortex/dns"
    # figure_directory_base = "/Users/Julien/julien_phd/post_processing/figures/taylor_green_vortex"
    figure_directory_base = "figures"
#=====================================================
#=====================================================
# DOFs: 96^3 | NSFR-Roe-Smag.010
#-----------------------------------------------------
reinit_inputs()
data_directory_base="/Users/Julien/NarvalFiles/2023_AIAA/"
date_for_runs="."
figure_subdirectory=date_for_runs
figure_title = "TGV at $Re_{\\infty}=1600$, $P5$, CFL $=0.10$"
figure_title = " "
figure_filename_postfix = "96_cDG_flux_nodes_vs_standard_DG_update"
legend_inside_input=True
#-----------------------------------------------------
subdirectories.append("2022-11-29_TGV_SPECTRA_96dofs/viscous_TGV_ILES_std_strong_DG_roe_dofs096_p5_procs512_timed")
filenames.append("turbulent_quantities.txt")
labels.append("Strong DG, Roe, uncoll.,\n $n_{quad}=2(P+1)$, $CFL=0.10$")
black_line_flag.append(False)
dashed_line_flag.append(False)
#-----------------------------------------------------
subdirectories.append("2022-11-29_TGV_SPECTRA_96dofs/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs096_p5_procs512_timed")
filenames.append("turbulent_quantities.txt")
labels.append("$c_{DG}$ NSFR coll., $CFL=0.10$")
black_line_flag.append(False)
dashed_line_flag.append(False)
# #-----------------------------------------------------
# subdirectories.append("2022-11-29_TGV_SPECTRA_96dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs096_p5_procs512_timed")
# filenames.append("turbulent_quantities.txt")
# labels.append("$c_{+}$ NSFR coll.")
# black_line_flag.append(False)
# dashed_line_flag.append(False)
# #-----------------------------------------------------
# subdirectories.append("2022-11-29_TGV_SPECTRA_96dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs096_p5_procs512_uncollocated")
# filenames.append("turbulent_quantities.txt")
# labels.append("$c_{+}$ NSFR uncoll.")
# black_line_flag.append(False)
# dashed_line_flag.append(False)
#-----------------------------------------------------
subdirectories.append("2022-11-29_TGV_SPECTRA_96dofs/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs096_p5_procs512_uncollocated")
filenames.append("turbulent_quantities-fail.txt")
labels.append("$c_{DG}$ NSFR uncoll.\n [without GL fix], $CFL=0.10$")
black_line_flag.append(False)
dashed_line_flag.append(False)
#-----------------------------------------------------
subdirectories.append("2022-11-29_TGV_SPECTRA_96dofs/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs096_p5_procs512_uncollocated")
filenames.append("turbulent_quantities.txt")
labels.append("$c_{DG}$ NSFR uncoll.\n [with GL fix], $CFL=0.10$")
black_line_flag.append(False)
dashed_line_flag.append(False)
#-----------------------------------------------------
subdirectories.append("2022-11-29_TGV_SPECTRA_96dofs/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs096_p5_procs512_uncollocated")
filenames.append("turbulent_quantities_before_gl_fix.txt")
labels.append("$c_{DG}$ NSFR uncoll.\n [without GL fix], $CFL=0.05$")
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
    tmax=10.0,
    legend_fontSize_input=14,
    solid_and_dashed_lines=False)
#-----------------------------------------------------
exit()

#=====================================================
# DOFs: 48^3 | NSFR-Roe-Smag.010
#-----------------------------------------------------
reinit_inputs()
data_directory_base="/Users/Julien/NarvalFiles/2023_AIAA/"
date_for_runs="."
figure_subdirectory=date_for_runs
figure_title = "TGV at $Re_{\\infty}=1600$, $P5$, CFL $=0.1$"
figure_title = " "
figure_filename_postfix = "96_cDG_vs_cPlus_vs_standard_DG"
legend_inside_input=True
#-----------------------------------------------------
subdirectories.append("2022-11-29_TGV_SPECTRA_96dofs/viscous_TGV_ILES_std_strong_DG_lax_friedrichs_dofs096_p5_procs512_timed")
filenames.append("turbulent_quantities.txt")
labels.append("Strong DG with $n_{quad}=2(P+1)$")
black_line_flag.append(False)
dashed_line_flag.append(False)
#-----------------------------------------------------
subdirectories.append("2022-11-29_TGV_SPECTRA_96dofs/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs096_p5_procs512_timed")
filenames.append("turbulent_quantities.txt")
labels.append("$c_{DG}$ NSFR")
black_line_flag.append(False)
dashed_line_flag.append(False)
#-----------------------------------------------------
subdirectories.append("2022-11-29_TGV_SPECTRA_96dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs096_p5_procs512_timed")
filenames.append("turbulent_quantities.txt")
labels.append("$c_{+}$ NSFR")
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
exit()

#=====================================================
# DOFs: 48^3 | NSFR-Roe-Smag.010
#-----------------------------------------------------
reinit_inputs()
data_directory_base="/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/"
date_for_runs="."
figure_subdirectory=date_for_runs
figure_title = "TGV at $Re_{\\infty}=1600$, $P5$, CFL $=0.1$, $c_{+}$ NSFR (IR)"
figure_title = " "
figure_filename_postfix = "48_cPlus_vs_cDG_SGS_models_with_roe"
legend_inside_input=True
#-----------------------------------------------------
subdirectories.append("viscous_TGV_ILES_cDG_IR_two_point_flux_dofs048_p5_procs64")
filenames.append("turbulent_quantities.txt")
labels.append("$c_{DG}$ NSFR")
black_line_flag.append(False)
dashed_line_flag.append(False)
#-----------------------------------------------------
subdirectories.append("viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_with_roe_dissipation_dofs048_p5_procs64_filter36timeslarger")
filenames.append("turbulent_quantities.txt")
labels.append("$c_{DG}$ NSFR-Roe-Smag.010")
black_line_flag.append(False)
dashed_line_flag.append(False)
#-----------------------------------------------------
subdirectories.append("viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs048_p5_procs64")
filenames.append("turbulent_quantities.txt")
labels.append("$c_{+}$ NSFR")
black_line_flag.append(False)
dashed_line_flag.append(False)
#-----------------------------------------------------
subdirectories.append("viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_with_roe_dissipation_dofs048_p5_procs64_filter36timeslarger")
filenames.append("turbulent_quantities.txt")
labels.append("$c_{+}$ NSFR-Roe-Smag.010")
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
    tmax=10.0,
    legend_fontSize_input=14,
    solid_and_dashed_lines=True)
#-----------------------------------------------------
exit()

#=====================================================
# DOFs: 24^3, 48^3, 96^3 | c: cPlus
#-----------------------------------------------------
reinit_inputs()
data_directory_base="/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/"
date_for_runs="2022-12-05"
figure_subdirectory=date_for_runs
figure_title = "TGV at $Re_{\\infty}=1600$, $P5$, CFL $=0.1$, $c_{+}$ NSFR (IR)"
figure_title = " "
figure_filename_postfix = "48_cPlus_l2r_sgs_tgv"
legend_inside_input=True
#-----------------------------------------------------
subdirectories.append("viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs048_p5_procs64")
filenames.append("turbulent_quantities.txt")
labels.append("$c_{+}$ NSFR")
black_line_flag.append(False)
dashed_line_flag.append(False)
#-----------------------------------------------------
subdirectories.append("viscous_TGV_ILES_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64")
filenames.append("turbulent_quantities.txt")
labels.append("$c_{+}$ NSFR-L2R")
black_line_flag.append(False)
dashed_line_flag.append(False)
#-----------------------------------------------------
subdirectories.append("viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_dofs048_p5_procs64_filter36timeslarger")
filenames.append("turbulent_quantities.txt")
labels.append("$c_{+}$ NSFR-SGS")
black_line_flag.append(False)
dashed_line_flag.append(False)
#-----------------------------------------------------
subdirectories.append("viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64_filter36timeslarger")
filenames.append("turbulent_quantities.txt")
labels.append("$c_{+}$ NSFR-L2R-SGS")
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
    tmax=10.0,
    legend_fontSize_input=14)
#-----------------------------------------------------
#=====================================================
# DOFs: 24^3, 48^3, 96^3 | c: cPlus
#-----------------------------------------------------
reinit_inputs()
data_directory_base="/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/"
date_for_runs="2022-12-05"
figure_subdirectory=date_for_runs
figure_title = "TGV at $Re_{\\infty}=1600$, $P5$, CFL $=0.1$, $c_{+}$ NSFR (IR)"
figure_title = " "
figure_filename_postfix = "48_cDG_l2r_sgs_tgv"
legend_inside_input=True
#-----------------------------------------------------
subdirectories.append("viscous_TGV_ILES_cDG_IR_two_point_flux_dofs048_p5_procs64")
filenames.append("turbulent_quantities.txt")
labels.append("$c_{DG}$ NSFR")
black_line_flag.append(False)
dashed_line_flag.append(False)
#-----------------------------------------------------
subdirectories.append("viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_dofs048_p5_procs64_filter36timeslarger")
filenames.append("turbulent_quantities.txt")
labels.append("$c_{DG}$ NSFR-SGS")
black_line_flag.append(False)
dashed_line_flag.append(False)
#-----------------------------------------------------
subdirectories.append("viscous_TGV_ILES_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64")
filenames.append("turbulent_quantities.txt")
labels.append("$c_{DG}$ NSFR-L2R")
black_line_flag.append(False)
dashed_line_flag.append(False)
#-----------------------------------------------------
subdirectories.append("viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64_filter36timeslarger")
filenames.append("turbulent_quantities.txt")
labels.append("$c_{DG}$ NSFR-L2R-SGS")
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
    tmax=10.0,
    legend_fontSize_input=14,
    solid_and_dashed_lines=False)
#-----------------------------------------------------
exit()

# #=====================================================
# # TP-Flux: KG | SGS: None | c: All | Riemann: none
# #-----------------------------------------------------
# reinit_inputs()
# date_for_runs="2022-09-06"
# figure_subdirectory=date_for_runs+"/"+"tpflux-KG_sgs-none_c-all_riemann-none"
# figure_title = "TGV at $Re_{\\infty}=1600$, $96^{3}$ DOFs ($P5$, $n_{el}=16^{3}$), CFL $=0.1$"#"Viscous TGV, $64^{3}$ DOFs, P3, CFL$=0.025$"
# figure_filename_postfix = "0"
# legend_inside_input=True
# #-----------------------------------------------------
# subdirectories.append(date_for_runs+"/"+"viscous_ILES_cDG_KG_two_point_flux_dofs096_p5")
# filenames.append("turbulent_quantities.txt")
# labels.append("$c_{DG}$ NSFR")
# black_line_flag.append(False)
# dashed_line_flag.append(False)
# #-----------------------------------------------------
# subdirectories.append(date_for_runs+"/"+"viscous_ILES_cPlus_KG_two_point_flux_dofs096_p5")
# filenames.append("turbulent_quantities.txt")
# labels.append("$c_{+}$ NSFR")
# black_line_flag.append(False)
# dashed_line_flag.append(False)
# #-----------------------------------------------------
# subdirectories.append(date_for_runs+"/"+"viscous_ILES_cSD_KG_two_point_flux_dofs096_p5")
# filenames.append("turbulent_quantities.txt")
# labels.append("$c_{SD}$ NSFR")
# black_line_flag.append(False)
# dashed_line_flag.append(False)
# #-----------------------------------------------------
# subdirectories.append(date_for_runs+"/"+"viscous_ILES_cHU_KG_two_point_flux_dofs096_p5")
# filenames.append("turbulent_quantities.txt")
# labels.append("$c_{HU}$ NSFR")
# black_line_flag.append(False)
# dashed_line_flag.append(False)
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
#     path_to_reference_result,
#     figure_filename_postfix,
#     figure_title,
#     log_axes_input,
#     legend_on_input,
#     legend_inside_input,
#     nlegendcols_input)
# #-----------------------------------------------------
# #=====================================================

# #=====================================================
# # TP-Flux: KG | SGS: None | c: All | Riemann: L2Roe
# #-----------------------------------------------------
# reinit_inputs()
# date_for_runs="2022-09-06"
# figure_subdirectory=date_for_runs+"/"+"tpflux-KG_sgs-none_c-all_riemann-l2roe"
# figure_title = "TGV at $Re_{\\infty}=1600$, $96^{3}$ DOFs ($P5$, $n_{el}=16^{3}$), CFL $=0.1$"#"Viscous TGV, $64^{3}$ DOFs, P3, CFL$=0.025$"
# figure_filename_postfix = "1"
# legend_inside_input=True
# #-----------------------------------------------------
# subdirectories.append(date_for_runs+"/"+"viscous_ILES_cDG_KG_two_point_flux_with_l2roe_dissipation_dofs096_p5")
# filenames.append("turbulent_quantities.txt")
# labels.append("$c_{DG}$ NSFR")
# black_line_flag.append(False)
# dashed_line_flag.append(False)
# #-----------------------------------------------------
# subdirectories.append(date_for_runs+"/"+"viscous_ILES_cPlus_KG_two_point_flux_with_l2roe_dissipation_dofs096_p5")
# filenames.append("turbulent_quantities.txt")
# labels.append("$c_{+}$ NSFR")
# black_line_flag.append(False)
# dashed_line_flag.append(False)
# #-----------------------------------------------------
# subdirectories.append(date_for_runs+"/"+"viscous_ILES_cSD_KG_two_point_flux_with_l2roe_dissipation_dofs096_p5")
# filenames.append("turbulent_quantities.txt")
# labels.append("$c_{SD}$ NSFR")
# black_line_flag.append(False)
# dashed_line_flag.append(False)
# #-----------------------------------------------------
# subdirectories.append(date_for_runs+"/"+"viscous_ILES_cHU_KG_two_point_flux_with_l2roe_dissipation_dofs096_p5")
# filenames.append("turbulent_quantities.txt")
# labels.append("$c_{HU}$ NSFR")
# black_line_flag.append(False)
# dashed_line_flag.append(False)
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
#     path_to_reference_result,
#     figure_filename_postfix,
#     figure_title,
#     log_axes_input,
#     legend_on_input,
#     legend_inside_input,
#     nlegendcols_input)
# #-----------------------------------------------------
# #=====================================================

# #=====================================================
# # TP-Flux: KG | SGS: Smag. | c: All | Riemann: none
# #-----------------------------------------------------
# reinit_inputs()
# date_for_runs="2022-09-06"
# figure_subdirectory=date_for_runs+"/"+"tpflux-KG_sgs-smag_c-all_riemann-none"
# figure_title = "TGV at $Re_{\\infty}=1600$, $96^{3}$ DOFs ($P5$, $n_{el}=16^{3}$), CFL $=0.1$"#"Viscous TGV, $64^{3}$ DOFs, P3, CFL$=0.025$"
# figure_filename_postfix = "2"
# legend_inside_input=True
# #-----------------------------------------------------
# subdirectories.append(date_for_runs+"/"+"viscous_LES_smagorinsky_cDG_KG_two_point_flux_dofs096_p5")
# filenames.append("turbulent_quantities.txt")
# labels.append("$c_{DG}$ NSFR")
# black_line_flag.append(False)
# dashed_line_flag.append(False)
# #-----------------------------------------------------
# subdirectories.append(date_for_runs+"/"+"viscous_LES_smagorinsky_cPlus_KG_two_point_flux_dofs096_p5")
# filenames.append("turbulent_quantities.txt")
# labels.append("$c_{+}$ NSFR")
# black_line_flag.append(False)
# dashed_line_flag.append(False)
# #-----------------------------------------------------
# subdirectories.append(date_for_runs+"/"+"viscous_LES_smagorinsky_cSD_KG_two_point_flux_dofs096_p5")
# filenames.append("turbulent_quantities.txt")
# labels.append("$c_{SD}$ NSFR")
# black_line_flag.append(False)
# dashed_line_flag.append(False)
# #-----------------------------------------------------
# subdirectories.append(date_for_runs+"/"+"viscous_LES_smagorinsky_cHU_KG_two_point_flux_dofs096_p5")
# filenames.append("turbulent_quantities.txt")
# labels.append("$c_{HU}$ NSFR")
# black_line_flag.append(False)
# dashed_line_flag.append(False)
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
#     path_to_reference_result,
#     figure_filename_postfix,
#     figure_title,
#     log_axes_input,
#     legend_on_input,
#     legend_inside_input,
#     nlegendcols_input)
# #-----------------------------------------------------
# #=====================================================

# #=====================================================
# # TP-Flux: KG | SGS: Smag. | c: All | Riemann: L2Roe
# #-----------------------------------------------------
# reinit_inputs()
# date_for_runs="2022-09-06"
# figure_subdirectory=date_for_runs+"/"+"tpflux-KG_sgs-smag_c-all_riemann-l2roe"
# figure_title = "TGV at $Re_{\\infty}=1600$, $96^{3}$ DOFs ($P5$, $n_{el}=16^{3}$), CFL $=0.1$"#"Viscous TGV, $64^{3}$ DOFs, P3, CFL$=0.025$"
# figure_filename_postfix = "3"
# legend_inside_input=True
# #-----------------------------------------------------
# subdirectories.append(date_for_runs+"/"+"viscous_LES_smagorinsky_cDG_KG_two_point_flux_with_l2roe_dissipation_dofs096_p5")
# filenames.append("turbulent_quantities.txt")
# labels.append("$c_{DG}$ NSFR")
# black_line_flag.append(False)
# dashed_line_flag.append(False)
# #-----------------------------------------------------
# subdirectories.append(date_for_runs+"/"+"viscous_LES_smagorinsky_cPlus_KG_two_point_flux_with_l2roe_dissipation_dofs096_p5")
# filenames.append("turbulent_quantities.txt")
# labels.append("$c_{+}$ NSFR")
# black_line_flag.append(False)
# dashed_line_flag.append(False)
# #-----------------------------------------------------
# subdirectories.append(date_for_runs+"/"+"viscous_LES_smagorinsky_cSD_KG_two_point_flux_with_l2roe_dissipation_dofs096_p5")
# filenames.append("turbulent_quantities.txt")
# labels.append("$c_{SD}$ NSFR")
# black_line_flag.append(False)
# dashed_line_flag.append(False)
# #-----------------------------------------------------
# subdirectories.append(date_for_runs+"/"+"viscous_LES_smagorinsky_cHU_KG_two_point_flux_with_l2roe_dissipation_dofs096_p5")
# filenames.append("turbulent_quantities.txt")
# labels.append("$c_{HU}$ NSFR")
# black_line_flag.append(False)
# dashed_line_flag.append(False)
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
#     path_to_reference_result,
#     figure_filename_postfix,
#     figure_title,
#     log_axes_input,
#     legend_on_input,
#     legend_inside_input,
#     nlegendcols_input)
# #-----------------------------------------------------
# #=====================================================

# #=====================================================
# # TP-Flux: KG | SGS: None | c: cDG | Riemann: L2Roe
# # DOFs: 96^3, 256^3
# #-----------------------------------------------------
# reinit_inputs()
# date_for_runs="2022-09-06"
# figure_subdirectory=date_for_runs+"/"+"convergence"
# figure_title = "TGV at $Re_{\\infty}=1600$, $c_{DG}$ NSFR, CFL $=0.1$"#"Viscous TGV, $64^{3}$ DOFs, P3, CFL$=0.025$"
# figure_filename_postfix = "4"
# legend_inside_input=True
# #-----------------------------------------------------
# subdirectories.append(date_for_runs+"/"+"viscous_ILES_cDG_KG_two_point_flux_with_l2roe_dissipation_dofs096_p5")
# filenames.append("turbulent_quantities.txt")
# labels.append("$96^{3}$ DOFs ($P5$, $n_{el}=16^{3}$)")
# black_line_flag.append(False)
# dashed_line_flag.append(False)
# #-----------------------------------------------------
# subdirectories.append(date_for_runs+"/"+"viscous_ILES_cDG_KG_two_point_flux_with_l2roe_dissipation_dofs0256_p3_procs4096")
# filenames.append("turbulent_quantities.txt")
# labels.append("$256^{3}$ DOFs ($P3$, $n_{el}=64^{3}$)")
# black_line_flag.append(False)
# dashed_line_flag.append(False)
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
#     path_to_reference_result,
#     figure_filename_postfix,
#     figure_title,
#     log_axes_input,
#     legend_on_input,
#     legend_inside_input,
#     nlegendcols_input)
# #-----------------------------------------------------
# #=====================================================

# #=====================================================
# # TP-Flux: KG | SGS: None | c: cDG | Riemann: L2Roe
# # DG: Weak and Strong with L2Roe + 0.5CFL KG+L2Roe
# # DOFs: 96^3
# #-----------------------------------------------------
# reinit_inputs()
# date_for_runs="2022-09-06"
# figure_subdirectory=date_for_runs+"/"+"strong_vs_weak"
# figure_title = "TGV at $Re_{\\infty}=1600$, $96^{3}$ DOFs ($P5$, $n_{el}=16^{3}$), CFL $=0.1$"#"Viscous TGV, $64^{3}$ DOFs, P3, CFL$=0.025$"
# figure_filename_postfix = "5"
# legend_inside_input=True
# #-----------------------------------------------------
# date_for_runs="2022-09-06"
# subdirectories.append(date_for_runs+"/"+"viscous_ILES_cDG_KG_two_point_flux_with_l2roe_dissipation_dofs096_p5")
# filenames.append("turbulent_quantities.txt")
# labels.append("$c_{DG}$ NSFR (KG + L$^{2}$Roe)\n CFL=0.1")
# black_line_flag.append(False)
# dashed_line_flag.append(False)
# #-----------------------------------------------------
# date_for_runs="2022-09-06"
# subdirectories.append(date_for_runs+"/"+"viscous_ILES_cDG_KG_two_point_flux_with_l2roe_dissipation_dofs0256_p3_procs4096")
# filenames.append("turbulent_quantities.txt")
# labels.append("$256^{3}$ DOFs ($P3$, $n_{el}=64^{3}$)")
# black_line_flag.append(False)
# dashed_line_flag.append(False)
# #-----------------------------------------------------
# date_for_runs="2022-10-04"
# subdirectories.append(date_for_runs+"/"+"viscous_ILES_strong_DG_l2roe_dofs096_p5")
# filenames.append("turbulent_quantities.txt")
# labels.append("Strong DG (L$^{2}$Roe)\n CFL=0.1")
# black_line_flag.append(False)
# dashed_line_flag.append(False)
# #-----------------------------------------------------
# date_for_runs="2022-10-04"
# subdirectories.append(date_for_runs+"/"+"viscous_ILES_weak_DG_l2roe_dofs096_p5")
# filenames.append("turbulent_quantities.txt")
# labels.append("Weak DG (L$^{2}$Roe)\n CFL=0.1")
# black_line_flag.append(False)
# dashed_line_flag.append(False)
#-----------------------------------------------------
# date_for_runs="2022-10-04"
# subdirectories.append(date_for_runs+"/"+"viscous_ILES_weak_DG_l2roe_dofs096_p5_master")
# filenames.append("turbulent_quantities.txt")
# labels.append("Weak DG (L$^{2}$Roe)\n CFL=0.1\n \\textbf{[master Branch]}")
# black_line_flag.append(False)
# dashed_line_flag.append(False)
# #-----------------------------------------------------
# date_for_runs="2022-10-04"
# subdirectories.append(date_for_runs+"/"+"viscous_ILES_weak_DG_l2roe_dofs096_p5_master_higher_penalty")
# filenames.append("turbulent_quantities.txt")
# labels.append("Weak DG (L$^{2}$Roe)\n CFL=0.1, C_{IP}=10\n \\textbf{[master Branch]}")
# black_line_flag.append(False)
# dashed_line_flag.append(False)
# #-----------------------------------------------------
# date_for_runs="2022-10-04"
# subdirectories.append(date_for_runs+"/"+"viscous_ILES_cDG_KG_two_point_flux_with_l2roe_dissipation_dofs096_p5_projected_initial_condition")
# filenames.append("turbulent_quantities.txt")
# labels.append("$c_{DG}$ NSFR (KG + L$^{2}$Roe)\n CFL=0.1, Projected IC")
# black_line_flag.append(False)
# dashed_line_flag.append(False)
# #-----------------------------------------------------
# date_for_runs="2022-10-04"
# subdirectories.append(date_for_runs+"/"+"viscous_ILES_cDG_KG_two_point_flux_with_l2roe_dissipation_dofs096_p5_uniform_density")
# filenames.append("turbulent_quantities.txt")
# labels.append("$c_{DG}$ NSFR (KG + L$^{2}$Roe)\n CFL=0.1, $\\rho_{0}=1$")
# black_line_flag.append(False)
# dashed_line_flag.append(False)
# #-----------------------------------------------------
# date_for_runs="2022-10-04"
# subdirectories.append(date_for_runs+"/"+"viscous_ILES_cDG_KG_two_point_flux_dofs064_p3_uniform_density_low_cfl")
# filenames.append("turbulent_quantities.txt")
# labels.append("$c_{DG}$ NSFR (KG), CFL=0.015\n $\\rho_{0}=1$, $\\Delta t=const.$ \n P3, $64^{3}$ DOFS, SSP-RK3")
# black_line_flag.append(False)
# dashed_line_flag.append(False)
# #-----------------------------------------------------
# date_for_runs="2022-10-04"
# subdirectories.append(date_for_runs+"/"+"viscous_ILES_cDG_KG_two_point_flux_with_l2roe_dissipation_dofs096_p5_constant_time_step")
# filenames.append("turbulent_quantities.txt")
# labels.append("$c_{DG}$ NSFR (KG), CFL=0.015\n $\\Delta t=const.$, RK4")
# black_line_flag.append(False)
# dashed_line_flag.append(False)
# #-----------------------------------------------------
# date_for_runs="2022-10-04"
# subdirectories.append(date_for_runs+"/"+"viscous_ILES_cDG_KG_two_point_with_l2roe_dissipation_dofs096_p5_adaptive_time_step_mpi_fix")
# filenames.append("turbulent_quantities.txt")
# labels.append("$c_{DG}$ NSFR (KG + L$^{2}$Roe), CFL=0.1\n $\\Delta t\\neq const.$, RK4\n \\textbf{MPI bug fix}")
# black_line_flag.append(False)
# dashed_line_flag.append(False)
# -----------------------------------------------------
# date_for_runs="2022-10-04"
# subdirectories.append(date_for_runs+"/"+"viscous_ILES_weak_DG_l2roe_br2_dofs096_p5_master")
# filenames.append("turbulent_quantities.txt")
# labels.append("Weak DG (L$^{2}$Roe+BR2)\n CFL=0.1\n \\textbf{[master Branch]}")
# black_line_flag.append(False)
# dashed_line_flag.append(False)
#-----------------------------------------------------
# date_for_runs="2022-10-04"
# subdirectories.append(date_for_runs+"/"+"viscous_ILES_weak_DG_lxf_dofs096_p5_master")
# filenames.append("turbulent_quantities.txt")
# labels.append("Weak DG (LxF)\n CFL=0.1\n \\textbf{[master Branch]}")
# black_line_flag.append(False)
# dashed_line_flag.append(False)
# #-----------------------------------------------------
# date_for_runs="2022-10-04"
# subdirectories.append(date_for_runs+"/"+"viscous_ILES_cDG_KG_two_point_flux_with_l2roe_dissipation_dofs096_p5")
# filenames.append("turbulent_quantities.txt")
# labels.append("$c_{DG}$ NSFR (KG + L$^{2}$Roe)\n CFL=0.05")
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
#     path_to_reference_result,
#     figure_filename_postfix,
#     figure_title,
#     log_axes_input,
#     legend_on_input,
#     legend_inside_input,
#     nlegendcols_input,
#     transparent_legend_input=False)
# #-----------------------------------------------------
# #=====================================================

#=====================================================
# DOFs: 96^3
#-----------------------------------------------------
reinit_inputs()
date_for_runs="2022-11-11"
figure_subdirectory=date_for_runs
figure_title = "TGV at $Re_{\\infty}=1600$, $96^{3}$ DOFs ($P5$, $n_{el}=16^{3}$), CFL $=0.1$"
figure_title = " "
figure_filename_postfix = "dofs096_tgv"
legend_inside_input=True
#-----------------------------------------------------
date_for_runs="2022-11-09_96dofs"
subdirectories.append(date_for_runs+"/"+"viscous_TGV_ILES_std_strong_DG_lax_friedrichs_dofs096_p5_procs512")
filenames.append("turbulent_quantities.txt")
labels.append("Strong DG")
black_line_flag.append(False)
dashed_line_flag.append(False)
#-----------------------------------------------------
date_for_runs="2022-11-09_96dofs"
subdirectories.append(date_for_runs+"/"+"viscous_TGV_ILES_cDG_IR_two_point_flux_dofs096_p5_procs512")
filenames.append("turbulent_quantities.txt")
labels.append("$c_{DG}$ NSFR")
black_line_flag.append(False)
dashed_line_flag.append(False)
#-----------------------------------------------------
date_for_runs="2022-11-09_96dofs"
subdirectories.append(date_for_runs+"/"+"viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs096_p5_procs512")
filenames.append("turbulent_quantities.txt")
labels.append("$c_{+}$ NSFR")
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
    legend_fontSize_input=14)
#-----------------------------------------------------
#=====================================================

#=====================================================
# DOFs: 48^3
#-----------------------------------------------------
reinit_inputs()
date_for_runs="2022-11-11"
figure_subdirectory=date_for_runs
figure_title = "TGV at $Re_{\\infty}=1600$, $48^{3}$ DOFs ($P5$, $n_{el}=8^{3}$), CFL $=0.1$"
figure_title = " "
figure_filename_postfix = "dofs048_tgv"
legend_inside_input=True
#-----------------------------------------------------
date_for_runs="2022-11-09_48dofs"
subdirectories.append(date_for_runs+"/"+"viscous_TGV_ILES_std_strong_DG_lax_friedrichs_dofs048_p5_procs64")
filenames.append("turbulent_quantities.txt")
labels.append("Strong DG")
black_line_flag.append(False)
dashed_line_flag.append(False)
#-----------------------------------------------------
date_for_runs="2022-11-09_48dofs"
subdirectories.append(date_for_runs+"/"+"viscous_TGV_ILES_cDG_IR_two_point_flux_dofs048_p5_procs64")
filenames.append("turbulent_quantities.txt")
labels.append("$c_{DG}$ NSFR")
black_line_flag.append(False)
dashed_line_flag.append(False)
#-----------------------------------------------------
date_for_runs="2022-11-09_48dofs"
subdirectories.append(date_for_runs+"/"+"viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs048_p5_procs64")
filenames.append("turbulent_quantities.txt")
labels.append("$c_{+}$ NSFR")
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
    legend_fontSize_input=14)
#-----------------------------------------------------
#=====================================================

#=====================================================
# DOFs: 48^3
#-----------------------------------------------------
reinit_inputs()
date_for_runs="2022-11-11"
figure_subdirectory=date_for_runs
figure_title = "TGV at $Re_{\\infty}=1600$, $24^{3}$ DOFs ($P5$, $n_{el}=4^{3}$), CFL $=0.1$"
figure_title = " "
figure_filename_postfix = "dofs024_tgv"
legend_inside_input=True
#-----------------------------------------------------
date_for_runs="2022-11-09_24dofs"
subdirectories.append(date_for_runs+"/"+"viscous_TGV_ILES_std_strong_DG_lax_friedrichs_dofs024_p5_procs16")
filenames.append("turbulent_quantities.txt")
labels.append("Strong DG")
black_line_flag.append(False)
dashed_line_flag.append(False)
#-----------------------------------------------------
date_for_runs="2022-11-09_24dofs"
subdirectories.append(date_for_runs+"/"+"viscous_TGV_ILES_cDG_IR_two_point_flux_dofs024_p5_procs16")
filenames.append("turbulent_quantities.txt")
labels.append("$c_{DG}$ NSFR")
black_line_flag.append(False)
dashed_line_flag.append(False)
#-----------------------------------------------------
date_for_runs="2022-11-09_24dofs"
subdirectories.append(date_for_runs+"/"+"viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs024_p5_procs16")
filenames.append("turbulent_quantities.txt")
labels.append("$c_{+}$ NSFR")
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
    legend_fontSize_input=14)
#-----------------------------------------------------
#=====================================================

#=====================================================
# DOFs: 24^3, 48^3, 96^3 | c: cDG
#-----------------------------------------------------
reinit_inputs()
date_for_runs="2022-11-11"
figure_subdirectory=date_for_runs
figure_title = "TGV at $Re_{\\infty}=1600$, $P5$, CFL $=0.1$, $c_{DG}$ NSFR (IR)"
figure_title = " "
figure_filename_postfix = "all_dofs_cDG_tgv"
legend_inside_input=True
#-----------------------------------------------------
date_for_runs="2022-11-09_24dofs"
subdirectories.append(date_for_runs+"/"+"viscous_TGV_ILES_cDG_IR_two_point_flux_dofs024_p5_procs16")
filenames.append("turbulent_quantities.txt")
labels.append("$24^{3}$ DOFs")
black_line_flag.append(False)
dashed_line_flag.append(False)
#-----------------------------------------------------
date_for_runs="2022-11-09_48dofs"
subdirectories.append(date_for_runs+"/"+"viscous_TGV_ILES_cDG_IR_two_point_flux_dofs048_p5_procs64")
filenames.append("turbulent_quantities.txt")
labels.append("$48^{3}$ DOFs")
black_line_flag.append(False)
dashed_line_flag.append(False)
#-----------------------------------------------------
date_for_runs="2022-11-09_96dofs"
subdirectories.append(date_for_runs+"/"+"viscous_TGV_ILES_cDG_IR_two_point_flux_dofs096_p5_procs512")
filenames.append("turbulent_quantities.txt")
labels.append("$96^{3}$ DOFs")
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
    legend_fontSize_input=14)
#-----------------------------------------------------
#=====================================================

#=====================================================
# DOFs: 24^3, 48^3, 96^3 | c: cPlus
#-----------------------------------------------------
reinit_inputs()
date_for_runs="2022-11-11"
figure_subdirectory=date_for_runs
figure_title = "TGV at $Re_{\\infty}=1600$, $P5$, CFL $=0.1$, $c_{+}$ NSFR (IR)"
figure_title = " "
figure_filename_postfix = "all_dofs_cPlus_tgv"
legend_inside_input=True
#-----------------------------------------------------
date_for_runs="2022-11-09_24dofs"
subdirectories.append(date_for_runs+"/"+"viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs024_p5_procs16")
filenames.append("turbulent_quantities.txt")
labels.append("$24^{3}$ DOFs")
black_line_flag.append(False)
dashed_line_flag.append(False)
#-----------------------------------------------------
date_for_runs="2022-11-09_48dofs"
subdirectories.append(date_for_runs+"/"+"viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs048_p5_procs64")
filenames.append("turbulent_quantities.txt")
labels.append("$48^{3}$ DOFs")
black_line_flag.append(False)
dashed_line_flag.append(False)
#-----------------------------------------------------
date_for_runs="2022-11-09_96dofs"
subdirectories.append(date_for_runs+"/"+"viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs096_p5_procs512")
filenames.append("turbulent_quantities.txt")
labels.append("$96^{3}$ DOFs")
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
    legend_fontSize_input=14)
#-----------------------------------------------------
#=====================================================
# DOFs: 24^3, 48^3, 96^3 | c: cPlus
#-----------------------------------------------------
reinit_inputs()
date_for_runs="2022-11-11"
figure_subdirectory=date_for_runs
figure_title = "TGV at $Re_{\\infty}=1600$, $P5$, CFL $=0.1$, $c_{+}$ NSFR (IR)"
figure_title = " "
figure_filename_postfix = "all_dofs_cDG_and_cPlus_tgv"
legend_inside_input=True
#-----------------------------------------------------
date_for_runs="2022-11-09_24dofs"
subdirectories.append(date_for_runs+"/"+"viscous_TGV_ILES_cDG_IR_two_point_flux_dofs024_p5_procs16")
filenames.append("turbulent_quantities.txt")
labels.append("$24^{3}$, $c_{DG}$ NSFR")
black_line_flag.append(False)
dashed_line_flag.append(False)
#-----------------------------------------------------
date_for_runs="2022-11-09_24dofs"
subdirectories.append(date_for_runs+"/"+"viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs024_p5_procs16")
filenames.append("turbulent_quantities.txt")
labels.append("$24^{3}$, $c_{+}$ NSFR")
black_line_flag.append(False)
dashed_line_flag.append(False)
#-----------------------------------------------------
date_for_runs="2022-11-09_48dofs"
subdirectories.append(date_for_runs+"/"+"viscous_TGV_ILES_cDG_IR_two_point_flux_dofs048_p5_procs64")
filenames.append("turbulent_quantities.txt")
labels.append("$48^{3}$, $c_{DG}$ NSFR")
black_line_flag.append(False)
dashed_line_flag.append(False)
#-----------------------------------------------------
date_for_runs="2022-11-09_48dofs"
subdirectories.append(date_for_runs+"/"+"viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs048_p5_procs64")
filenames.append("turbulent_quantities.txt")
labels.append("$48^{3}$, $c_{+}$ NSFR")
black_line_flag.append(False)
dashed_line_flag.append(False)
#-----------------------------------------------------
date_for_runs="2022-11-09_96dofs"
subdirectories.append(date_for_runs+"/"+"viscous_TGV_ILES_cDG_IR_two_point_flux_dofs096_p5_procs512")
filenames.append("turbulent_quantities.txt")
labels.append("$96^{3}$, $c_{DG}$ NSFR")
black_line_flag.append(False)
dashed_line_flag.append(False)
#-----------------------------------------------------
date_for_runs="2022-11-09_96dofs"
subdirectories.append(date_for_runs+"/"+"viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs096_p5_procs512")
filenames.append("turbulent_quantities.txt")
labels.append("$96^{3}$, $c_{+}$ NSFR")
black_line_flag.append(False)
dashed_line_flag.append(False)

clr_input_store = ['tab:blue','tab:blue','tab:red','tab:red','tab:green','tab:green']#,'tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
mrkr_input_store = []
lnstl_input_store = ['solid','dashed','solid','dashed','solid','dashed']
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
    clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store,
    legend_fontSize_input=14)
#-----------------------------------------------------
#=====================================================