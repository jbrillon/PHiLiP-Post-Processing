#-----------------------------------------------------
# Import public libraries
import numpy as np # NumPy: contains basic numerical routines
#-----------------------------------------------------
import sys
# load tools
sys.path.append("/Users/Julien/PHiLiP-Post-Processing/src/tools");
from assemble_mpi_flow_field_files_and_reorder import assemble_mpi_flow_field_files_and_reorder
from generate_spectra_files import generate_spectra_file_from_flow_field_file
# load submodules
sys.path.append("/Users/Julien/PHiLiP-Post-Processing/submodules/quickplotlib/lib"); import quickplotlib as qp
#-----------------------------------------------------
#=====================================================
# Helper functions
#=====================================================
#-----------------------------------------------------
def append_to_plot(x_,y_,label_):
    global x,y,labels
    labels.append(label_);x.append(x_);y.append(y_)
#-----------------------------------------------------
def batch_append_to_plot(paths_,labels_,filename):
    global x,y,labels
    for i,path in enumerate(paths_):
        spectra = np.loadtxt(path+filename)
        append_to_plot(spectra[:,0],spectra[:,1],labels_[i])
#-----------------------------------------------------

#=====================================================
# Global variables
#=====================================================
global x,y,labels
x=[];y=[];labels=[];

# # =====================================================
# # 24 DOF check
# # =====================================================
# title_label = "TGV at $t=8.0$"
# figure_filename = "spectra_t8"
# # load files
# # 
# # filename_without_extension="/Users/Julien/julien_phd/post_processing/data/taylor_green_vortex/2022-11-09_24dofs/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs024_p5_procs16/flow_field_files/velocity_vorticity-0_reordered_equidistant"
# # filename_without_extension="/Users/Julien/julien_phd/post_processing/data/taylor_green_vortex/2022-11-09_48dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs048_p5_procs64/flow_field_files/velocity_vorticity-0_reordered_equidistant"
# # filename_without_extension="/Users/Julien/julien_phd/post_processing/data/taylor_green_vortex/2022-11-09_48dofs/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs048_p5_procs64/flow_field_files/velocity_vorticity-0_reordered_equidistant"
# # filename_without_extension="/Users/Julien/julien_phd/post_processing/data/taylor_green_vortex/2022-11-09_96dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs096_p5_procs512/flow_field_files/velocity_vorticity-0_reordered"
# # filename_without_extension="/Users/Julien/julien_phd/post_processing/data/taylor_green_vortex/2022-11-09_96dofs/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs096_p5_procs512/flow_field_files/velocity_vorticity-0_reordered"

# # generate_spectra_file_from_flow_field_file("/Users/Julien/DHIT-Flow-Setup/dofs048_p5_velocity/velocity_equidistant_nodes","fld",n_skiprows=0,use_TurboGenPy=True)
# # spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/dofs048_p5_velocity/velocity_equidistant_nodes_spectra.fld")
# # append_to_plot(spectra[:,0],spectra[:,1],"PHiLiP input")

# spectra = np.loadtxt("/Users/Julien/PHiLiP-Post-Processing/cases/taylor_green_vortex/data/mastellone2016_dns_spectra_t8.txt",skiprows=1,delimiter=',')
# append_to_plot(spectra[:,0],spectra[:,1],"Spectral t=8.0")

# spectra = np.loadtxt("/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_24dofs/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs024_p5_procs16/flow_field_files/velocity_vorticity-0_reordered_spectra.dat")
# append_to_plot(spectra[:,0],spectra[:,1],"$c_{DG}$ $24^3$")

# spectra = np.loadtxt("/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_24dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs024_p5_procs16/flow_field_files/velocity_vorticity-0_reordered_spectra.dat")
# append_to_plot(spectra[:,0],spectra[:,1],"$c_{+}$ $24^3$")

# spectra = np.loadtxt("/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs048_p5_procs64/flow_field_files/velocity_vorticity-0_reordered_spectra.dat")
# append_to_plot(spectra[:,0],spectra[:,1],"$c_{DG}$ $48^3$")

# spectra = np.loadtxt("/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs048_p5_procs64/flow_field_files/velocity_vorticity-0_reordered_spectra.dat")
# append_to_plot(spectra[:,0],spectra[:,1],"$c_{+}$ $48^3$")

# spectra = np.loadtxt("/Users/Julien/julien_phd/post_processing/data/taylor_green_vortex/2022-11-09_96dofs/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs096_p5_procs512/flow_field_files/velocity_vorticity-0_reordered_equidistant_spectra.dat")
# append_to_plot(spectra[:,0],spectra[:,1],"$c_{DG}$ $96^3$")

# spectra = np.loadtxt("/Users/Julien/julien_phd/post_processing/data/taylor_green_vortex/2022-11-09_96dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs096_p5_procs512/flow_field_files/velocity_vorticity-0_reordered_equidistant_spectra.dat")
# append_to_plot(spectra[:,0],spectra[:,1],"$c_{+}$ $96^3$")

# spectra = np.loadtxt("/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_ILES_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/flow_field_files/velocity_vorticity-0_reordered_spectra.dat")
# append_to_plot(spectra[:,0],spectra[:,1],"$c_{DG}$ + L$^{2}$Roe $96^3$")

# # spectra = np.loadtxt("/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_ILES_cHU_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/flow_field_files/velocity_vorticity-0_reordered_spectra.dat")
# # append_to_plot(spectra[:,0],spectra[:,1],"$c_{HU}$ + L$^{2}$Roe $96^3$")

# # spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/dofs048_p5_velocity/flow_field_files/velocity_vorticity-0_reordered_spectra.dat",skiprows=0,dtype=np.float64)
# # append_to_plot(spectra[:,0],spectra[:,1],"PHiLiP t=0.00")

# # spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/dofs048_p5_velocity/flow_field_files/velocity_vorticity-1_reordered_spectra.dat",skiprows=0,dtype=np.float64)
# # append_to_plot(spectra[:,0],spectra[:,1],"PHiLiP t=0.50")

# # spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/dofs048_p5_velocity/flow_field_files/velocity_vorticity-2_reordered_spectra.dat",skiprows=0,dtype=np.float64)
# # append_to_plot(spectra[:,0],spectra[:,1],"PHiLiP t=0.75")

# # generate_spectra_file_from_flow_field_file("/Users/Julien/DHIT-Flow-Setup/dofs048_p5_velocity/velocity_equidistant_nodes","fld",n_skiprows=0,use_TurboGenPy=True)
# # spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/dofs048_p5_velocity/velocity_equidistant_nodes_spectra.fld")
# # append_to_plot(spectra[:,0],spectra[:,1],"$48^3$ DOF")

# # generate_spectra_file_from_flow_field_file("/Users/Julien/DHIT-Flow-Setup/dofs128_p3_velocity/velocity_equidistant_nodes","fld",n_skiprows=0,use_TurboGenPy=True)
# # spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/dofs128_p3_velocity/velocity_equidistant_nodes_spectra.fld")
# # append_to_plot(spectra[:,0],spectra[:,1],"$128^3$ DOF")

# #=====================================================
# # Plotting function -- Spectra
# #=====================================================
# qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
#     title_label=title_label,
#     fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
#     nlegendcols=2,
#     xlimits=[2e0,1e2],ylimits=[1e-6,1e0],
#     markers=False,legend_on=True,legend_labels_tex=labels,
#     which_lines_black=[0]
#     # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
#     )
# qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
#     title_label=title_label,
#     fig_directory="figures",figure_filename=figure_filename+"_semilogy",log_axes="y",figure_filetype="pdf",
#     nlegendcols=2,
#     xlimits=[2e0,80],ylimits=[1e-6,1e0],
#     markers=False,legend_on=True,legend_labels_tex=labels,
#     which_lines_black=[0]
#     # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
#     )


# # t = 9 SPECTRA

# title_label = "TGV at $t=9.0$"
# figure_filename = "spectra_t9"
# x=[];y=[];labels=[];
# # spectra = np.loadtxt("/Users/Julien/PHiLiP-Post-Processing/cases/taylor_green_vortex/data/van_rees_dns_spectra.txt",skiprows=2,delimiter=',')
# # append_to_plot(spectra[:,0],spectra[:,1],"Spectral t=9")

# spectra = np.loadtxt("/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_24dofs/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs024_p5_procs16/flow_field_files/velocity_vorticity-1_reordered_spectra.dat")
# append_to_plot(spectra[:,0],spectra[:,1],"$c_{DG}$ $24^3$")

# spectra = np.loadtxt("/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_24dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs024_p5_procs16/flow_field_files/velocity_vorticity-1_reordered_spectra.dat")
# append_to_plot(spectra[:,0],spectra[:,1],"$c_{+}$ $24^3$")

# spectra = np.loadtxt("/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs048_p5_procs64/flow_field_files/velocity_vorticity-1_reordered_spectra.dat")
# append_to_plot(spectra[:,0],spectra[:,1],"$c_{DG}$ $48^3$")

# spectra = np.loadtxt("/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs048_p5_procs64/flow_field_files/velocity_vorticity-1_reordered_spectra.dat")
# append_to_plot(spectra[:,0],spectra[:,1],"$c_{+}$ $48^3$")

# spectra = np.loadtxt("/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_ILES_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/flow_field_files/velocity_vorticity-1_reordered_spectra.dat")
# append_to_plot(spectra[:,0],spectra[:,1],"$c_{DG}$ + L$^{2}$Roe $96^3$")

# # spectra = np.loadtxt("/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_ILES_cHU_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/flow_field_files/velocity_vorticity-1_reordered_spectra.dat")
# # append_to_plot(spectra[:,0],spectra[:,1],"$c_{HU}$ + L$^{2}$Roe $96^3$")

# qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
#     title_label=title_label,
#     fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
#     nlegendcols=2,
#     xlimits=[2e0,1e2],ylimits=[1e-6,1e0],
#     markers=False,legend_on=True,legend_labels_tex=labels,
#     which_lines_black=[0]
#     # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
#     )
# qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
#     title_label=title_label,
#     fig_directory="figures",figure_filename=figure_filename+"_semilogy",log_axes="y",figure_filetype="pdf",
#     nlegendcols=2,
#     xlimits=[2e0,80],ylimits=[1e-6,1e0],
#     markers=False,legend_on=True,legend_labels_tex=labels,
#     which_lines_black=[0]
#     # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
#     )

# reference batch paths for TGV results
# batch_paths=[ \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_24dofs/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs024_p5_procs16/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_24dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs024_p5_procs16/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_ILES_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_ILES_cHU_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_ILES_cSD_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_ILES_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_ILES_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cHU_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_ILES_cSD_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cHU_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cSD_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cHU_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cSD_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ]

# batch_labels = [ \
# "$24^3$ DOFs, $c_{DG}$ IR"\
# ,"$24^3$ DOFs, $c_{+}$ IR"\
# ,"$48^3$ DOFs, $c_{DG}$ IR"\
# ,"$48^3$ DOFs, $c_{+}$ IR"\
# ,"$48^3$ DOFs, $c_{DG}$ IR+L$^{2}$Roe"\
# ,"$48^3$ DOFs, $c_{+}$ IR+L$^{2}$Roe"\
# ,"$48^3$ DOFs, $c_{DG}$ IR+Smag.SGS"\
# ,"$48^3$ DOFs, $c_{+}$ IR+Smag.SGS"\
# ,"$48^3$ DOFs, $c_{DG}$ IR+L$^{2}$Roe+Smag.SGS"\
# ,"$48^3$ DOFs, $c_{+}$ IR+L$^{2}$Roe+Smag.SGS"\
# ,"$96^3$ DOFs, $c_{DG}$ IR"\
# ,"$96^3$ DOFs, $c_{+}$ IR"\
# ,"$96^3$ DOFs, $c_{HU}$ IR"\
# ,"$96^3$ DOFs, $c_{SD}$ IR"\
# ,"$96^3$ DOFs, $c_{DG}$ IR+L$^{2}$Roe"\
# ,"$96^3$ DOFs, $c_{+}$ IR+L$^{2}$Roe"\
# ,"$96^3$ DOFs, $c_{DG}$ IR+Smag.SGS"\
# ,"$96^3$ DOFs, $c_{+}$ IR+Smag.SGS"\
# ,"$96^3$ DOFs, $c_{DG}$ IR+L$^{2}$Roe+Smag.SGS"\
# ,"$96^3$ DOFs, $c_{+}$ IR+L$^{2}$Roe+Smag.SGS"\
# ,"$96^3$ DOFs, $c_{HU}$ IR+L$^{2}$Roe"\
# ,"$96^3$ DOFs, $c_{SD}$ IR+L$^{2}$Roe"\
# ,"$96^3$ DOFs, $c_{HU}$ IR+Smag.SGS"\
# ,"$96^3$ DOFs, $c_{SD}$ IR+Smag.SGS"\
# ,"$96^3$ DOFs, $c_{HU}$ IR+L$^{2}$Roe+Smag.SGS"\
# ,"$96^3$ DOFs, $c_{SD}$ IR+L$^{2}$Roe+Smag.SGS"\
# ]

title_label = "$48^3$ DOFs, TGV at $t=8.0$"
figure_filename = "spectra_t8_48"
x=[];y=[];labels=[];
# reference result
spectra = np.loadtxt("/Users/Julien/PHiLiP-Post-Processing/cases/taylor_green_vortex/data/mastellone2016_dns_spectra_t8.txt",skiprows=1,delimiter=',')
append_to_plot(spectra[:,0],spectra[:,1],"DNS [Mastellone]")
batch_labels = [ \
# "$24^3$ DOFs, $c_{DG}$ IR"\
# ,"$24^3$ DOFs, $c_{+}$ IR"\
"$c_{DG}$ NSFR"\
,"$c_{+}$ NSFR"\
,"$c_{DG}$ NSFR + L$^{2}$Roe"\
,"$c_{+}$ NSFR + L$^{2}$Roe"\
# ,"$48^3$ DOFs, $c_{DG}$ IR+Smag.SGS"\
# ,"$48^3$ DOFs, $c_{+}$ IR+Smag.SGS"\
# ,"$48^3$ DOFs, $c_{DG}$ IR+L$^{2}$Roe+Smag.SGS"\
# ,"$48^3$ DOFs, $c_{+}$ IR+L$^{2}$Roe+Smag.SGS"\
# ,"$96^3$ DOFs, $c_{DG}$ IR"\
# ,"$96^3$ DOFs, $c_{+}$ IR"\
# ,"$96^3$ DOFs, $c_{HU}$ IR"\
# ,"$96^3$ DOFs, $c_{SD}$ IR"\
# ,"$96^3$ DOFs, $c_{DG}$ IR+L$^{2}$Roe"\
# ,"$96^3$ DOFs, $c_{+}$ IR+L$^{2}$Roe"\
# ,"$96^3$ DOFs, $c_{DG}$ IR+Smag.SGS"\
# ,"$96^3$ DOFs, $c_{+}$ IR+Smag.SGS"\
# ,"$96^3$ DOFs, $c_{DG}$ IR+L$^{2}$Roe+Smag.SGS"\
# ,"$96^3$ DOFs, $c_{+}$ IR+L$^{2}$Roe+Smag.SGS"\
# ,"$96^3$ DOFs, $c_{HU}$ IR+L$^{2}$Roe"\
# ,"$96^3$ DOFs, $c_{SD}$ IR+L$^{2}$Roe"\
# ,"$96^3$ DOFs, $c_{HU}$ IR+Smag.SGS"\
# ,"$96^3$ DOFs, $c_{SD}$ IR+Smag.SGS"\
# ,"$96^3$ DOFs, $c_{HU}$ IR+L$^{2}$Roe+Smag.SGS"\
# ,"$96^3$ DOFs, $c_{SD}$ IR+L$^{2}$Roe+Smag.SGS"\
]
batch_paths=[ \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_24dofs/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs024_p5_procs16/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_24dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs024_p5_procs16/"\
"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs048_p5_procs64/"\
,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs048_p5_procs64/"\
,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_ILES_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64/"\
,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_ILES_cHU_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_ILES_cSD_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_ILES_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_ILES_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cHU_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_ILES_cSD_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cHU_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cSD_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cHU_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cSD_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
]
batch_append_to_plot(batch_paths, batch_labels, "flow_field_files/velocity_vorticity-0_reordered_spectra.dat")
qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
    # title_label=title_label,
    fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
    nlegendcols=1,
    xlimits=[2e0,1e2],ylimits=[1e-6,1e-1],
    markers=False,legend_on=True,legend_labels_tex=labels,
    which_lines_black=[0],
    which_lines_markers=[0],transparent_legend=True,legend_border_on=False,grid_lines_on=False
    # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
    )
title_label = "$48^3$ DOFs, TGV at $t=9.0$"
figure_filename = "spectra_t9_48"
x=[];y=[];labels=[];
spectra = np.loadtxt("/Users/Julien/PHiLiP-Post-Processing/cases/taylor_green_vortex/data/carton2014_dns_spectra_t9.txt",skiprows=1,delimiter=',')
append_to_plot(spectra[:,0],spectra[:,1],"DNS [Carton]")
batch_append_to_plot(batch_paths, batch_labels, "flow_field_files/velocity_vorticity-1_reordered_spectra.dat")
qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
    # title_label=title_label,
    fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
    nlegendcols=1,
    xlimits=[2e0,1e2],ylimits=[1e-6,1e-1],
    markers=False,legend_on=True,legend_labels_tex=labels,
    which_lines_black=[0],
    # which_lines_markers=[0],
    transparent_legend=True,legend_border_on=False,grid_lines_on=False
    # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
    )

title_label = "TGV at $t=8.0$"
figure_filename = "spectra_t8_all"
x=[];y=[];labels=[];
# reference result
spectra = np.loadtxt("/Users/Julien/PHiLiP-Post-Processing/cases/taylor_green_vortex/data/mastellone2016_dns_spectra_t8.txt",skiprows=1,delimiter=',')
append_to_plot(spectra[:,0],spectra[:,1],"DNS [Mastellone]")
clr_input_store = ['k','tab:blue','tab:blue','tab:red','tab:red','tab:green','tab:green']#,'tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
mrkr_input_store = ['None','None','None','None','None','None','None']
lnstl_input_store = ['solid','solid','dashed','solid','dashed','solid','dashed']
batch_labels = [ \
"$24^3$, $c_{DG}$ NSFR"\
,"$24^3$, $c_{+}$ NSFR"\
,"$48^3$, $c_{DG}$ NSFR"\
,"$48^3$, $c_{+}$ NSFR"\
# ,"$48^3$ DOFs, $c_{DG}$ IR+L$^{2}$Roe"\
# ,"$48^3$ DOFs, $c_{+}$ IR+L$^{2}$Roe"\
# ,"$48^3$ DOFs, $c_{DG}$ IR+Smag.SGS"\
# ,"$48^3$ DOFs, $c_{+}$ IR+Smag.SGS"\
# ,"$48^3$ DOFs, $c_{DG}$ IR+L$^{2}$Roe+Smag.SGS"\
# ,"$48^3$ DOFs, $c_{+}$ IR+L$^{2}$Roe+Smag.SGS"\
,"$96^3$, $c_{DG}$ NSFR"\
,"$96^3$, $c_{+}$ NSFR"\
# ,"$96^3$ DOFs, $c_{HU}$ IR"\
# ,"$96^3$ DOFs, $c_{SD}$ IR"\
# ,"$96^3$ DOFs, $c_{DG}$ IR+L$^{2}$Roe"\
# ,"$96^3$ DOFs, $c_{+}$ IR+L$^{2}$Roe"\
# ,"$96^3$ DOFs, $c_{DG}$ IR+Smag.SGS"\
# ,"$96^3$ DOFs, $c_{+}$ IR+Smag.SGS"\
# ,"$96^3$ DOFs, $c_{DG}$ IR+L$^{2}$Roe+Smag.SGS"\
# ,"$96^3$ DOFs, $c_{+}$ IR+L$^{2}$Roe+Smag.SGS"\
# ,"$96^3$ DOFs, $c_{HU}$ IR+L$^{2}$Roe"\
# ,"$96^3$ DOFs, $c_{SD}$ IR+L$^{2}$Roe"\
# ,"$96^3$ DOFs, $c_{HU}$ IR+Smag.SGS"\
# ,"$96^3$ DOFs, $c_{SD}$ IR+Smag.SGS"\
# ,"$96^3$ DOFs, $c_{HU}$ IR+L$^{2}$Roe+Smag.SGS"\
# ,"$96^3$ DOFs, $c_{SD}$ IR+L$^{2}$Roe+Smag.SGS"\
]
batch_paths=[ \
"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_24dofs/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs024_p5_procs16/"\
,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_24dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs024_p5_procs16/"\
,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs048_p5_procs64/"\
,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_ILES_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64/"\
,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs096_p5_procs512/"\
,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_ILES_cHU_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_ILES_cSD_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_ILES_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_ILES_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cHU_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_ILES_cSD_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cHU_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cSD_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cHU_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cSD_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
]
batch_append_to_plot(batch_paths, batch_labels, "flow_field_files/velocity_vorticity-0_reordered_spectra.dat")
qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
    # title_label=title_label,
    fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
    nlegendcols=1,
    xlimits=[2e0,1e2],ylimits=[1e-6,1e-1],
    markers=False,legend_on=True,legend_labels_tex=labels,
    which_lines_black=[0],
    which_lines_markers=[0],transparent_legend=True,legend_border_on=False,grid_lines_on=False,
    clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store
    # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
    )
title_label = "TGV at $t=9.0$"
figure_filename = "spectra_t9_all"
x=[];y=[];labels=[];
spectra = np.loadtxt("/Users/Julien/PHiLiP-Post-Processing/cases/taylor_green_vortex/data/carton2014_dns_spectra_t9.txt",skiprows=1,delimiter=',')
append_to_plot(spectra[:,0],spectra[:,1],"DNS [Carton]")
batch_append_to_plot(batch_paths, batch_labels, "flow_field_files/velocity_vorticity-1_reordered_spectra.dat")
qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
    # title_label=title_label,
    fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
    nlegendcols=1,
    xlimits=[2e0,1e2],ylimits=[1e-6,1e-1],
    markers=False,legend_on=True,legend_labels_tex=labels,
    which_lines_black=[0],
    # which_lines_markers=[0],
    transparent_legend=True,legend_border_on=False,grid_lines_on=False,
    clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store
    # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
    )
# title_label = "$96^3$ DOFs, TGV at $t=8.0$"
# figure_filename = "spectra_t8_96"
# x=[];y=[];labels=[];
# # reference result
# spectra = np.loadtxt("/Users/Julien/PHiLiP-Post-Processing/cases/taylor_green_vortex/data/mastellone2016_dns_spectra_t8.txt",skiprows=1,delimiter=',')
# append_to_plot(spectra[:,0],spectra[:,1],"DNS [Mastellone]")
# # - results
# batch_paths = [ \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs096_p5_procs512/", \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_dofs096_p5_procs512/", \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_ILES_cHU_IR_two_point_flux_dofs096_p5_procs512/", \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_ILES_cSD_IR_two_point_flux_dofs096_p5_procs512/", \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_ILES_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/", \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_ILES_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/", \
# # "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_dofs096_p5_procs512/", \
# # "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_dofs096_p5_procs512/", \
# # "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/", \
# # "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/", \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cHU_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/", \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_ILES_cSD_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/" \
# # "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cHU_IR_two_point_flux_dofs096_p5_procs512/", \
# # "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cSD_IR_two_point_flux_dofs096_p5_procs512/", \
# # "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cHU_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/", \
# # "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cSD_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/" \
# ]
# batch_labels = [ \
# "$c_{DG}$ IR", \
# "$c_{+}$ IR", \
# "$c_{HU}$ IR", \
# "$c_{SD}$ IR", \
# "$c_{DG}$ IR+L$^{2}$Roe", \
# "$c_{+}$ IR+L$^{2}$Roe", \
# # "$c_{DG}$ IR+Smag.SGS", \
# # "$c_{+}$ IR+Smag.SGS", \
# # "$c_{DG}$ IR+L$^{2}$Roe+Smag.SGS", \
# # "$c_{+}$ IR+L$^{2}$Roe+Smag.SGS", \
# "$c_{HU}$ IR+L$^{2}$Roe", \
# "$c_{SD}$ IR+L$^{2}$Roe" \
# # "$c_{HU}$ IR+Smag.SGS", \
# # "$c_{SD}$ IR+Smag.SGS", \
# # "$c_{HU}$ IR+L$^{2}$Roe+Smag.SGS", \
# # "$c_{SD}$ IR+L$^{2}$Roe+Smag.SGS", \
# ]
# batch_append_to_plot(batch_paths, batch_labels, "flow_field_files/velocity_vorticity-0_reordered_spectra.dat")
# qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
#     title_label=title_label,
#     fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
#     nlegendcols=1,
#     xlimits=[2e0,1e2],ylimits=[1e-6,1e-1],
#     markers=False,legend_on=True,legend_labels_tex=labels,
#     which_lines_black=[0],
#     which_lines_markers=[0],transparent_legend=True,legend_border_on=False,grid_lines_on=False
#     # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
#     )
# qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
#     title_label=title_label,
#     fig_directory="figures",figure_filename=figure_filename+"_zoom",log_axes="both",figure_filetype="pdf",
#     nlegendcols=1,
#     # xlimits=[2e1,8e1],ylimits=[1e-6,1e-3],
#     xlimits=[3e1,6e1],ylimits=[1e-5,2e-4],
#     markers=False,legend_on=True,legend_labels_tex=labels,
#     which_lines_black=[0],
#     which_lines_markers=[0],transparent_legend=True,legend_border_on=False,grid_lines_on=False
#     # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
#     )
# title_label = "$96^3$ DOFs, TGV at $t=9.0$"
# figure_filename = "spectra_t9_96"
# x=[];y=[];labels=[];
# # reference result
# # spectra = np.loadtxt("/Users/Julien/PHiLiP-Post-Processing/cases/taylor_green_vortex/data/mastellone2016_dns_spectra_t8.txt",skiprows=1,delimiter=',')
# # append_to_plot(spectra[:,0],spectra[:,1],"DNS [Mastellone] t=8.0")
# # - results
# # same as above
# batch_append_to_plot(batch_paths, batch_labels, "flow_field_files/velocity_vorticity-1_reordered_spectra.dat")
# qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
#     title_label=title_label,
#     fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
#     nlegendcols=1,
#     xlimits=[2e0,1e2],ylimits=[1e-6,1e-1],
#     markers=False,legend_on=True,legend_labels_tex=labels,
#     # which_lines_black=[0],
#     transparent_legend=True,legend_border_on=False,grid_lines_on=False
#     # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
#     )
# qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
#     title_label=title_label,
#     fig_directory="figures",figure_filename=figure_filename+"_zoom",log_axes="both",figure_filetype="pdf",
#     nlegendcols=1,
#     # xlimits=[2e1,8e1],ylimits=[1e-6,1e-3],
#     xlimits=[3e1,6e1],ylimits=[1e-5,2e-4],
#     markers=False,legend_on=True,legend_labels_tex=labels,
#     # which_lines_black=[0],
#     transparent_legend=True,legend_border_on=False,grid_lines_on=False
#     # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
#     )

# title_label = "$96^3$ DOFs, TGV at $t=8.0$"
# figure_filename = "spectra_t8_96_c"
# x=[];y=[];labels=[];
# # reference result
# spectra = np.loadtxt("/Users/Julien/PHiLiP-Post-Processing/cases/taylor_green_vortex/data/mastellone2016_dns_spectra_t8.txt",skiprows=1,delimiter=',')
# append_to_plot(spectra[:,0],spectra[:,1],"Spectral")
# # - results
# batch_paths = [ \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs096_p5_procs512/", \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_dofs096_p5_procs512/", \
# # "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_ILES_cHU_IR_two_point_flux_dofs096_p5_procs512/", \
# # "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_ILES_cSD_IR_two_point_flux_dofs096_p5_procs512/", \
# # "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_ILES_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/", \
# # "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_ILES_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/", \
# # "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_dofs096_p5_procs512/", \
# # "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_dofs096_p5_procs512/", \
# # "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/", \
# # "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/", \
# # "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cHU_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/", \
# # "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_ILES_cSD_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/" \
# # "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cHU_IR_two_point_flux_dofs096_p5_procs512/", \
# # "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cSD_IR_two_point_flux_dofs096_p5_procs512/", \
# # "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cHU_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/", \
# # "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cSD_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/" \
# ]
# batch_labels = [ \
# "$c_{DG}$ IR", \
# "$c_{+}$ IR", \
# # "$c_{HU}$ IR", \
# # "$c_{SD}$ IR", \
# # "$c_{DG}$ IR+L$^{2}$Roe", \
# # "$c_{+}$ IR+L$^{2}$Roe", \
# # "$c_{DG}$ IR+Smag.SGS", \
# # "$c_{+}$ IR+Smag.SGS", \
# # "$c_{DG}$ IR+L$^{2}$Roe+Smag.SGS", \
# # "$c_{+}$ IR+L$^{2}$Roe+Smag.SGS", \
# # "$c_{HU}$ IR+L$^{2}$Roe", \
# # "$c_{SD}$ IR+L$^{2}$Roe" \
# # "$c_{HU}$ IR+Smag.SGS", \
# # "$c_{SD}$ IR+Smag.SGS", \
# # "$c_{HU}$ IR+L$^{2}$Roe+Smag.SGS", \
# # "$c_{SD}$ IR+L$^{2}$Roe+Smag.SGS", \
# ]
# batch_append_to_plot(batch_paths, batch_labels, "flow_field_files/velocity_vorticity-0_reordered_spectra.dat")
# qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
#     title_label=title_label,
#     fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
#     nlegendcols=1,
#     xlimits=[2e0,1e2],ylimits=[1e-6,1e-1],
#     markers=False,legend_on=True,legend_labels_tex=labels,
#     which_lines_black=[0],
#     which_lines_markers=[0],transparent_legend=True,legend_border_on=False,grid_lines_on=False
#     # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
#     )
# qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
#     title_label=title_label,
#     fig_directory="figures",figure_filename=figure_filename+"_zoom",log_axes="both",figure_filetype="pdf",
#     nlegendcols=1,
#     # xlimits=[2e1,8e1],ylimits=[1e-6,1e-3],
#     xlimits=[3e1,6e1],ylimits=[1e-5,2e-4],
#     markers=False,legend_on=True,legend_labels_tex=labels,
#     which_lines_black=[0],
#     which_lines_markers=[0],transparent_legend=True,legend_border_on=False,grid_lines_on=False
#     # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
#     )
# title_label = "$96^3$ DOFs, TGV at $t=9.0$"
# figure_filename = "spectra_t9_96_c"
# x=[];y=[];labels=[];
# # reference result
# # spectra = np.loadtxt("/Users/Julien/PHiLiP-Post-Processing/cases/taylor_green_vortex/data/mastellone2016_dns_spectra_t8.txt",skiprows=1,delimiter=',')
# # append_to_plot(spectra[:,0],spectra[:,1],"Spectral t=8.0")
# # - results
# # same as above
# batch_append_to_plot(batch_paths, batch_labels, "flow_field_files/velocity_vorticity-1_reordered_spectra.dat")
# qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
#     title_label=title_label,
#     fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
#     nlegendcols=1,
#     xlimits=[2e0,1e2],ylimits=[1e-6,1e-1],
#     markers=False,legend_on=True,legend_labels_tex=labels,
#     # which_lines_black=[0],
#     transparent_legend=True,legend_border_on=False,grid_lines_on=False
#     # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
#     )
# qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
#     title_label=title_label,
#     fig_directory="figures",figure_filename=figure_filename+"_zoom",log_axes="both",figure_filetype="pdf",
#     nlegendcols=1,
#     # xlimits=[2e1,8e1],ylimits=[1e-6,1e-3],
#     xlimits=[3e1,6e1],ylimits=[1e-5,2e-4],
#     markers=False,legend_on=True,legend_labels_tex=labels,
#     # which_lines_black=[0],
#     transparent_legend=True,legend_border_on=False,grid_lines_on=False
#     # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
#     )


title_label = "$96^3$ DOFs, TGV at $t=8.0$"
figure_filename = "spectra_t8_96_riemann"
x=[];y=[];labels=[];
# reference result
spectra = np.loadtxt("/Users/Julien/PHiLiP-Post-Processing/cases/taylor_green_vortex/data/mastellone2016_dns_spectra_t8.txt",skiprows=1,delimiter=',')
append_to_plot(spectra[:,0],spectra[:,1],"DNS [Mastellone]")
# - results
batch_paths = [ \
"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs096_p5_procs512/", \
"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_dofs096_p5_procs512/", \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_ILES_cHU_IR_two_point_flux_dofs096_p5_procs512/", \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_ILES_cSD_IR_two_point_flux_dofs096_p5_procs512/", \
"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_ILES_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/", \
"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_ILES_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/", \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_dofs096_p5_procs512/", \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_dofs096_p5_procs512/", \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/", \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/", \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cHU_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/", \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_ILES_cSD_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/" \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cHU_IR_two_point_flux_dofs096_p5_procs512/", \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cSD_IR_two_point_flux_dofs096_p5_procs512/", \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cHU_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/", \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cSD_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/" \
]
batch_labels = [ \
"$c_{DG}$ NSFR", \
"$c_{+}$ NSFR", \
# "$c_{HU}$ IR", \
# "$c_{SD}$ IR", \
"$c_{DG}$ NSFR + L$^{2}$Roe", \
"$c_{+}$ NSFR + L$^{2}$Roe", \
# "$c_{DG}$ IR+Smag.SGS", \
# "$c_{+}$ IR+Smag.SGS", \
# "$c_{DG}$ IR+L$^{2}$Roe+Smag.SGS", \
# "$c_{+}$ IR+L$^{2}$Roe+Smag.SGS", \
# "$c_{HU}$ IR+L$^{2}$Roe", \
# "$c_{SD}$ IR+L$^{2}$Roe" \
# "$c_{HU}$ IR+Smag.SGS", \
# "$c_{SD}$ IR+Smag.SGS", \
# "$c_{HU}$ IR+L$^{2}$Roe+Smag.SGS", \
# "$c_{SD}$ IR+L$^{2}$Roe+Smag.SGS", \
]
batch_append_to_plot(batch_paths, batch_labels, "flow_field_files/velocity_vorticity-0_reordered_spectra.dat")
qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
    # title_label=title_label,
    fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
    nlegendcols=1,
    xlimits=[2e0,1e2],ylimits=[1e-6,1e-1],
    markers=False,legend_on=True,legend_labels_tex=labels,
    which_lines_black=[0],
    which_lines_markers=[0],transparent_legend=True,legend_border_on=False,grid_lines_on=False
    # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
    )
qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
    # title_label=title_label,
    fig_directory="figures",figure_filename=figure_filename+"_zoom",log_axes="both",figure_filetype="pdf",
    nlegendcols=1,
    # xlimits=[2e1,8e1],ylimits=[1e-6,1e-3],
    xlimits=[3e1,6e1],ylimits=[1e-5,2e-4],
    markers=False,legend_on=True,legend_labels_tex=labels,
    which_lines_black=[0],
    which_lines_markers=[0],transparent_legend=True,legend_border_on=False,grid_lines_on=False
    # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
    )
title_label = "$96^3$ DOFs, TGV at $t=9.0$"
figure_filename = "spectra_t9_96_riemann"
x=[];y=[];labels=[];
# reference result
spectra = np.loadtxt("/Users/Julien/PHiLiP-Post-Processing/cases/taylor_green_vortex/data/carton2014_dns_spectra_t9.txt",skiprows=1,delimiter=',')
append_to_plot(spectra[:,0],spectra[:,1],"DNS [Carton]")
# - results
# same as above
batch_append_to_plot(batch_paths, batch_labels, "flow_field_files/velocity_vorticity-1_reordered_spectra.dat")
qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
    # title_label=title_label,
    fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
    nlegendcols=1,
    xlimits=[2e0,1e2],ylimits=[1e-6,1e-1],
    markers=False,legend_on=True,legend_labels_tex=labels,
    # which_lines_black=[0],
    transparent_legend=True,legend_border_on=False,grid_lines_on=False
    # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
    )
qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
    # title_label=title_label,
    fig_directory="figures",figure_filename=figure_filename+"_zoom",log_axes="both",figure_filetype="pdf",
    nlegendcols=1,
    xlimits=[3e1,6e1],ylimits=[1e-5,2e-4],
    markers=False,legend_on=True,legend_labels_tex=labels,
    which_lines_black=[0],
    transparent_legend=True,legend_border_on=False,grid_lines_on=False
    # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
    )

title_label = "$48^3$ DOFs, TGV at $t=8.0$"
figure_filename = "spectra_t8_48_cPlus_sgs_filters"
x=[];y=[];labels=[];
# reference result
spectra = np.loadtxt("/Users/Julien/PHiLiP-Post-Processing/cases/taylor_green_vortex/data/mastellone2016_dns_spectra_t8.txt",skiprows=1,delimiter=',')
append_to_plot(spectra[:,0],spectra[:,1],"DNS [Mastellone]")
# - results
batch_paths = [ \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs096_p5_procs512/", \
"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs048_p5_procs64/",\
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_ILES_cHU_IR_two_point_flux_dofs096_p5_procs512/", \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_ILES_cSD_IR_two_point_flux_dofs096_p5_procs512/", \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_ILES_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/", \
"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64/", \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_dofs096_p5_procs512/", \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_dofs096_p5_procs512/", \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/", \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/", \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cHU_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/", \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_ILES_cSD_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/" \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cHU_IR_two_point_flux_dofs096_p5_procs512/", \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cSD_IR_two_point_flux_dofs096_p5_procs512/", \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cHU_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/", \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cSD_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/" \
"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_dofs048_p5_procs64_filter216timeslarger/", \
"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_dofs048_p5_procs64_filter36timeslarger/", \
"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_dofs048_p5_procs64_filter6timeslarger/" \
]
batch_labels = [ \
# "$c_{DG}$ IR", \
"$c_{+}$ NSFR", \
# "$c_{HU}$ IR", \
# "$c_{SD}$ IR", \
# "$c_{DG}$ IR+L$^{2}$Roe", \
"$c_{+}$ NSFR + L$^{2}$Roe", \
# "$c_{DG}$ IR+Smag.SGS", \
# "$c_{+}$ IR+Smag.SGS", \
# "$c_{DG}$ IR+L$^{2}$Roe+Smag.SGS", \
# "$c_{+}$ IR+L$^{2}$Roe+Smag.SGS", \
# "$c_{HU}$ IR+L$^{2}$Roe", \
# "$c_{SD}$ IR+L$^{2}$Roe" \
# "$c_{HU}$ IR+Smag.SGS", \
# "$c_{SD}$ IR+Smag.SGS", \
# "$c_{HU}$ IR+L$^{2}$Roe+Smag.SGS", \
# "$c_{SD}$ IR+L$^{2}$Roe+Smag.SGS", \
"$c_{+}$ NSFR + SGS $\\Delta=V_{cell}$", \
"$c_{+}$ NSFR + SGS $\\Delta=\\frac{V_{cell}}{p+1}$", \
"$c_{+}$ NSFR + SGS $\\Delta=\\frac{V_{cell}}{(p+1)^{2}}$", \
]
batch_append_to_plot(batch_paths, batch_labels, "flow_field_files/velocity_vorticity-0_reordered_spectra.dat")
qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
    # title_label=title_label,
    fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
    nlegendcols=1,
    xlimits=[2e0,1e2],ylimits=[1e-6,1e-1],
    markers=False,legend_on=True,legend_labels_tex=labels,
    which_lines_black=[0],
    which_lines_markers=[0],transparent_legend=True,legend_border_on=False,grid_lines_on=False
    # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
    )
# qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
#     title_label=title_label,
#     fig_directory="figures",figure_filename=figure_filename+"_zoom",log_axes="both",figure_filetype="pdf",
#     nlegendcols=1,
#     # xlimits=[2e1,8e1],ylimits=[1e-6,1e-3],
#     xlimits=[3e1,6e1],ylimits=[1e-5,2e-4],
#     markers=False,legend_on=True,legend_labels_tex=labels,
#     which_lines_black=[0],
#     which_lines_markers=[0],transparent_legend=True,legend_border_on=False,grid_lines_on=False
#     # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
#     )
title_label = "$48^3$ DOFs, TGV at $t=9.0$"
figure_filename = "spectra_t9_48_cPlus_sgs_filters"
x=[];y=[];labels=[];
# reference result
spectra = np.loadtxt("/Users/Julien/PHiLiP-Post-Processing/cases/taylor_green_vortex/data/carton2014_dns_spectra_t9.txt",skiprows=1,delimiter=',')
append_to_plot(spectra[:,0],spectra[:,1],"DNS [Carton]")
# - results
# same as above
batch_append_to_plot(batch_paths, batch_labels, "flow_field_files/velocity_vorticity-1_reordered_spectra.dat")
qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
    # title_label=title_label,
    fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
    nlegendcols=1,
    xlimits=[2e0,1e2],ylimits=[1e-6,1e-1],
    markers=False,legend_on=True,legend_labels_tex=labels,
    which_lines_black=[0],
    transparent_legend=True,legend_border_on=False,grid_lines_on=False
    # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
    )
# qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
#     title_label=title_label,
#     fig_directory="figures",figure_filename=figure_filename+"_zoom",log_axes="both",figure_filetype="pdf",
#     nlegendcols=1,
#     xlimits=[3e1,6e1],ylimits=[1e-5,2e-4],
#     markers=False,legend_on=True,legend_labels_tex=labels,
#     # which_lines_black=[0],
#     transparent_legend=True,legend_border_on=False,grid_lines_on=False
#     # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
#     )

title_label = "$48^3$ DOFs, TGV at $t=8.0$"
figure_filename = "spectra_t8_48_cPlus"
x=[];y=[];labels=[];
# reference result
spectra = np.loadtxt("/Users/Julien/PHiLiP-Post-Processing/cases/taylor_green_vortex/data/mastellone2016_dns_spectra_t8.txt",skiprows=1,delimiter=',')
append_to_plot(spectra[:,0],spectra[:,1],"DNS [Mastellone]")
# - results
batch_paths = [ \
"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs048_p5_procs64/",\
"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64/", \
"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_dofs048_p5_procs64_filter36timeslarger/", \
"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64_filter36timeslarger/"
]
batch_labels = [ \
"$c_{+}$ NSFR", \
"$c_{+}$ NSFR-L2R", \
"$c_{+}$ NSFR-SGS", \
"$c_{+}$ NSFR-L2R-SGS", \
]
batch_append_to_plot(batch_paths, batch_labels, "flow_field_files/velocity_vorticity-0_reordered_spectra.dat")
qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
    # title_label=title_label,
    fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
    nlegendcols=1,
    xlimits=[2e0,1e2],ylimits=[1e-6,1e-1],
    markers=False,legend_on=True,legend_labels_tex=labels,
    which_lines_black=[0],
    # which_lines_markers=[0],
    transparent_legend=True,legend_border_on=False,grid_lines_on=False
    # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
    )
# qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
#     title_label=title_label,
#     fig_directory="figures",figure_filename=figure_filename+"_zoom",log_axes="both",figure_filetype="pdf",
#     nlegendcols=1,
#     # xlimits=[2e1,8e1],ylimits=[1e-6,1e-3],
#     xlimits=[3e1,6e1],ylimits=[1e-5,2e-4],
#     markers=False,legend_on=True,legend_labels_tex=labels,
#     which_lines_black=[0],
#     which_lines_markers=[0],transparent_legend=True,legend_border_on=False,grid_lines_on=False
#     # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
#     )
title_label = "$48^3$ DOFs, TGV at $t=9.0$"
figure_filename = "spectra_t9_48_cPlus"
x=[];y=[];labels=[];
# reference result
spectra = np.loadtxt("/Users/Julien/PHiLiP-Post-Processing/cases/taylor_green_vortex/data/carton2014_dns_spectra_t9.txt",skiprows=1,delimiter=',')
append_to_plot(spectra[:,0],spectra[:,1],"DNS [Carton]")
# - results
# same as above
batch_append_to_plot(batch_paths, batch_labels, "flow_field_files/velocity_vorticity-1_reordered_spectra.dat")
qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
    # title_label=title_label,
    fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
    nlegendcols=1,
    xlimits=[2e0,1e2],ylimits=[1e-6,1e-1],
    markers=False,legend_on=True,legend_labels_tex=labels,
    which_lines_black=[0],
    transparent_legend=True,legend_border_on=False,grid_lines_on=False
    # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
    )
# qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
#     title_label=title_label,
#     fig_directory="figures",figure_filename=figure_filename+"_zoom",log_axes="both",figure_filetype="pdf",
#     nlegendcols=1,
#     xlimits=[3e1,6e1],ylimits=[1e-5,2e-4],
#     markers=False,legend_on=True,legend_labels_tex=labels,
#     # which_lines_black=[0],
#     transparent_legend=True,legend_border_on=False,grid_lines_on=False
#     # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
#     )

title_label = "$48^3$ DOFs, TGV at $t=8.0$"
figure_filename = "spectra_t8_48_cDG"
x=[];y=[];labels=[];
# reference result
spectra = np.loadtxt("/Users/Julien/PHiLiP-Post-Processing/cases/taylor_green_vortex/data/mastellone2016_dns_spectra_t8.txt",skiprows=1,delimiter=',')
append_to_plot(spectra[:,0],spectra[:,1],"DNS [Mastellone]")
# - results
batch_paths = [ \
"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs048_p5_procs64/",\
"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_ILES_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64/", \
"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_dofs048_p5_procs64_filter36timeslarger/", \
"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64_filter36timeslarger/"
]
batch_labels = [ \
"$c_{DG}$ NSFR", \
"$c_{DG}$ NSFR-L2R", \
"$c_{DG}$ NSFR-SGS", \
"$c_{DG}$ NSFR-L2R-SGS", \
]
batch_append_to_plot(batch_paths, batch_labels, "flow_field_files/velocity_vorticity-0_reordered_spectra.dat")
qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
    # title_label=title_label,
    fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
    nlegendcols=1,
    xlimits=[2e0,1e2],ylimits=[1e-6,1e-1],
    markers=False,legend_on=True,legend_labels_tex=labels,
    which_lines_black=[0],
    # which_lines_markers=[0],
    transparent_legend=True,legend_border_on=False,grid_lines_on=False
    # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
    )
# qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
#     title_label=title_label,
#     fig_directory="figures",figure_filename=figure_filename+"_zoom",log_axes="both",figure_filetype="pdf",
#     nlegendcols=1,
#     # xlimits=[2e1,8e1],ylimits=[1e-6,1e-3],
#     xlimits=[3e1,6e1],ylimits=[1e-5,2e-4],
#     markers=False,legend_on=True,legend_labels_tex=labels,
#     which_lines_black=[0],
#     which_lines_markers=[0],transparent_legend=True,legend_border_on=False,grid_lines_on=False
#     # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
#     )
title_label = "$48^3$ DOFs, TGV at $t=9.0$"
figure_filename = "spectra_t9_48_cDG"
x=[];y=[];labels=[];
# reference result
spectra = np.loadtxt("/Users/Julien/PHiLiP-Post-Processing/cases/taylor_green_vortex/data/carton2014_dns_spectra_t9.txt",skiprows=1,delimiter=',')
append_to_plot(spectra[:,0],spectra[:,1],"DNS [Carton]")
# - results
# same as above
batch_append_to_plot(batch_paths, batch_labels, "flow_field_files/velocity_vorticity-1_reordered_spectra.dat")
qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
    # title_label=title_label,
    fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
    nlegendcols=1,
    xlimits=[2e0,1e2],ylimits=[1e-6,1e-1],
    markers=False,legend_on=True,legend_labels_tex=labels,
    which_lines_black=[0],
    transparent_legend=True,legend_border_on=False,grid_lines_on=False
    # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
    )
# qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
#     title_label=title_label,
#     fig_directory="figures",figure_filename=figure_filename+"_zoom",log_axes="both",figure_filetype="pdf",
#     nlegendcols=1,
#     xlimits=[3e1,6e1],ylimits=[1e-5,2e-4],
#     markers=False,legend_on=True,legend_labels_tex=labels,
#     # which_lines_black=[0],
#     transparent_legend=True,legend_border_on=False,grid_lines_on=False
#     # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
#     )

title_label = "$48^3$ DOFs, TGV at $t=8.0$"
figure_filename = "spectra_t8_48_cPlus_vs_cDG_L2R_SGS"
x=[];y=[];labels=[];
clr_input_store = ['k','tab:blue','tab:blue','tab:red','tab:red','tab:green','tab:green']#,'tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
mrkr_input_store = ['None','None','None','None','None','None','None']
lnstl_input_store = ['solid','dashed','solid','dashed','solid','dashed','solid']
# reference result
spectra = np.loadtxt("/Users/Julien/PHiLiP-Post-Processing/cases/taylor_green_vortex/data/mastellone2016_dns_spectra_t8.txt",skiprows=1,delimiter=',')
append_to_plot(spectra[:,0],spectra[:,1],"DNS [Mastellone]")
# - results
batch_paths = [ \
"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs048_p5_procs64/",\
"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs048_p5_procs64/",\
"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64_filter36timeslarger/",\
"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64_filter36timeslarger/"
]
batch_labels = [ \
"$c_{DG}$ NSFR", \
"$c_{+}$ NSFR", \
"$c_{DG}$ NSFR-L2R-SGS", \
"$c_{+}$ NSFR-L2R-SGS"
]
batch_append_to_plot(batch_paths, batch_labels, "flow_field_files/velocity_vorticity-0_reordered_spectra.dat")
qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
    # title_label=title_label,
    fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
    nlegendcols=1,
    xlimits=[2e0,1e2],ylimits=[1e-6,1e-1],
    markers=False,legend_on=True,legend_labels_tex=labels,
    which_lines_black=[0],
    # which_lines_markers=[0],
    transparent_legend=True,legend_border_on=False,grid_lines_on=False,
    clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store
    # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
    )
# qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
#     title_label=title_label,
#     fig_directory="figures",figure_filename=figure_filename+"_zoom",log_axes="both",figure_filetype="pdf",
#     nlegendcols=1,
#     # xlimits=[2e1,8e1],ylimits=[1e-6,1e-3],
#     xlimits=[3e1,6e1],ylimits=[1e-5,2e-4],
#     markers=False,legend_on=True,legend_labels_tex=labels,
#     which_lines_black=[0],
#     which_lines_markers=[0],transparent_legend=True,legend_border_on=False,grid_lines_on=False
#     # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
#     )
title_label = "$48^3$ DOFs, TGV at $t=9.0$"
figure_filename = "spectra_t9_48_cPlus_vs_cDG_L2R_SGS"
x=[];y=[];labels=[];
# reference result
spectra = np.loadtxt("/Users/Julien/PHiLiP-Post-Processing/cases/taylor_green_vortex/data/carton2014_dns_spectra_t9.txt",skiprows=1,delimiter=',')
append_to_plot(spectra[:,0],spectra[:,1],"DNS [Carton]")
# - results
# same as above
batch_append_to_plot(batch_paths, batch_labels, "flow_field_files/velocity_vorticity-1_reordered_spectra.dat")
qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
    # title_label=title_label,
    fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
    nlegendcols=1,
    xlimits=[2e0,1e2],ylimits=[1e-6,1e-1],
    markers=False,legend_on=True,legend_labels_tex=labels,
    which_lines_black=[0],
    transparent_legend=True,legend_border_on=False,grid_lines_on=False,
    clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store
    # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
    )

title_label = "$24^3$ DOFs, TGV at $t=8.0$"
figure_filename = "spectra_t8_24_cDG"
x=[];y=[];labels=[];
# reference result
spectra = np.loadtxt("/Users/Julien/PHiLiP-Post-Processing/cases/taylor_green_vortex/data/mastellone2016_dns_spectra_t8.txt",skiprows=1,delimiter=',')
append_to_plot(spectra[:,0],spectra[:,1],"DNS [Mastellone]")
# - results
batch_paths = [ \
"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_24dofs/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs024_p5_procs16/",\
"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_24dofs/viscous_TGV_ILES_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs024_p5_procs16/", \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_24dofs/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_dofs024_p5_procs16_filter36timeslarger/", \
"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_24dofs/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs024_p5_procs16_filter36timeslarger/"
]
batch_labels = [ \
"$c_{DG}$ NSFR", \
"$c_{DG}$ NSFR-L2R", \
# "$c_{DG}$ NSFR + SGS", \
"$c_{DG}$ NSFR-L2R-SGS", \
]
batch_append_to_plot(batch_paths, batch_labels, "flow_field_files/velocity_vorticity-0_reordered_spectra.dat")
qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
    # title_label=title_label,
    fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
    nlegendcols=1,
    xlimits=[2e0,1e2],ylimits=[1e-6,1e-1],
    markers=False,legend_on=True,legend_labels_tex=labels,
    which_lines_black=[0],
    which_lines_markers=[0],transparent_legend=True,legend_border_on=False,grid_lines_on=False
    # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
    )
# qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
#     title_label=title_label,
#     fig_directory="figures",figure_filename=figure_filename+"_zoom",log_axes="both",figure_filetype="pdf",
#     nlegendcols=1,
#     # xlimits=[2e1,8e1],ylimits=[1e-6,1e-3],
#     xlimits=[3e1,6e1],ylimits=[1e-5,2e-4],
#     markers=False,legend_on=True,legend_labels_tex=labels,
#     which_lines_black=[0],
#     which_lines_markers=[0],transparent_legend=True,legend_border_on=False,grid_lines_on=False
#     # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
#     )
title_label = "$24^3$ DOFs, TGV at $t=9.0$"
figure_filename = "spectra_t9_24_cDG"
x=[];y=[];labels=[];
# reference result
spectra = np.loadtxt("/Users/Julien/PHiLiP-Post-Processing/cases/taylor_green_vortex/data/carton2014_dns_spectra_t9.txt",skiprows=1,delimiter=',')
append_to_plot(spectra[:,0],spectra[:,1],"DNS [Carton]")
# - results
# same as above
batch_append_to_plot(batch_paths, batch_labels, "flow_field_files/velocity_vorticity-1_reordered_spectra.dat")
qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
    # title_label=title_label,
    fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
    nlegendcols=1,
    xlimits=[2e0,1e2],ylimits=[1e-6,1e-1],
    markers=False,legend_on=True,legend_labels_tex=labels,
    which_lines_black=[0],
    transparent_legend=True,legend_border_on=False,grid_lines_on=False
    # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
    )
# qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
#     title_label=title_label,
#     fig_directory="figures",figure_filename=figure_filename+"_zoom",log_axes="both",figure_filetype="pdf",
#     nlegendcols=1,
#     xlimits=[3e1,6e1],ylimits=[1e-5,2e-4],
#     markers=False,legend_on=True,legend_labels_tex=labels,
#     # which_lines_black=[0],
#     transparent_legend=True,legend_border_on=False,grid_lines_on=False
#     # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
#     )

title_label = "$24^3$ DOFs, TGV at $t=8.0$"
figure_filename = "spectra_t8_24_cPlus"
x=[];y=[];labels=[];
# reference result
spectra = np.loadtxt("/Users/Julien/PHiLiP-Post-Processing/cases/taylor_green_vortex/data/mastellone2016_dns_spectra_t8.txt",skiprows=1,delimiter=',')
append_to_plot(spectra[:,0],spectra[:,1],"DNS [Mastellone]")
# - results
batch_paths = [ \
"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_24dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs024_p5_procs16/",\
"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_24dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs024_p5_procs16/", \
# "/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_24dofs/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_dofs024_p5_procs16_filter36timeslarger/", \
"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_24dofs/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs024_p5_procs16_filter36timeslarger/"
]
batch_labels = [ \
"$c_{+}$ NSFR", \
"$c_{+}$ NSFR-L2R", \
# "$c_{+}$ NSFR + SGS", \
"$c_{+}$ NSFR-L2R-SGS", \
]
batch_append_to_plot(batch_paths, batch_labels, "flow_field_files/velocity_vorticity-0_reordered_spectra.dat")
qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
    # title_label=title_label,
    fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
    nlegendcols=1,
    xlimits=[2e0,1e2],ylimits=[1e-6,1e-1],
    markers=False,legend_on=True,legend_labels_tex=labels,
    which_lines_black=[0],
    which_lines_markers=[0],transparent_legend=True,legend_border_on=False,grid_lines_on=False
    # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
    )
# qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
#     title_label=title_label,
#     fig_directory="figures",figure_filename=figure_filename+"_zoom",log_axes="both",figure_filetype="pdf",
#     nlegendcols=1,
#     # xlimits=[2e1,8e1],ylimits=[1e-6,1e-3],
#     xlimits=[3e1,6e1],ylimits=[1e-5,2e-4],
#     markers=False,legend_on=True,legend_labels_tex=labels,
#     which_lines_black=[0],
#     which_lines_markers=[0],transparent_legend=True,legend_border_on=False,grid_lines_on=False
#     # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
#     )
title_label = "$24^3$ DOFs, TGV at $t=9.0$"
figure_filename = "spectra_t9_24_cPlus"
x=[];y=[];labels=[];
# reference result
spectra = np.loadtxt("/Users/Julien/PHiLiP-Post-Processing/cases/taylor_green_vortex/data/carton2014_dns_spectra_t9.txt",skiprows=1,delimiter=',')
append_to_plot(spectra[:,0],spectra[:,1],"DNS [Carton]")
# - results
# same as above
batch_append_to_plot(batch_paths, batch_labels, "flow_field_files/velocity_vorticity-1_reordered_spectra.dat")
qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
    # title_label=title_label,
    fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
    nlegendcols=1,
    xlimits=[2e0,1e2],ylimits=[1e-6,1e-1],
    markers=False,legend_on=True,legend_labels_tex=labels,
    which_lines_black=[0],
    transparent_legend=True,legend_border_on=False,grid_lines_on=False
    # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
    )
# qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
#     title_label=title_label,
#     fig_directory="figures",figure_filename=figure_filename+"_zoom",log_axes="both",figure_filetype="pdf",
#     nlegendcols=1,
#     xlimits=[3e1,6e1],ylimits=[1e-5,2e-4],
#     markers=False,legend_on=True,legend_labels_tex=labels,
#     # which_lines_black=[0],
#     transparent_legend=True,legend_border_on=False,grid_lines_on=False
#     # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
#     )

batch_paths=[ \
"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_24dofs/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs024_p5_procs16/"\
,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_24dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs024_p5_procs16/"\
,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs048_p5_procs64/"\
,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_ILES_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64/"\
,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs096_p5_procs512/"\
,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_ILES_cHU_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_ILES_cSD_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_ILES_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_ILES_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cHU_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_ILES_cSD_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cHU_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cSD_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cHU_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cSD_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
]

title_label = "TGV at $t=8.0$"
figure_filename = "spectra_t8_all_cDG_l2roe"
x=[];y=[];labels=[];
# reference result
spectra = np.loadtxt("/Users/Julien/PHiLiP-Post-Processing/cases/taylor_green_vortex/data/mastellone2016_dns_spectra_t8.txt",skiprows=1,delimiter=',')
append_to_plot(spectra[:,0],spectra[:,1],"DNS [Mastellone]")
clr_input_store = ['k','tab:blue','tab:blue','tab:red','tab:red','tab:green','tab:green']#,'tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
mrkr_input_store = ['None','None','None','None','None','None','None']
lnstl_input_store = ['solid','dashed','solid','dashed','solid','dashed','solid']
batch_labels = [ \
"$24^3$, $c_{DG}$ NSFR"\
,"$24^3$, $c_{DG}$ NSFR-L2R"\
,"$48^3$, $c_{DG}$ NSFR"\
,"$48^3$, $c_{DG}$ NSFR-L2R"\
,"$96^3$, $c_{DG}$ NSFR"\
,"$96^3$, $c_{DG}$ NSFR-L2R"\
]
batch_paths=[ \
"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_24dofs/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs024_p5_procs16/"\
,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_24dofs/viscous_TGV_ILES_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs024_p5_procs16/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_24dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs024_p5_procs16/"\
,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs048_p5_procs64/"\
,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_ILES_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64/"\
,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_ILES_cHU_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_ILES_cSD_IR_two_point_flux_dofs096_p5_procs512/"\
,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_ILES_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_ILES_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cHU_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_ILES_cSD_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cHU_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cSD_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cHU_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cSD_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
]
batch_append_to_plot(batch_paths, batch_labels, "flow_field_files/velocity_vorticity-0_reordered_spectra.dat")
qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
    # title_label=title_label,
    fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
    nlegendcols=1,
    xlimits=[2e0,1e2],ylimits=[1e-6,1e-1],
    markers=False,legend_on=True,legend_labels_tex=labels,
    which_lines_black=[0],
    # which_lines_markers=[0],
    transparent_legend=True,legend_border_on=False,grid_lines_on=False,
    clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store
    # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
    )
title_label = "TGV at $t=9.0$"
figure_filename = "spectra_t9_all_cDG_l2roe"
x=[];y=[];labels=[];
spectra = np.loadtxt("/Users/Julien/PHiLiP-Post-Processing/cases/taylor_green_vortex/data/carton2014_dns_spectra_t9.txt",skiprows=1,delimiter=',')
append_to_plot(spectra[:,0],spectra[:,1],"DNS [Carton]")
batch_append_to_plot(batch_paths, batch_labels, "flow_field_files/velocity_vorticity-1_reordered_spectra.dat")
qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
    # title_label=title_label,
    fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
    nlegendcols=1,
    xlimits=[2e0,1e2],ylimits=[1e-6,1e-1],
    markers=False,legend_on=True,legend_labels_tex=labels,
    which_lines_black=[0],
    # which_lines_markers=[0],
    transparent_legend=True,legend_border_on=False,grid_lines_on=False,
    clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store
    # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
    )

title_label = "TGV at $t=8.0$"
figure_filename = "spectra_t8_all_cPlus_l2roe"
x=[];y=[];labels=[];
# reference result
spectra = np.loadtxt("/Users/Julien/PHiLiP-Post-Processing/cases/taylor_green_vortex/data/mastellone2016_dns_spectra_t8.txt",skiprows=1,delimiter=',')
append_to_plot(spectra[:,0],spectra[:,1],"DNS [Mastellone]")
clr_input_store = ['k','tab:blue','tab:blue','tab:red','tab:red','tab:green','tab:green']#,'tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
mrkr_input_store = ['None','None','None','None','None','None','None']
lnstl_input_store = ['solid','dashed','solid','dashed','solid','dashed','solid']
batch_labels = [ \
"$24^3$, $c_{+}$ NSFR"\
,"$24^3$, $c_{+}$ NSFR-L2R"\
,"$48^3$, $c_{+}$ NSFR"\
,"$48^3$, $c_{+}$ NSFR-L2R"\
,"$96^3$, $c_{+}$ NSFR"\
,"$96^3$, $c_{+}$ NSFR-L2R"\
]
batch_paths=[ \
"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_24dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs024_p5_procs16/"\
,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_24dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs024_p5_procs16/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_24dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs024_p5_procs16/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs048_p5_procs64/"\
,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_ILES_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64/"\
,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_ILES_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_48dofs/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_ILES_cDG_IR_two_point_flux_dofs096_p5_procs512/"\
,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_ILES_cHU_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_ILES_cSD_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_ILES_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_ILES_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cDG_cPlus/viscous_TGV_LES_smagorinsky_cPlus_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cHU_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_ILES_cSD_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cHU_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cSD_IR_two_point_flux_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cHU_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
# ,"/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_TGV_SPECTRA_96dofs_cSD_cHU/viscous_TGV_LES_smagorinsky_cSD_IR_two_point_flux_with_l2roe_dissipation_dofs096_p5_procs512/"\
]
batch_append_to_plot(batch_paths, batch_labels, "flow_field_files/velocity_vorticity-0_reordered_spectra.dat")
qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
    # title_label=title_label,
    fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
    nlegendcols=1,
    xlimits=[2e0,1e2],ylimits=[1e-6,1e-1],
    markers=False,legend_on=True,legend_labels_tex=labels,
    which_lines_black=[0],
    # which_lines_markers=[0],
    transparent_legend=True,legend_border_on=False,grid_lines_on=False,
    clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store
    # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
    )
title_label = "TGV at $t=9.0$"
figure_filename = "spectra_t9_all_cPlus_l2roe"
x=[];y=[];labels=[];
spectra = np.loadtxt("/Users/Julien/PHiLiP-Post-Processing/cases/taylor_green_vortex/data/carton2014_dns_spectra_t9.txt",skiprows=1,delimiter=',')
append_to_plot(spectra[:,0],spectra[:,1],"DNS [Carton]")
batch_append_to_plot(batch_paths, batch_labels, "flow_field_files/velocity_vorticity-1_reordered_spectra.dat")
qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
    # title_label=title_label,
    fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
    nlegendcols=1,
    xlimits=[2e0,1e2],ylimits=[1e-6,1e-1],
    markers=False,legend_on=True,legend_labels_tex=labels,
    which_lines_black=[0],
    # which_lines_markers=[0],
    transparent_legend=True,legend_border_on=False,grid_lines_on=False,
    clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store
    # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]
    )