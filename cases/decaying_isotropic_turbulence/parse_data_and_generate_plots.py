#-----------------------------------------------------
# Import public libraries
import numpy as np # NumPy: contains basic numerical routines
from scipy.interpolate import interp1d
#-----------------------------------------------------
import sys
# load tools
sys.path.append("../../src/tools");
from assemble_mpi_flow_field_files_and_reorder import assemble_mpi_flow_field_files_and_reorder
from generate_spectra_files import generate_spectra_file_from_flow_field_file
# load submodules
sys.path.append("../../submodules/quickplotlib/lib"); import quickplotlib as qp
#-----------------------------------------------------
#=====================================================
# Global variables
#=====================================================
global x,y,labels
x=[];y=[];labels=[];
#=====================================================
# Helper functions
#=====================================================
#-----------------------------------------------------
def reinit_inputs():
    global x,y,labels
    x=[];y=[];labels=[];
    return
#-----------------------------------------------------
def append_to_plot(x_,y_,label_):
    global x,y,labels
    labels.append(label_);x.append(x_);y.append(y_)
#-----------------------------------------------------
#=====================================================

#=====================================================
# Store reference spectras
#=====================================================
# (1) Original spectra file
spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/energy.prf",skiprows=1,dtype=np.float64)
input_spectra_to_f77_code = 1.0*spectra

# - Misra and Lund non-dimensionalization:
M = 5.08 # [cm] (mesh size from experiment)
u_rms = 22.2 # [cm/s] (rms velocity from experiment)
L_ref = 11.0*M/(2.0*np.pi) # cm
U_ref = np.sqrt(3.0/2.0)*u_rms # cm/s
energy_ref = U_ref*U_ref*L_ref # cm3/s2

# (2) CBC Experiment data; t=0
spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/cbc_experiment_table3.txt",skiprows=3,usecols=(0,1),dtype=np.float64)
# -- nondimensionalize experiment values
spectra[:,0] *= L_ref # non-dimensionalize wavenumber
spectra[:,1] /= energy_ref # non-dimensionalize energy
cbc_spectra_t0 = 1.0*spectra # store

# (2) CBC Experiment data; t=1
spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/cbc_experiment_table3.txt",skiprows=3,usecols=(0,2),dtype=np.float64)
# -- nondimensionalize experiment values
spectra[:,0] *= L_ref # non-dimensionalize wavenumber
spectra[:,1] /= energy_ref # non-dimensionalize energy
cbc_spectra_t1 = 1.0*spectra # store

# (2) CBC Experiment data; t=1
spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/cbc_experiment_table3.txt",skiprows=2,max_rows=17,usecols=(0,3),dtype=np.float64)
# -- nondimensionalize experiment values
spectra[:,0] *= L_ref # non-dimensionalize wavenumber
spectra[:,1] /= energy_ref # non-dimensionalize energy
cbc_spectra_t2 = 1.0*spectra # store

# # =====================================================
# # 24 DOF
# # =====================================================
# title_label = "DHIT, $24^{3}$ DOF, (P5, $N_{el}=4^{3}$)"
# figure_filename = "spectra_024dof"

# append_to_plot(input_spectra_to_f77_code[:,0],input_spectra_to_f77_code[:,1],"Input to f77 code")
# append_to_plot(cbc_spectra_t0[:,0],cbc_spectra_t0[:,1],"CBC t=0")
# append_to_plot(cbc_spectra_t1[:,0],cbc_spectra_t1[:,1],"CBC t=1")
# append_to_plot(cbc_spectra_t2[:,0],cbc_spectra_t2[:,1],"CBC t=2")

# # generate_spectra_file_from_flow_field_file("/Users/Julien/DHIT-Flow-Setup/dofs024_p5_velocity/flow_field_files/velocity_vorticity-0_reordered_equidistant","dat",n_skiprows=0,use_TurboGenPy=True)
# spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/dofs024_p5_velocity/flow_field_files/velocity_vorticity-0_reordered_equidistant_spectra.dat")
# append_to_plot(spectra[:,0],spectra[:,1],"PHiLiP t=0.0")

# # generate_spectra_file_from_flow_field_file("/Users/Julien/DHIT-Flow-Setup/dofs024_p5_velocity/flow_field_files/velocity_vorticity-1_reordered_equidistant","dat",n_skiprows=0,use_TurboGenPy=True)
# spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/dofs024_p5_velocity/flow_field_files/velocity_vorticity-1_reordered_equidistant_spectra.dat")
# append_to_plot(spectra[:,0],spectra[:,1],"PHiLiP t=1.0")

# # generate_spectra_file_from_flow_field_file("/Users/Julien/DHIT-Flow-Setup/dofs024_p5_velocity/flow_field_files/velocity_vorticity-2_reordered_equidistant","dat",n_skiprows=0,use_TurboGenPy=True)
# spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/dofs024_p5_velocity/flow_field_files/velocity_vorticity-2_reordered_equidistant_spectra.dat")
# append_to_plot(spectra[:,0],spectra[:,1],"PHiLiP t=2.0")

# qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
#     title_label=title_label,
#     fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
#     xlimits=[8e-1,3e2],ylimits=[1e-6,6e-1],
#     markers=False,legend_on=True,legend_labels_tex=labels,
#     which_lines_only_markers=[1,2,3],
#     which_lines_dashed=[0])


# # =====================================================
# # 128 DOF check
# # =====================================================
# reinit_inputs()
# title_label = "DHIT, $128^{3}$ DOF, (P3, $N_{el}=32^{3}$)"
# figure_filename = "spectra_128"

# append_to_plot(input_spectra_to_f77_code[:,0],input_spectra_to_f77_code[:,1],"Input to f77 code")
# append_to_plot(cbc_spectra_t0[:,0],cbc_spectra_t0[:,1],"CBC t=0")
# append_to_plot(cbc_spectra_t1[:,0],cbc_spectra_t1[:,1],"CBC t=1")

# # generate_spectra_file_from_flow_field_file("/Users/Julien/DHIT-Flow-Setup/dofs048_p5_velocity/velocity_equidistant_nodes","fld",n_skiprows=0,use_TurboGenPy=True)
# # spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/dofs048_p5_velocity/velocity_equidistant_nodes_spectra.fld")
# # append_to_plot(spectra[:,0],spectra[:,1],"PHiLiP input")

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
# spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/dofs128_p3_velocity/velocity_equidistant_nodes_spectra.fld")
# append_to_plot(spectra[:,0],spectra[:,1],"$128^3$ DOF eq t=0")

# # NOTE: CHANGE THIS TO EQUIDISTANT ONCE FILE IS AVAILABLE
# spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/dofs128_p3_velocity/flow_field_files/velocity_vorticity-0_reordered_spectra.dat")
# append_to_plot(spectra[:,0],spectra[:,1],"$128^3$ DOF gl t=0")

# # NOTE: CHANGE THIS TO EQUIDISTANT ONCE FILE IS AVAILABLE
# spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/dofs128_p3_velocity/flow_field_files/velocity_vorticity-1_reordered_spectra.dat")
# append_to_plot(spectra[:,0],spectra[:,1],"$128^3$ DOF gl t=0.5")

# qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
#     title_label=title_label,
#     fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
#     xlimits=[8e-1,3e2],ylimits=[1e-6,6e-1],
#     markers=False,legend_on=True,legend_labels_tex=labels,
#     which_lines_only_markers=[1,2],
#     which_lines_dashed=[0])

# =====================================================
# Initialization
# =====================================================
reinit_inputs()
title_label = "DHIT Spectra at $t=0$, $c_{DG}$ Ismail-Roe"
figure_filename = "spectra_initialization"

# append_to_plot(input_spectra_to_f77_code[:,0],input_spectra_to_f77_code[:,1],"Input to f77 code")
append_to_plot(cbc_spectra_t0[:,0],cbc_spectra_t0[:,1],"CBC Experiment")

spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/dofs128_p3_velocity/velocity_equidistant_nodes_spectra.fld")
append_to_plot(spectra[:,0],spectra[:,1],"$128^3$ DOFs")

spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/dofs048_p5_velocity/velocity_equidistant_nodes_spectra.fld")
append_to_plot(spectra[:,0],spectra[:,1],"$48^3$ DOFs")

qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
    fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
    black_lines=False,
    xlimits=[1e0,2e2],ylimits=[1e-4,3e-1],
    markers=False,legend_on=True,legend_labels_tex=labels,
    transparent_legend=True,legend_border_on=False,grid_lines_on=False,
    which_lines_only_markers=[0],
    which_lines_black=[0])

# fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
#     nlegendcols=1,
#     xlimits=[1e0,2e2],ylimits=[1e-4,3e-1],
#     markers=False,legend_on=True,legend_labels_tex=labels,
#     which_lines_black=[0],
#     # which_lines_markers=[0],
#     transparent_legend=True,legend_border_on=False,grid_lines_on=False,
#     clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store)

# =====================================================
# 128 DOF check
# =====================================================
reinit_inputs()
title_label = "DHIT Spectra at $t=0.5$"
title_label = " "
figure_filename = "spectra_comparison"

# append_to_plot(input_spectra_to_f77_code[:,0],input_spectra_to_f77_code[:,1],"Input to f77 code")
# append_to_plot(cbc_spectra_t0[:,0],cbc_spectra_t0[:,1],"CBC t=0")
# append_to_plot(cbc_spectra_t1[:,0],cbc_spectra_t1[:,1],"CBC t=1")

# interpolated CBC spectra
# cbc_spectra_t0_polyfit = interp1d(cbc_spectra_t0[:,0], cbc_spectra_t0[:,1], kind='cubic')
# cbc_spectra_t1_polyfit = interp1d(cbc_spectra_t1[:,0], cbc_spectra_t1[:,1], kind='cubic')
# xnew = np.linspace(cbc_spectra_t0[0,0],cbc_spectra_t0[-1,0],1000)
# append_to_plot(xnew,cbc_spectra_t0_polyfit(xnew),"CBC t=0 fitted")
# append_to_plot(xnew,cbc_spectra_t1_polyfit(xnew),"CBC t=1 fitted")
# append_to_plot(xnew,0.5*(1.0*cbc_spectra_t0_polyfit(xnew)+1.0*cbc_spectra_t1_polyfit(xnew)),"CBC t=0.5 (interpolated)")

# generate_spectra_file_from_flow_field_file("/Users/Julien/DHIT-Flow-Setup/dofs048_p5_velocity/velocity_equidistant_nodes","fld",n_skiprows=0,use_TurboGenPy=True)
# spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/dofs048_p5_velocity/velocity_equidistant_nodes_spectra.fld")
# append_to_plot(spectra[:,0],spectra[:,1],"PHiLiP input")

# spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/dofs048_p5_velocity/flow_field_files/velocity_vorticity-0_reordered_spectra.dat",skiprows=0,dtype=np.float64)
# append_to_plot(spectra[:,0],spectra[:,1],"PHiLiP t=0.00")

spectra = np.loadtxt("/Users/Julien/NarvalFiles/2023_AIAA/2022-11-24_DHIT_128dofs/viscous_DHIT_ILES_cDG_IR_two_point_flux_dofs0128_p3_procs512/flow_field_files/velocity_vorticity-1_reordered_equidistant_spectra.dat",skiprows=0,dtype=np.float64)
# spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/dofs128_p3_velocity/flow_field_files/velocity_vorticity-1_reordered_equidistant_spectra.dat",skiprows=0,dtype=np.float64)
# append_to_plot(spectra[:,0],spectra[:,1],"$128^3$ DOF, (P3, $N_{el}=32^{3}$)")
append_to_plot(spectra[:,0],spectra[:,1],"$128^3$, $c_{DG}$ NSFR")


# spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/dofs048_p5_velocity/flow_field_files/velocity_vorticity-1_reordered_equidistant_spectra.dat",skiprows=0,dtype=np.float64)
# append_to_plot(spectra[:,0],spectra[:,1],"$48^3$ DOF, (P5, $N_{el}=8^{3}$)")

spectra = np.loadtxt("/Users/Julien/NarvalFiles/2023_AIAA/2022-11-24_DHIT_048dofs/viscous_DHIT_ILES_cDG_IR_two_point_flux_dofs048_p5_procs64/flow_field_files/velocity_vorticity-1_reordered_equidistant_spectra.dat",skiprows=0,dtype=np.float64)
# spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/dofs048_p5_velocity/flow_field_files/velocity_vorticity-1_reordered_equidistant_spectra.dat",skiprows=0,dtype=np.float64)
# append_to_plot(spectra[:,0],spectra[:,1],"$48^3$ DOF, (P5, $N_{el}=8^{3}$)")
append_to_plot(spectra[:,0],spectra[:,1],"$48^3$, $c_{DG}$ NSFR")

# generate_spectra_file_from_flow_field_file("/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_DHIT_048dofs/viscous_DHIT_ILES_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64/flow_field_files/velocity_vorticity-1_reordered_equidistant","dat",n_skiprows=0,use_TurboGenPy=True)
# spectra = np.loadtxt("/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_DHIT_048dofs/viscous_DHIT_ILES_cDG_IR_two_point_flux_with_l2roe_dissipation_dofs048_p5_procs64/flow_field_files/velocity_vorticity-1_reordered_equidistant_spectra.dat",skiprows=0,dtype=np.float64)
# append_to_plot(spectra[:,0],spectra[:,1],"$48^3$ DOFs P5 $c_{DG}$ IR+L$^{2}$Roe")

# generate_spectra_file_from_flow_field_file("/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_DHIT_048dofs/viscous_DHIT_ILES_cPlus_IR_two_point_flux_dofs048_p5_procs64/flow_field_files/velocity_vorticity-1_reordered_equidistant","dat",n_skiprows=0,use_TurboGenPy=True)
spectra = np.loadtxt("/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_DHIT_048dofs/viscous_DHIT_ILES_cPlus_IR_two_point_flux_dofs048_p5_procs64/flow_field_files/velocity_vorticity-1_reordered_equidistant_spectra.dat",skiprows=0,dtype=np.float64)
append_to_plot(spectra[:,0],spectra[:,1],"$48^3$, $c_{+}$ NSFR")

# generate_spectra_file_from_flow_field_file("/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_DHIT_048dofs/viscous_DHIT_LES_smagorinsky_cDG_IR_two_point_flux_dofs048_p5_procs64/flow_field_files/velocity_vorticity-1_reordered_equidistant","dat",n_skiprows=0,use_TurboGenPy=True)
# spectra = np.loadtxt("/Users/Julien/NarvalFiles/2023_AIAA/2022-11-29_DHIT_048dofs/viscous_DHIT_LES_smagorinsky_cDG_IR_two_point_flux_dofs048_p5_procs64/flow_field_files/velocity_vorticity-1_reordered_equidistant_spectra.dat",skiprows=0,dtype=np.float64)
# append_to_plot(spectra[:,0],spectra[:,1],"$48^3$ DOFs P5 $c_{DG}$ IR+Smag.SGS")
clr_input_store = ['k','tab:blue','tab:red','tab:red','tab:red','tab:green','tab:green']#,'tab:orange','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
mrkr_input_store = ['None','None','None','None','None','None','None']
lnstl_input_store = ['solid','dashed','dashed','dashed','solid','dashed','solid']
qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
    fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
    nlegendcols=1,
    xlimits=[1e0,2e2],ylimits=[1e-4,3e-1],
    markers=False,legend_on=True,legend_labels_tex=labels,
    which_lines_black=[0],
    # which_lines_markers=[0],
    transparent_legend=True,legend_border_on=False,grid_lines_on=False,
    clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store)
    # which_lines_only_markers=[1,2,3],which_lines_dashed=[0]


#=====================================================
# Plotting function -- Spectra
#=====================================================



exit()
#=====================================================
# Plotting function -- Turbulent Quantities
#=====================================================
sys.path.append("/Users/Julien/PHiLiP-Post-Processing/src");
from plot_unsteady_integrated_turbulent_flow_quantities import get_dissipation_discrete

unsteady_data = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/dofs048_p5_velocity/turbulent_quantities.txt",skiprows=1)
qp.plotfxn(xdata=unsteady_data[:,0],ydata=unsteady_data[:,1],xlabel="Time",ylabel="Integrated KE",
    title_label="DHIT $48^{3}$DOFs, $p5$, cDG Ismail-Roe",
    fig_directory="figures",figure_filename="kinetic_energy_vs_time",log_axes=None,figure_filetype="pdf",
    # xlimits=[8e-1,3e2],ylimits=[1e-6,5e-1],
    markers=False,legend_on=False)
qp.plotfxn(xdata=unsteady_data[:,0],ydata=get_dissipation_discrete(unsteady_data[:,0],unsteady_data[:,1]),xlabel="Time",ylabel="Dissipation Rate",
    title_label="DHIT $48^{3}$DOFs, $p5$, cDG Ismail-Roe",
    fig_directory="figures",figure_filename="dissipation_rate_vs_time",log_axes=None,figure_filetype="pdf",
    # xlimits=[8e-1,3e2],
    ylimits=[-100,100],
    markers=False,legend_on=False)