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

#=====================================================
# Global variables
#=====================================================
global x,y,labels
x=[];y=[];labels=[];


    # title_label = "DHIT Initialization Check: TurboGenPy"#\n P3, $N_{el}=32^{3}$ ($128^{3}$ DOF)"
    # figure_filename = "spectra_turbogenpy"

    # # Original spectra file
    # spectra = np.loadtxt("/Users/Julien/TurboGenPy/cbc_spectrum.txt",skiprows=0,max_rows=21,usecols=(0,1),dtype=np.float64)
    # x.append(spectra[:,0]*100.0)
    # y.append(spectra[:,1]/100.0/100.0/100.0)
    # labels.append("CBC Exp. t=0.0")

    # # Original spectra file
    # spectra = np.loadtxt("/Users/Julien/TurboGenPy/tkespec_cbc_32.32.32_5000_modes.txt",skiprows=0,dtype=np.float64)
    # x.append(spectra[:,0])
    # y.append(spectra[:,1])
    # labels.append("TurboGenPY Spec")

    # # computed from the outputted velocity field
    # file = "/Users/Julien/TurboGenPy/vel_cbc_32.32.32_5000_modes.txt"
    # load_velocity_field_compute_spectra_append_to_plot(
    #     file=file,
    #     label="Computed",
    #     n_skiprows=1,
    #     use_TurboGenPy=True)

    # #=====================================================
    # # Plotting function
    # #=====================================================
    # qp.plotfxn(xdata=x,ydata=y,xlabel="$k$",ylabel="$E(k)$",
    #     title_label=title_label,
    #     fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
    #     xlimits=[1e1,1e3],ylimits=[1e-6,1e-3],
    #     markers=True,legend_on=True,legend_labels_tex=labels)
    # exit()

#=====================================================
# Add reference spectras
#=====================================================
# (1) Original spectra file
spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/energy.prf",skiprows=1,dtype=np.float64)
append_to_plot(spectra[:,0],spectra[:,1],"Input to f77 code")

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
# add to plot
append_to_plot(spectra[:,0],spectra[:,1],"CBC t=0")

# (2) CBC Experiment data; t=1
spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/cbc_experiment_table3.txt",skiprows=3,usecols=(0,2),dtype=np.float64)
# -- nondimensionalize experiment values
spectra[:,0] *= L_ref # non-dimensionalize wavenumber
spectra[:,1] /= energy_ref # non-dimensionalize energy
# add to plot
append_to_plot(spectra[:,0],spectra[:,1],"CBC t=1")

# (2) CBC Experiment data; t=1
spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/cbc_experiment_table3.txt",skiprows=2,max_rows=17,usecols=(0,3),dtype=np.float64)
# -- nondimensionalize experiment values
spectra[:,0] *= L_ref # non-dimensionalize wavenumber
spectra[:,1] /= energy_ref # non-dimensionalize energy
# add to plot
append_to_plot(spectra[:,0],spectra[:,1],"CBC t=2")

# =====================================================
# 24 DOF check
# =====================================================
title_label = "DHIT, $24^{3}$ DOF, (P5, $N_{el}=4^{3}$)"
figure_filename = "spectra_024dof_after_bug_fix"
# load files

# generate_spectra_file_from_flow_field_file("/Users/Julien/DHIT-Flow-Setup/dofs048_p5_velocity/velocity_equidistant_nodes","fld",n_skiprows=0,use_TurboGenPy=True)
# spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/dofs048_p5_velocity/velocity_equidistant_nodes_spectra.fld")
# append_to_plot(spectra[:,0],spectra[:,1],"PHiLiP input")

generate_spectra_file_from_flow_field_file("/Users/Julien/DHIT-Flow-Setup/dofs024_p5_velocity/flow_field_files/velocity_vorticity-0_reordered_equidistant","dat",n_skiprows=0,use_TurboGenPy=True)
spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/dofs024_p5_velocity/flow_field_files/velocity_vorticity-0_reordered_equidistant_spectra.dat")
append_to_plot(spectra[:,0],spectra[:,1],"PHiLiP t=0.0")

generate_spectra_file_from_flow_field_file("/Users/Julien/DHIT-Flow-Setup/dofs024_p5_velocity/flow_field_files/velocity_vorticity-1_reordered_equidistant","dat",n_skiprows=0,use_TurboGenPy=True)
spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/dofs024_p5_velocity/flow_field_files/velocity_vorticity-1_reordered_equidistant_spectra.dat")
append_to_plot(spectra[:,0],spectra[:,1],"PHiLiP t=1.0")

generate_spectra_file_from_flow_field_file("/Users/Julien/DHIT-Flow-Setup/dofs024_p5_velocity/flow_field_files/velocity_vorticity-2_reordered_equidistant","dat",n_skiprows=0,use_TurboGenPy=True)
spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/dofs024_p5_velocity/flow_field_files/velocity_vorticity-2_reordered_equidistant_spectra.dat")
append_to_plot(spectra[:,0],spectra[:,1],"PHiLiP t=2.0")

# spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/dofs048_p5_velocity/flow_field_files/velocity_vorticity-0_reordered_spectra.dat",skiprows=0,dtype=np.float64)
# append_to_plot(spectra[:,0],spectra[:,1],"PHiLiP t=0.00")

# spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/dofs048_p5_velocity/flow_field_files/velocity_vorticity-1_reordered_spectra.dat",skiprows=0,dtype=np.float64)
# append_to_plot(spectra[:,0],spectra[:,1],"PHiLiP t=0.50")

# spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/dofs048_p5_velocity/flow_field_files/velocity_vorticity-2_reordered_spectra.dat",skiprows=0,dtype=np.float64)
# append_to_plot(spectra[:,0],spectra[:,1],"PHiLiP t=0.75")

# generate_spectra_file_from_flow_field_file("/Users/Julien/DHIT-Flow-Setup/dofs048_p5_velocity/velocity_equidistant_nodes","fld",n_skiprows=0,use_TurboGenPy=True)
# spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/dofs048_p5_velocity/velocity_equidistant_nodes_spectra.fld")
# append_to_plot(spectra[:,0],spectra[:,1],"$48^3$ DOF")

# generate_spectra_file_from_flow_field_file("/Users/Julien/DHIT-Flow-Setup/dofs128_p3_velocity/velocity_equidistant_nodes","fld",n_skiprows=0,use_TurboGenPy=True)
# spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/dofs128_p3_velocity/velocity_equidistant_nodes_spectra.fld")
# append_to_plot(spectra[:,0],spectra[:,1],"$128^3$ DOF")

#=====================================================
# Plotting function -- Spectra
#=====================================================
qp.plotfxn(xdata=x,ydata=y,xlabel="$k$",ylabel="$E(k)$",
    title_label=title_label,
    fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
    xlimits=[8e-1,3e2],ylimits=[1e-6,6e-1],
    markers=False,legend_on=True,legend_labels_tex=labels,
    which_lines_only_markers=[1,2,3],
    which_lines_dashed=[0])


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