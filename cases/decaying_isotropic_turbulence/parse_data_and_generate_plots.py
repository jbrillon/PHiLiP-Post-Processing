#-----------------------------------------------------
# Import public libraries
import numpy as np # NumPy: contains basic numerical routines
from scipy.interpolate import interp1d
#-----------------------------------------------------
import os;CURRENT_PATH = os.path.split(os.path.realpath(__file__))[0]+"/";
import sys
# load tools
sys.path.append(CURRENT_PATH+"../../src/tools");
from assemble_mpi_flow_field_files_and_reorder import assemble_mpi_flow_field_files_and_reorder
from generate_spectra_files import generate_spectra_file_from_flow_field_file
# load submodules
sys.path.append(CURRENT_PATH+"../../submodules/quickplotlib/lib"); import quickplotlib as qp
#-----------------------------------------------------
from sys import platform
if platform == "linux" or platform == "linux2":
    # linux
    filesystem="/media/julien/Samsung_T5/"
elif platform == "darwin":
    # OS X
    filesystem="/Volumes/Samsung_T5/"
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
def get_dissipation_discrete(time,kinetic_energy,smoothing=False):
    if(smoothing):
        # smoothing if noisy
        dissipation_rate_val = splrep(time,kinetic_energy,k=5,s=0.0000001)
        dissipation_rate_val = -splev(time,dissipation_rate_val,der=1)
        return dissipation_rate_val
    else:
        dissipation_rate_val = -np.gradient(kinetic_energy,time)
        return dissipation_rate_val
#=====================================================


filename=filesystem+"NarvalFiles/2023_JCP/DHIT/viscous_DHIT_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs128_p3_CFL-0.2_procs512/turbulent_quantities.txt"
time, kinetic_energy, enstrophy, vorticity_based_dissipation, pressure_dilatation_based_dissipation, strain_rate_based_dissipation, deviatoric_strain_rate_based_dissipation = np.loadtxt(filename,skiprows=1,dtype=np.float64,unpack=True)

normalized_kinetic_energy = kinetic_energy/kinetic_energy[0]
# normalized_kinetic_energy=get_dissipation_discrete(time,kinetic_energy,smoothing=False)
labels=[]
xdata=[]
ydata=[]
xdata.append(time)
ydata.append(normalized_kinetic_energy)
# ydata.append(kinetic_energy)
labels.append("NSFR")

# compute reference curve 1
index_of_reference_curve = 1
x_ref_curve = np.linspace(time[0],time[-1],100)
order_for_ref_curve = -1.2
ref_curve_label = "$\\left(t^{*}\\right)^{-1.2}$"
shift = 1.0
y_ref_curve = (x_ref_curve**(order_for_ref_curve))/np.exp(shift)
# clr_input_store.insert(index_of_reference_curve,"k")#"tab:gray"
# lnstl_input_store.insert(index_of_reference_curve,"dotted")

# xdata.append(x_ref_curve)
# ydata.append(y_ref_curve)
# labels.append(ref_curve_label)

title_label = "DHIT, $128^{3}$ DOF, P3, $c_{DG}$ NSFR.IR-GLL, CFL$=0.2$"
qp.plotfxn(xdata=xdata,ydata=ydata,
    xlabel='Nondimensional Time, $t^{*}$',
    ylabel='Nondimensional Kinetic Energy, $K^{*}$',
    title_label=title_label,
    fig_directory="figures",
    figure_filename='kinetic_energy_vs_time',log_axes="both",figure_filetype="pdf",
    # xlimits=[8e-1,3e2],ylimits=[1e-6,6e-1],
    # xlimits=[0.0,2.0],
    # xlimits=[np.amin(time),np.amax(time)],
    # ylimits=[np.amin(normalized_kinetic_energy),np.amax(normalized_kinetic_energy)],
    markers=False,
    legend_on=True,legend_labels_tex=labels,
    # which_lines_only_markers=[0],
    which_lines_dashed=[1],
    # nlegendcols=2,
    transparent_legend=True,
    legend_border_on=False,
    grid_lines_on=False)

exit()

labels=[]
xdata=[]
ydata=[]
xdata.append(time)
ydata.append(enstrophy)
# ydata.append(kinetic_energy)
labels.append("NSFR")

qp.plotfxn(xdata=xdata,ydata=ydata,
    xlabel='Nondimensional Time, $t^{*}$',
    ylabel='Nondimensional Enstrophy, $\\zeta^{*}$',
    title_label=title_label,
    fig_directory="figures",
    figure_filename='enstrophy_vs_time',log_axes="both",figure_filetype="pdf",
    # xlimits=[8e-1,3e2],ylimits=[1e-6,6e-1],
    # xlimits=[0.0,2.0],
    # xlimits=[np.amin(time),np.amax(time)],
    # ylimits=[np.amin(normalized_kinetic_energy),np.amax(normalized_kinetic_energy)],
    markers=False,
    legend_on=True,legend_labels_tex=labels,
    # which_lines_only_markers=[0],
    which_lines_dashed=[1],
    # nlegendcols=2,
    transparent_legend=True,
    legend_border_on=False,
    grid_lines_on=False)

labels=[]
xdata=[]
ydata=[]
xdata.append(time)
ydata.append(pressure_dilatation_based_dissipation)
# ydata.append(kinetic_energy)
labels.append("NSFR")

qp.plotfxn(xdata=xdata,ydata=ydata,
    xlabel='Nondimensional Time, $t^{*}$',
    ylabel='Nondimensional Pressure Dilatation',
    title_label=title_label,
    fig_directory="figures",
    figure_filename='pressure_dilatation_vs_time',log_axes="both",figure_filetype="pdf",
    # xlimits=[8e-1,3e2],ylimits=[1e-6,6e-1],
    # xlimits=[0.0,2.0],
    # xlimits=[np.amin(time),np.amax(time)],
    # ylimits=[np.amin(normalized_kinetic_energy),np.amax(normalized_kinetic_energy)],
    markers=False,
    legend_on=True,legend_labels_tex=labels,
    # which_lines_only_markers=[0],
    which_lines_dashed=[1],
    # nlegendcols=2,
    transparent_legend=True,
    legend_border_on=False,
    grid_lines_on=False)

exit()


#=====================================================
# Store reference spectras
#=====================================================
# (1) Original spectra file
# spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/energy.prf",skiprows=1,dtype=np.float64)
# input_spectra_to_f77_code = 1.0*spectra

# - Misra and Lund non-dimensionalization:
M = 5.08 # [cm] (mesh size from experiment)
u_rms = 22.2 # [cm/s] (rms velocity from experiment)
L_ref = 11.0*M/(2.0*np.pi) # cm
U_ref = np.sqrt(3.0/2.0)*u_rms # cm/s
energy_ref = U_ref*U_ref*L_ref # cm3/s2

# (2) CBC Experiment data; t=0
spectra = np.loadtxt("./data/cbc_experiment_table3.txt",skiprows=3,usecols=(0,1),dtype=np.float64)
# -- nondimensionalize experiment values
spectra[:,0] *= L_ref # non-dimensionalize wavenumber
spectra[:,1] /= energy_ref # non-dimensionalize energy
cbc_spectra_t0 = 1.0*spectra # store

# (2) CBC Experiment data; t=1
spectra = np.loadtxt("./data/cbc_experiment_table3.txt",skiprows=3,usecols=(0,2),dtype=np.float64)
# -- nondimensionalize experiment values
spectra[:,0] *= L_ref # non-dimensionalize wavenumber
spectra[:,1] /= energy_ref # non-dimensionalize energy
cbc_spectra_t1 = 1.0*spectra # store

# (2) CBC Experiment data; t=2
spectra = np.loadtxt("./data/cbc_experiment_table3.txt",skiprows=2,max_rows=17,usecols=(0,3),dtype=np.float64)
# -- nondimensionalize experiment values
spectra[:,0] *= L_ref # non-dimensionalize wavenumber
spectra[:,1] /= energy_ref # non-dimensionalize energy
cbc_spectra_t2 = 1.0*spectra # store

# (3) Vermeire, Nadarajah, and Tucker "Implicit large eddy simulation using the high-order correction procedure via reconstruction scheme", 2016.
vermeire_spectra_p1_126dofs = np.loadtxt("./data/vermeire2016_ref_data/P1.dat",dtype=np.float64)
vermeire_spectra_p2_96dofs = np.loadtxt("./data/vermeire2016_ref_data/P2.dat",dtype=np.float64)
vermeire_spectra_p3_84dofs = np.loadtxt("./data/vermeire2016_ref_data/P3.dat",dtype=np.float64)
vermeire_spectra_p4_80dofs = np.loadtxt("./data/vermeire2016_ref_data/P4.dat",dtype=np.float64)
vermeire_spectra_p5_78dofs = np.loadtxt("./data/vermeire2016_ref_data/P5.dat",dtype=np.float64)
fds_spectra = np.loadtxt("./data/jefferson-loveday_tucker_2010_ref_data/fds_t2_spectra.dat",skiprows=1,delimiter=",",dtype=np.float64)

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


# =====================================================
# 128 DOF check
# =====================================================
reinit_inputs()
title_label = "DHIT, $128^{3}$ DOF, P3, $c_{DG}$ NSFR.IR-GLL, CFL$=0.2$"
figure_filename = "spectra_128"

# append_to_plot(input_spectra_to_f77_code[:,0],input_spectra_to_f77_code[:,1],"Input to f77 code")
append_to_plot(cbc_spectra_t0[:,0],cbc_spectra_t0[:,1],"CBC $t^{*}=0$")
append_to_plot(cbc_spectra_t1[:,0],cbc_spectra_t1[:,1],"CBC $t^{*}=1$")
append_to_plot(cbc_spectra_t2[:,0],cbc_spectra_t2[:,1],"CBC $t^{*}=2$")

# # generate_spectra_file_from_flow_field_file("/Users/Julien/DHIT-Flow-Setup/dofs048_p5_velocity/velocity_equidistant_nodes","fld",n_skiprows=0,use_TurboGenPy=True)
# # spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/dofs048_p5_velocity/velocity_equidistant_nodes_spectra.fld")
# # append_to_plot(spectra[:,0],spectra[:,1],"PHiLiP input")



spectra = np.loadtxt(filesystem+"NarvalFiles/2023_JCP/DHIT/viscous_DHIT_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs128_p3_CFL-0.2_procs512/flow_field_files/velocity_vorticity-0_reordered_spectra_smoothing_on.dat",skiprows=0,dtype=np.float64)
append_to_plot(spectra[:,0],spectra[:,1],"$t^{*}=0$")

# spectra = np.loadtxt("/Volumes/Samsung_T5/NarvalFiles/2023_JCP/DHIT/viscous_DHIT_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs128_p3_CFL-0.2_procs512/flow_field_files/velocity_vorticity-1_reordered_spectra.dat",skiprows=0,dtype=np.float64)
# append_to_plot(spectra[:,0],spectra[:,1],"$c_{DG}$ NSFR.IR-GLL, t=0.5")

spectra = np.loadtxt(filesystem+"NarvalFiles/2023_JCP/DHIT/viscous_DHIT_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs128_p3_CFL-0.2_procs512/flow_field_files/velocity_vorticity-2_reordered_spectra_smoothing_on.dat",skiprows=0,dtype=np.float64)
append_to_plot(spectra[:,0],spectra[:,1],"$t^{*}=1$")

spectra = np.loadtxt(filesystem+"NarvalFiles/2023_JCP/DHIT/viscous_DHIT_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs128_p3_CFL-0.2_procs512/flow_field_files/velocity_vorticity-4_reordered_spectra_smoothing_on.dat",skiprows=0,dtype=np.float64)
append_to_plot(spectra[:,0],spectra[:,1],"$t^{*}=2$")

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

qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
    title_label=title_label,
    fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
    # xlimits=[8e-1,3e2],ylimits=[1e-6,6e-1],
    xlimits=[2.0,3e2],ylimits=[1e-5,1e-1],
    markers=False,legend_on=True,legend_labels_tex=labels,
    which_lines_only_markers=[0,1,2],
    which_lines_dashed=[0,1,2])

qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
    title_label=title_label,
    fig_directory="figures",figure_filename=figure_filename+"-cut",log_axes="both",figure_filetype="pdf",
    # xlimits=[8e-1,3e2],ylimits=[1e-6,6e-1],
    xlimits=[1e0,3.4e1],ylimits=[1e-3,1e-1],
    markers=False,legend_on=True,legend_labels_tex=labels,
    which_lines_only_markers=[0,1,2],
    which_lines_dashed=[0,1,2],
    nlegendcols=2)

# =====================================================
# 128 DOF check
# =====================================================
reinit_inputs()
title_label = "DHIT, $128^{3}$ DOF, P3, $c_{DG}$ NSFR.IR-GLL, CFL$=0.2$"
figure_filename = "spectra_128_no_smoothing_oversampled"

# append_to_plot(input_spectra_to_f77_code[:,0],input_spectra_to_f77_code[:,1],"Input to f77 code")
append_to_plot(cbc_spectra_t0[:,0],cbc_spectra_t0[:,1],"CBC $t^{*}=0$")
append_to_plot(cbc_spectra_t1[:,0],cbc_spectra_t1[:,1],"CBC $t^{*}=1$")
append_to_plot(cbc_spectra_t2[:,0],cbc_spectra_t2[:,1],"CBC $t^{*}=2$")

append_to_plot(fds_spectra[:,0],fds_spectra[:,1],"FDS $t^{*}=2$")
# append_to_plot(vermeire_spectra_p1_126dofs[:,0],vermeire_spectra_p1_126dofs[:,1],"Vermeire P1 $126^3$ DOFs $t^{*}=2$")
append_to_plot(vermeire_spectra_p5_78dofs[:,0],vermeire_spectra_p5_78dofs[:,1],"[Vermeire] CPR P5 $78^3$ DOFs $t^{*}=2$")

# spectra = np.loadtxt(filesystem+"NarvalFiles/2023_JCP/DHIT/viscous_DHIT_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs128_p3_CFL-0.2_procs512/flow_field_files/velocity_vorticity-0_reordered_spectra.dat",skiprows=0,dtype=np.float64)
# append_to_plot(spectra[:,0],spectra[:,1],"$t^{*}=0$ before")

spectra = np.loadtxt(filesystem+"NarvalFiles/2023_JCP/DHIT/viscous_DHIT_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs128_p3_CFL-0.2_procs512_oversampled_nquad12/flow_field_files/velocity_vorticity-0_reordered_spectra_no_smoothing.dat",skiprows=0,dtype=np.float64)
append_to_plot(spectra[:,0],spectra[:,1],"NSFR $t^{*}=0$")

# spectra = np.loadtxt(filesystem+"NarvalFiles/2023_JCP/DHIT/viscous_DHIT_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs128_p3_CFL-0.2_procs512/flow_field_files/velocity_vorticity-2_reordered_spectra.dat",skiprows=0,dtype=np.float64)
# append_to_plot(spectra[:,0],spectra[:,1],"$t^{*}=1$ before")

spectra = np.loadtxt(filesystem+"NarvalFiles/2023_JCP/DHIT/viscous_DHIT_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs128_p3_CFL-0.2_procs512_oversampled_nquad12/flow_field_files/velocity_vorticity-2_reordered_spectra_no_smoothing.dat",skiprows=0,dtype=np.float64)
append_to_plot(spectra[:,0],spectra[:,1],"NSFR $t^{*}=1$")

# spectra = np.loadtxt(filesystem+"NarvalFiles/2023_JCP/DHIT/viscous_DHIT_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs128_p3_CFL-0.2_procs512/flow_field_files/velocity_vorticity-4_reordered_spectra.dat",skiprows=0,dtype=np.float64)
# append_to_plot(spectra[:,0],spectra[:,1],"$t^{*}=2$ before")

spectra = np.loadtxt(filesystem+"NarvalFiles/2023_JCP/DHIT/viscous_DHIT_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs128_p3_CFL-0.2_procs512_oversampled_nquad12/flow_field_files/velocity_vorticity-4_reordered_spectra_no_smoothing.dat",skiprows=0,dtype=np.float64)
append_to_plot(spectra[:,0],spectra[:,1],"NSFR $t^{*}=2$")
# spectra = np.loadtxt(filesystem+"NarvalFiles/2023_JCP/DHIT/viscous_DHIT_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs128_p3_CFL-0.2_procs512_oversampled_nquad12/flow_field_files/velocity_vorticity-4_reordered_spectra_no_smoothing.dat",skiprows=0,dtype=np.float64)
# append_to_plot(spectra[:,0],spectra[:,1],"$t^{*}=2$ nq12")

qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
    title_label=title_label,
    fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
    # xlimits=[8e-1,3e2],ylimits=[1e-6,6e-1],
    # xlimits=[2.0,3e2],ylimits=[1e-5,1e-1],
    # xlimits=[1e0,0.5*96],ylimits=[1e-5,1e-1],
    xlimits=[1e0,3.0e1],ylimits=[1e-5,1e-1],
    markers=False,legend_on=True,legend_labels_tex=labels,
    which_lines_only_markers=[0,1,2,3],
    which_lines_black=[3],
    which_lines_dashed=[0,1,2,3,4])

qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
    title_label=title_label,
    fig_directory="figures",figure_filename=figure_filename+"-cut",log_axes="both",figure_filetype="pdf",
    # xlimits=[8e-1,3e2],ylimits=[1e-6,6e-1],
    xlimits=[1e0,3.4e1],ylimits=[1e-3,1e-1],
    markers=False,legend_on=True,legend_labels_tex=labels,
    which_lines_only_markers=[0,1,2],
    which_lines_dashed=[0,1,2],
    nlegendcols=2)

exit()

# =====================================================
# Comparison to Vermeire et al. 2016 at t=2
# =====================================================
reinit_inputs()
title_label = "DHIT at $t^{*}=2$, $c_{DG}$ NSFR.IR-GLL, CFL$=0.2$"
figure_filename = "spectra_128_t2_comparison"

# append_to_plot(input_spectra_to_f77_code[:,0],input_spectra_to_f77_code[:,1],"Input to f77 code")
append_to_plot(cbc_spectra_t2[:,0],cbc_spectra_t2[:,1],"CBC Experiment")

# # generate_spectra_file_from_flow_field_file("/Users/Julien/DHIT-Flow-Setup/dofs048_p5_velocity/velocity_equidistant_nodes","fld",n_skiprows=0,use_TurboGenPy=True)
# # spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/dofs048_p5_velocity/velocity_equidistant_nodes_spectra.fld")
# # append_to_plot(spectra[:,0],spectra[:,1],"PHiLiP input")

append_to_plot(vermeire_spectra_p1_126dofs[:,0],vermeire_spectra_p1_126dofs[:,1],"Vermeire P1 ($126^3$ DOFs)")
append_to_plot(vermeire_spectra_p2_96dofs[:,0],vermeire_spectra_p2_96dofs[:,1],"Vermeire P2 ($96^3$ DOFs)")
append_to_plot(vermeire_spectra_p3_84dofs[:,0],vermeire_spectra_p3_84dofs[:,1],"Vermeire P3 ($84^3$ DOFs)")
append_to_plot(vermeire_spectra_p4_80dofs[:,0],vermeire_spectra_p4_80dofs[:,1],"Vermeire P4 ($80^3$ DOFs)")
append_to_plot(vermeire_spectra_p5_78dofs[:,0],vermeire_spectra_p5_78dofs[:,1],"Vermeire P5 ($78^3$ DOFs)")

spectra = np.loadtxt(filesystem+"NarvalFiles/2023_JCP/DHIT/viscous_DHIT_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs128_p3_CFL-0.2_procs512/flow_field_files/velocity_vorticity-4_reordered_spectra_smoothing_on.dat",skiprows=0,dtype=np.float64)
append_to_plot(spectra[:,0],spectra[:,1],"P3 ($128^3$ DOFs) smoothing on")
spectra = np.loadtxt(filesystem+"NarvalFiles/2023_JCP/DHIT/viscous_DHIT_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs128_p3_CFL-0.2_procs512/flow_field_files/velocity_vorticity-4_reordered_spectra.dat",skiprows=0,dtype=np.float64)
append_to_plot(spectra[:,0],spectra[:,1],"P3 ($128^3$ DOFs) smoothing off")

qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
    title_label=title_label,
    fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
    # xlimits=[8e-1,3e2],ylimits=[1e-6,6e-1],
    xlimits=[1e0,3.4e1],ylimits=[1e-4,1e-1],
    markers=False,legend_on=True,legend_labels_tex=labels,
    which_lines_only_markers=[0],
    which_lines_dashed=[0,1,2,3,4,5],
    transparent_legend=True,
    legend_border_on=False,grid_lines_on=False)

qp.plotfxn(xdata=x,ydata=y,xlabel="Nondimensional Wavenumber, $k^{*}$",ylabel="Nondimensional TKE Spectra, $E^{*}(k^{*},t^{*})$",
    title_label=title_label,
    fig_directory="figures",figure_filename=figure_filename+"-cut",log_axes="both",figure_filetype="pdf",
    # xlimits=[8e-1,3e2],ylimits=[1e-6,6e-1],
    xlimits=[1e0,3.4e1],ylimits=[1e-3,1e-1],
    markers=False,legend_on=True,legend_labels_tex=labels,
    which_lines_only_markers=[0],
    which_lines_dashed=[0,1,2,3,4,5],
    nlegendcols=2,
    transparent_legend=True,
    legend_border_on=False,grid_lines_on=False)

# qp.plotfxn(xdata=time_store,#[time,time],
#             ydata=kinetic_energy_store,#[kinetic_energy,kolmogorov_slope],
#                 ylabel='Nondimensional Kinetic Energy, $K^{*}$',#=\\frac{1}{\\rho_{\\infty}V_{\\infty}^{2}|\\Omega|}\\int_{\\Omega}\\rho(u\\cdot\\u)d\\Omega$',
#                 xlabel='Nondimensional Time, $t^{*}$',
#                 figure_filename=,
#                 title_label=figure_title,
#                 markers=False,
#                 legend_labels_tex=labels_store,
#                 black_lines=False,
#                 xlimits=[0,tmax],
#                 ylimits=[0.0,0.14],
#                 log_axes=log_axes_input,
#                 which_lines_black=which_lines_black_input,
#                 which_lines_dashed=which_lines_dashed_input,
#                 legend_on=legend_on_input,
#                 legend_inside=legend_inside_input,
#                 nlegendcols=nlegendcols_input,
#                 figure_size=(6,6),
#                 transparent_legend=True,#transparent_legend_input,
#                 legend_border_on=False,
#                 grid_lines_on=False,
#                 clr_input=clr_input_store,mrkr_input=mrkr_input_store,lnstl_input=lnstl_input_store,
#                 legend_fontSize=legend_fontSize_input,
#                 legend_location="lower left")#,

exit()
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