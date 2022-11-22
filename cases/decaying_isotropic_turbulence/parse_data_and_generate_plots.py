#-----------------------------------------------------
# Import public libraries
import numpy as np # NumPy: contains basic numerical routines
#-----------------------------------------------------
import sys
# load tools
sys.path.append("/Users/Julien/PHiLiP-Post-Processing/src/tools");
from assemble_mpi_flow_field_files_and_reorder import assemble_mpi_flow_field_files_and_reorder
# load submodules
sys.path.append("/Users/Julien/PHiLiP-Post-Processing/submodules/quickplotlib/lib"); import quickplotlib as qp
sys.path.append("/Users/Julien/PHiLiP-Post-Processing/submodules/Energy_Spectrum"); import Energy_Spectrum as es
sys.path.append("/Users/Julien/PHiLiP-Post-Processing/submodules/TurboGenPY"); from tkespec import compute_tke_spectrum
#-----------------------------------------------------
#=====================================================
# Helper functions
#=====================================================
def get_velocity_components_as_3d_arrays_from_velocity_field(velocity_field):
    npts = 32 # number of points per direction
    u = np.zeros((npts,npts,npts))
    v = np.zeros((npts,npts,npts))
    w = np.zeros((npts,npts,npts))
    ig = 0 # global index
    for k in range(0,npts):
        for j in range(0,npts):
            for i in range(0,npts):
                u[i,j,k] = velocity_field[ig,0]
                v[i,j,k] = velocity_field[ig,1]
                w[i,j,k] = velocity_field[ig,2]
                ig += 1
    return [u,v,w]

def load_velocity_field_compute_spectra_append_to_plot(file,label,n_skiprows=1):
    global x,y,labels
    print("Loading file: %s" % file)
    data = np.loadtxt(file,skiprows=n_skiprows,usecols=(3,4,5),dtype=np.float64)
    print("done.")

    velocity_field = [data[:,0],data[:,1],data[:,2]] # u,v,w

    # code from farshad
    # if TGV: use Farshad's code
    # spectra = es.compute_Ek_spectrum(velocity_field=velocity_field)
    # x.append(spectra[0])
    # y.append(spectra[1])

    # code from turbogenpy
    # if DHIT: use TurboGenPY code
    u,v,w = get_velocity_components_as_3d_arrays_from_velocity_field(data)
    lx = 9.0*2.0*np.pi / 100.0
    ly = 9.0*2.0*np.pi / 100.0
    lz = 9.0*2.0*np.pi / 100.0
    smoothing=True # False by default
    knyquist, wavenumbers, tkespec = compute_tke_spectrum(u, v, w, lx, ly, lz, True)
    x.append(wavenumbers)
    y.append(tkespec)

    
    labels.append(label)
#-----------------------------------------------------

#=====================================================
# Global variables
#=====================================================
global x,y,labels
x=[];y=[];labels=[];


title_label = "DHIT Initialization Check: TurboGenPy"#\n P3, $N_{el}=32^{3}$ ($128^{3}$ DOF)"
figure_filename = "spectra_turbogenpy"

# Original spectra file
spectra = np.loadtxt("/Users/Julien/TurboGenPy/cbc_spectrum.txt",skiprows=0,max_rows=21,usecols=(0,1),dtype=np.float64)
x.append(spectra[:,0]*100.0)
y.append(spectra[:,1]/100.0/100.0/100.0)
labels.append("CBC Exp. t=0.0")

# Original spectra file
spectra = np.loadtxt("/Users/Julien/TurboGenPy/tkespec_cbc_32.32.32_5000_modes.txt",skiprows=0,dtype=np.float64)
x.append(spectra[:,0])
y.append(spectra[:,1])
labels.append("TurboGenPY Spec")

# computed from the outputted velocity field
file = "/Users/Julien/TurboGenPy/vel_cbc_32.32.32_5000_modes.txt"
load_velocity_field_compute_spectra_append_to_plot(
    file=file,
    label="Computed",
    n_skiprows=1)

#=====================================================
# Plotting function
#=====================================================
qp.plotfxn(xdata=x,ydata=y,xlabel="$k$",ylabel="$E(k)$",
    title_label=title_label,
    fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
    xlimits=[1e1,1e3],ylimits=[1e-6,1e-3],
    markers=True,legend_on=True,legend_labels_tex=labels)
exit()

#=====================================================
# Add reference spectras
#=====================================================
# Original spectra file
spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/energy.prf",skiprows=1,dtype=np.float64)
x.append(spectra[:,0])
y.append(spectra[:,1])
labels.append("Input to box.for")

# CBC Experiment data
spectra = np.loadtxt("/Users/Julien/DHIT-Flow-Setup/cbc_experiment_table3.txt",skiprows=3,usecols=(0,1),dtype=np.float64)
# Misra and Lund non-dimensionalization:
M = 5.08 # [cm] (mesh size from experiment)
u_rms = 22.2 # [cm/s] (rms velocity from experiment)
L_ref = 11.0*M/(2.0*np.pi) # cm
U_ref = np.sqrt(3.0/2.0)*u_rms # cm/s
energy_ref = U_ref*U_ref*L_ref # cm3/s2
# -- nondimensionalize experiment values
spectra[:,0] *= L_ref # non-dimensionalize wavenumber
spectra[:,1] /= energy_ref # non-dimensionalize energy
# add to plot
x.append(spectra[:,0])
y.append(spectra[:,1])
labels.append("CBC Experiment (scaled)")


#=====================================================
# 24 DOF check
#=====================================================
# title_label = "DHIT Initialization Check\n P5, $N_{el}=4^{3}$ ($24^{3}$ DOF)"
# figure_filename = "spectra_24dof"
# load_velocity_field_append_to_plot(
#     file="/Users/Julien/DHIT-Flow-Setup/velocity_equidistant_nodes.fld",
#     label="PHiLiP Input",
#     n_skiprows=0)

# load_velocity_field_append_to_plot(
#     file="/Users/Julien/DHIT-Flow-Setup/philip_outputs/test/velocity_vorticity-0_reordered_for_spectra.dat",
#     label="PHiLiP Output",
#     n_skiprows=0)

# # for calling it -- as a first test
# file_path_and_prefix = "/Users/Julien/DHIT-Flow-Setup/philip_outputs/test/velocity_vorticity-0"
# file_extension = "dat"
# poly_degree=5
# nElements_per_direction=4
# nValues_per_row=6
# num_procs=4
# assemble_mpi_flow_field_files_and_reorder(file_path_and_prefix,file_extension,num_procs,nValues_per_row,nElements_per_direction,poly_degree)
# file = file_path_and_prefix+"_reordered"+"."+file_extension
# load_velocity_field_append_to_plot(
#     file=file,
#     label="PHiLiP Output 2",
#     n_skiprows=1)


#=====================================================
# 128 DOF
#=====================================================
title_label = "DHIT Initialization Check\n P3, $N_{el}=32^{3}$ ($128^{3}$ DOF)"
figure_filename = "spectra_128dof"

# philip input
load_velocity_field_compute_spectra_append_to_plot(
    file="/Users/Julien/DHIT-Flow-Setup/dofs128_p3_velocity/velocity_equidistant_nodes.fld",
    label="PHiLiP Input",
    n_skiprows=0)

# philip output
file_path = "/Users/Julien/PHiLiP-Post-Processing/cases/decaying_isotropic_turbulence/data/2022-11-11_DHIT_test_128dofs/viscous_DHIT_ILES_cDG_IR_two_point_flux_dofs0128_p3_procs512/"
prefix = "velocity_vorticity-0"
file_extension = "dat"
file_path_and_prefix = file_path+prefix
do_assemble_and_reorder_mpi_files = False # Set as true for first time
if(do_assemble_and_reorder_mpi_files):
    poly_degree=3
    nElements_per_direction=32
    nValues_per_row=7 # with vorticity
    num_procs=512
    assemble_mpi_flow_field_files_and_reorder(file_path_and_prefix,file_extension,num_procs,nValues_per_row,nElements_per_direction,poly_degree)
file = file_path_and_prefix+"_reordered"+"."+file_extension
load_velocity_field_compute_spectra_append_to_plot(
    file=file,
    label="PHiLiP Output",
    n_skiprows=1)


#=====================================================
# Plotting function
#=====================================================
qp.plotfxn(xdata=x,ydata=y,xlabel="$k$",ylabel="$E(k)$",
    title_label=title_label,
    fig_directory="figures",figure_filename=figure_filename,log_axes="both",figure_filetype="pdf",
    xlimits=[8e-1,3e2],ylimits=[1e-6,1e-1],
    markers=True,legend_on=True,legend_labels_tex=labels)


