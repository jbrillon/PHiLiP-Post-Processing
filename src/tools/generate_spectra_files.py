#-----------------------------------------------------
# Import public libraries
import numpy as np # NumPy: contains basic numerical routines
#-----------------------------------------------------
import sys
# load tools
sys.path.append("/Users/Julien/PHiLiP-Post-Processing/src/tools");
from assemble_mpi_flow_field_files_and_reorder import assemble_mpi_flow_field_files_and_reorder
# load submodules
sys.path.append("/Users/Julien/PHiLiP-Post-Processing/submodules/Energy_Spectrum"); import Energy_Spectrum as es
sys.path.append("/Users/Julien/PHiLiP-Post-Processing/submodules/TurboGenPY"); from tkespec import compute_tke_spectrum
#-----------------------------------------------------
#=====================================================
# Helper functions
#=====================================================
def get_velocity_components_as_3d_arrays_from_velocity_field(velocity_field):
    # determinue number of points per direction
    npts = int(velocity_field[0].shape[0])
    npts = int(round(npts**(1.0/3.0)))
    # build velocity components as 3D arrays
    u = np.zeros((npts,npts,npts))
    v = np.zeros((npts,npts,npts))
    w = np.zeros((npts,npts,npts))
    ig = 0 # global index
    for k in range(0,npts):
        for j in range(0,npts):
            for i in range(0,npts):
                u[i,j,k] = velocity_field[0][ig]
                v[i,j,k] = velocity_field[1][ig]
                w[i,j,k] = velocity_field[2][ig]
                ig += 1
    return [u,v,w]
#-----------------------------------------------------
def get_fluctuating_velocity_field(velocity_field):
    print(" - Getting fluctuating velocity field...")
    components_string=["x","y","z"]
    for i in range(0,3):
        mean_velocity = np.average(velocity_field[i])
        velocity_field[i] -= mean_velocity
        print(" - * Mean velocity in %s-direction: %1.8e" % (components_string[i],mean_velocity))
    print(" - done.")
    return velocity_field
#-----------------------------------------------------  
def get_tke_spectra(velocity_fluctuation_field,use_TurboGenPy):
    print(" - Computing spectra...")
    # compute tke spectra
    if(use_TurboGenPy):
        # code from turbogenpy
        u,v,w = get_velocity_components_as_3d_arrays_from_velocity_field(velocity_fluctuation_field)

        # TO DO: Modify this function call
        # lx = 9.0*2.0*np.pi / 100.0 # [m]
        # ly = 9.0*2.0*np.pi / 100.0 # [m]
        # lz = 9.0*2.0*np.pi / 100.0 # [m]
        lx = 2.0*np.pi
        ly = 2.0*np.pi
        lz = 2.0*np.pi
        use_smoothing=True # False by default
        knyquist, wavenumbers, tkespec = compute_tke_spectrum(u, v, w, lx, ly, lz, use_smoothing)
        
        # write spectra
        npts_spectra = int(wavenumbers.shape[0])
        data_out = np.zeros((npts_spectra,2),dtype=np.float64)
        data_out[:,0] = wavenumbers
        data_out[:,1] = tkespec
    else:
        # code from farshad (default)
        spectra = es.compute_Ek_spectrum(velocity_field=velocity_fluctuation_field)
        # write spectra
        npts_spectra = int(spectra[0].shape[0])
        data_out = np.zeros((npts_spectra,2),dtype=np.float64)
        data_out[:,0] = spectra[0]
        data_out[:,1] = spectra[1]
    print(" - done.")
    return data_out
#-----------------------------------------------------
def generate_spectra_file_from_flow_field_file(
    file_without_extension,
    file_extension,
    n_skiprows=1,
    use_TurboGenPy=False):

    print("Generating spectra file from flow field file...")
    file = file_without_extension+"."+file_extension
    print(" - Loading file: %s" % file)
    data = np.loadtxt(file,skiprows=n_skiprows,usecols=(3,4,5),dtype=np.float64)
    print(" - done.")

    file_out = file_without_extension+"_spectra."+file_extension

    velocity_field = [data[:,0],data[:,1],data[:,2]] # u,v,w
    spectra = get_tke_spectra(get_fluctuating_velocity_field(velocity_field),use_TurboGenPy)
    np.savetxt(file_out,spectra)

    print("done.")
    return
#-----------------------------------------------------
def batch_assemble_mpi_flow_field_files_reorder_generate_spectra(
    file_path=[],
    n_different_files_in_path=[],
    file_prefix=[],
    poly_degree=[],
    nElements_per_direction=[],
    nValues_per_row=[],
    num_procs=[],
    file_extension="dat"):
    #-----------------------------------------------------
    # Safeguard for when empty args are passed
    #-----------------------------------------------------
    if(file_path==[] or file_prefix==[] or poly_degree==[] or nElements_per_direction==[] or num_procs==[]):
        print("batch_assemble_mpi_flow_field_files_and_reorder error: an empty essential argument was passed")
        print("aborting...")
        return
    else:
        n_file_paths = int(len(file_path))
    #-----------------------------------------------------
    # If n_different_files_in_path is empty, assume only
    # one file to assemble per path
    #-----------------------------------------------------
    if(n_different_files_in_path==[]):
        # if empty, assume only one file to assemble per path
        n_files_to_assemble = int(len(file_path))
        for i in range(0,n_files_to_assemble):
            n_different_files_in_path.append(1)
    #-----------------------------------------------------
    # Assemble the files
    #-----------------------------------------------------
    for i in range(0,n_file_paths):
        # loop for multiple prefixes per path
        n_prefixes = n_different_files_in_path[i]
        for j in range(0,n_prefixes):
            # get prefix
            if (n_prefixes>1):
                prefix = file_prefix[i][j]
            else:
                prefix = file_prefix[i]
            
            # get file path and prefix
            file_path_and_prefix = file_path[i] + prefix
            
            # assemble
            assemble_mpi_flow_field_files_and_reorder(
                file_path_and_prefix,
                file_extension,
                num_procs[i],
                nValues_per_row[i],
                nElements_per_direction[i],
                poly_degree[i])
            
            # generate the spectra file
            velocity_file_for_spectra_without_extension = file_path_and_prefix+"_reordered"
            generate_spectra_file_from_flow_field_file(
                velocity_file_for_spectra_without_extension,
                file_extension,
                n_skiprows=1,
                use_TurboGenPy=True)
    return
#-----------------------------------------------------