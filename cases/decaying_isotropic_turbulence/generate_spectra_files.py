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
#-----------------------------------------------------
#=====================================================
# Helper functions
#=====================================================
def generate_spectra_file_from_flow_field_file(file,label,n_skiprows=1):
    print("Loading file: %s" % s)
    data = np.loadtxt(file,skiprows=n_skiprows,dtype=np.float64)
    print("done.")

    print()
    velocity_field = [data[:,3],data[:,4],data[:,5]] # u,v,w
    spectra = es.compute_Ek_spectrum(velocity_field=velocity_field)
    x.append(spectra[0])
    y.append(spectra[1])
    labels.append(label)
#-----------------------------------------------------
