import os;CURRENT_PATH = os.path.split(os.path.realpath(__file__))[0]+"/";
import sys
sys.path.append(CURRENT_PATH+"../../src/tools");
from batch_assemble_mpi_flow_field_files_reorder_generate_spectra import batch_assemble_mpi_flow_field_files_reorder_generate_spectra_from_txt

file = sys.argv[1] # to pass file from terminal
print("Specified input file for paths: %s" % file)
batch_assemble_mpi_flow_field_files_reorder_generate_spectra_from_txt(file)