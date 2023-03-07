import sys
sys.path.append("../../src/tools");
from batch_assemble_mpi_flow_field_files_reorder_generate_spectra import batch_assemble_mpi_flow_field_files_reorder_generate_spectra_from_txt

file = input('Specify input file: ')
batch_assemble_mpi_flow_field_files_reorder_generate_spectra_from_txt(file)