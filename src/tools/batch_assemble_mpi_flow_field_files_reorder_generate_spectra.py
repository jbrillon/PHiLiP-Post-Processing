import sys; sys.path.append("../../src/tools");
from generate_spectra_files import batch_assemble_mpi_flow_field_files_reorder_generate_spectra

#=====================================================
# Input variables
#=====================================================
global file_path_store, file_prefix_store, n_different_files_in_path_store, \
poly_degree_store, nElements_per_direction_store, nValues_per_row_store, \
num_procs_store
# Initialize
file_path_store = []
file_prefix_store = []
n_different_files_in_path_store=[]
poly_degree_store = []
nElements_per_direction_store = []
nValues_per_row_store = []
num_procs_store = []
#=====================================================
# Helper function
#=====================================================
def add_to_batch(
    file_path=" ",
    nValues_per_row=7,
    output_flow_field_file_directory_name="flow_field_files/"):
    
    global file_path_store, file_prefix_store, n_different_files_in_path_store, \
    poly_degree_store, nElements_per_direction_store, nValues_per_row_store, \
    num_procs_store

    file1 = open(file_path + 'parameters_for_assembling_mpi_files.txt', 'r')
    flow_case=file1.readline().rstrip('\n')
    poly_degree=int(file1.readline().rstrip('\n'))
    nElements_per_direction=int(file1.readline().rstrip('\n'))
    num_procs=int(file1.readline().rstrip('\n'))
    file1.close()
    # number of different files in path based on the flow case
    if(flow_case=="DHIT"):
        n_files_per_dir = 2
        # n_files_per_dir = 2 # TESTING
        print("n_files_per_dir = %i" % n_files_per_dir)
    elif(flow_case=="TGV"):
        n_files_per_dir = 2

    # generate the standard prefix files
    prefix_base = "velocity_vorticity-"
    standard_prefix_files = []
    for i in range(0,n_files_per_dir):
        standard_prefix_files.append(prefix_base+str(i))

    file_path_store.append(file_path+output_flow_field_file_directory_name)
    file_prefix_store.append(standard_prefix_files)
    n_different_files_in_path_store.append(n_files_per_dir)
    poly_degree_store.append(poly_degree)
    nElements_per_direction_store.append(nElements_per_direction)
    nValues_per_row_store.append(nValues_per_row)
    num_procs_store.append(num_procs)
    return
#-----------------------------------------------------
def batch_assemble_mpi_flow_field_files_reorder_generate_spectra_from_txt(input_file):
    global file_path_store, file_prefix_store, n_different_files_in_path_store, \
    poly_degree_store, nElements_per_direction_store, nValues_per_row_store, \
    num_procs_store

    # =========================================================
    #                   PATHS FOR BATCH ASSEMBLY
    # =========================================================
    file1 = open(input_file, 'r')
    paths = file1.readlines()
    for path in paths:
        add_to_batch(file_path=path.rstrip('\n'))
    file1.close()
    # =========================================================
    # CALL THE BATCH FUNCTION
    # =========================================================
    batch_assemble_mpi_flow_field_files_reorder_generate_spectra(
        file_path=file_path_store,
        n_different_files_in_path=n_different_files_in_path_store,
        file_prefix=file_prefix_store,
        poly_degree=poly_degree_store,
        nElements_per_direction=nElements_per_direction_store,
        nValues_per_row=nValues_per_row_store,
        num_procs=num_procs_store,
        file_extension="dat")
    return
#-----------------------------------------------------