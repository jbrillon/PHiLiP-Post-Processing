import os;CURRENT_PATH = os.path.split(os.path.realpath(__file__))[0]+"/";
import sys; sys.path.append(CURRENT_PATH+"../../src/tools");
from generate_spectra_files import batch_assemble_mpi_flow_field_files_reorder_generate_spectra
from sys import platform
if platform == "linux" or platform == "linux2":
    # linux
    filesystem="/media/julien/Samsung_T5/"
elif platform == "darwin":
    # OS X
    filesystem="/Volumes/Samsung_T5/"
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
        n_files_per_dir = 5
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
        # add_to_batch(file_path=filesystem+path.rstrip('\n'))
        add_to_batch(file_path=path.rstrip('\n')) # for narval
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
# Converting velocity field to different nodes
#-----------------------------------------------------
from get_DOF_vars import get_DOF_vars
from convert_equidistant_to_gauss_lobatto_nodes import convert_equidistant_to_gauss_lobatto_nodes
#-----------------------------------------------------
def batch_convert_velocity_field_at_equidistant_nodes_to_gLL_nodes_from_txt(input_file):
    global file_path_store, file_prefix_store, n_different_files_in_path_store, \
    poly_degree_store, nElements_per_direction_store, nValues_per_row_store, \
    num_procs_store
    # # =========================================================
    # #                   PATHS FOR BATCH ASSEMBLY
    # # =========================================================
    file1 = open(input_file, 'r')
    paths = file1.readlines()
    # paths = ["NarvalFiles/2023_JCP/robustness/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs048_p5_procs64/"] # for testing
    # paths = ["NarvalFiles/2023_JCP/robustness/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs024_p5_procs16/"]
    # paths = ["NarvalFiles/2023_JCP/filtered_dns_viscous_tgv/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p7_procs1024/"]
    filesystem = "./data_for_spectra_fix/" # FOR NARVAL
    for path in paths:
        # add_to_batch(file_path=filesystem+path.rstrip('\n'))
        # add_to_batch(file_path=filesystem+path) # for testing
        add_to_batch(file_path=filesystem+path.rstrip('\n'))
    file1.close()

    # =========================================================
    # CALL THE BATCH FUNCTION
    # =========================================================
    batch_convert_velocity_field_at_equidistant_nodes_to_gLL_nodes(
        file_path=file_path_store,
        n_different_files_in_path=n_different_files_in_path_store,
        file_prefix=file_prefix_store,
        poly_degree=poly_degree_store,
        nElements_per_direction=nElements_per_direction_store,
        nValues_per_row=nValues_per_row_store,
        num_procs=num_procs_store,
        file_extension="dat")

#-----------------------------------------------------
def batch_convert_velocity_field_at_equidistant_nodes_to_gLL_nodes(
    file_path=[],
    n_different_files_in_path=[],
    file_prefix=[],
    poly_degree=[],
    nElements_per_direction=[],
    nValues_per_row=[],
    num_procs=[],# not needed but could overwrite to be fixed as 8
    file_extension="dat"):
    #-----------------------------------------------------
    # Safeguard for when empty args are passed
    #-----------------------------------------------------
    if(file_path==[] or file_prefix==[] or poly_degree==[] or nElements_per_direction==[] or num_procs==[]):
        print("batch_convert_velocity_field_at_equidistant_nodes_to_gLL_nodes: an empty essential argument was passed")
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
            prefix = file_prefix[i][j]
            
            # get file path and prefix
            file_path_and_prefix = file_path[i] + prefix
            
            # # assemble
            # assemble_mpi_flow_field_files_and_reorder(
            #     file_path_and_prefix,
            #     file_extension,
            #     num_procs[i],
            #     nValues_per_row[i],
            #     nElements_per_direction[i],
            #     poly_degree[i])
            
            # generate the spectra file
            velocity_file_for_spectra_without_extension = file_path_and_prefix+"_reordered"
            velocity_file_for_spectra_with_extension = velocity_file_for_spectra_without_extension+"."+file_extension
            converted_velocity_file_for_spectra_with_extension = velocity_file_for_spectra_without_extension+"_gll_nodes"+"."+file_extension
            #-----------------------------------------------------
            # Get DOF variables
            #-----------------------------------------------------
            nElements,nQuadPoints_per_element,nQuadPoints,nDOF,reduced_nQuadPoints,reduced_nDOF = get_DOF_vars(nElements_per_direction[i],poly_degree[i])
            
            #-----------------------------------------------------
            # Convert from equisdistant to Gauss Lobatto nodes
            #-----------------------------------------------------
            convert_equidistant_to_gauss_lobatto_nodes(
                velocity_file_for_spectra_without_extension+"."+file_extension,
                nElements_per_direction[i],
                nQuadPoints_per_element,
                nValues_per_row[i],
                poly_degree[i],
                nDOF,
                output_filename=converted_velocity_file_for_spectra_with_extension,
                test_reading=False)
            
    return
#-----------------------------------------------------
# For generating input files
#-----------------------------------------------------
def batch_generate_philip_input_files_from_txt(input_file):
    global file_path_store, file_prefix_store, n_different_files_in_path_store, \
    poly_degree_store, nElements_per_direction_store, nValues_per_row_store, \
    num_procs_store
    # # =========================================================
    # #                   PATHS FOR BATCH ASSEMBLY
    # # =========================================================
    file1 = open(input_file, 'r')
    paths = file1.readlines()
    # paths = ["NarvalFiles/2023_JCP/robustness/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs048_p5_procs64/"] # for testing
    # paths = ["NarvalFiles/2023_JCP/robustness/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs024_p5_procs16/"]
    # paths = ["NarvalFiles/2023_JCP/filtered_dns_viscous_tgv/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p7_procs1024/"]
    filesystem = "./data_for_spectra_fix/" # FOR NARVAL
    for path in paths:
        # add_to_batch(file_path=filesystem+path.rstrip('\n'))
        # add_to_batch(file_path=filesystem+path) # for testing
        add_to_batch(file_path=filesystem+path.rstrip('\n'))
    file1.close()

    # =========================================================
    # CALL THE BATCH FUNCTION
    # =========================================================
    batch_generate_philip_input_files(
        file_path=file_path_store,
        n_different_files_in_path=n_different_files_in_path_store,
        file_prefix=file_prefix_store,
        poly_degree=poly_degree_store,
        nElements_per_direction=nElements_per_direction_store,
        nValues_per_row=nValues_per_row_store,
        num_procs=num_procs_store,
        file_extension="dat")
#-----------------------------------------------------
from generate_philip_input_files_from_velocity_field_for_spectra_fix import generate_philip_input_files
def batch_generate_philip_input_files(
    file_path=[],
    n_different_files_in_path=[],
    file_prefix=[],
    poly_degree=[],
    nElements_per_direction=[],
    nValues_per_row=[],
    num_procs=[],# not needed but could overwrite to be fixed as 8
    file_extension="dat"):
    #-----------------------------------------------------
    # Safeguard for when empty args are passed
    #-----------------------------------------------------
    if(file_path==[] or file_prefix==[] or poly_degree==[] or nElements_per_direction==[] or num_procs==[]):
        print("batch_generate_philip_input_files: an empty essential argument was passed")
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
        print("Generating for path: %s" % file_path[i])
        for j in range(0,n_prefixes):
            #-----------------------------------------------------
            # Get DOF variables
            #-----------------------------------------------------
            nElements,nQuadPoints_per_element,nQuadPoints,nDOF,reduced_nQuadPoints,reduced_nDOF = get_DOF_vars(nElements_per_direction[i],poly_degree[i])

            input_vel_field_filename_=file_path[i]+"velocity_vorticity-"+str(j)+"_reordered_gll_nodes.dat"
            print(" - Generating for prefix: %s" % j)
            
            generate_philip_input_files(
                nElements_per_direction[i],
                nQuadPoints_per_element,
                6,# nValues_per_row
                nDOF,
                8,# number of processors to output to
                output_dir=file_path[i],
                input_vel_field_filename=input_vel_field_filename_,
                prefix_string=str(j))
            print(" - done.")
        print("done.")
            
    return
#-----------------------------------------------------