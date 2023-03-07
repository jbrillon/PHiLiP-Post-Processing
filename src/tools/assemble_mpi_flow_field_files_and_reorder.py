import numpy as np
import sys; sys.path.append("../../src/tools");
from get_DOF_vars import get_DOF_vars
from assemble_mpi_files import assemble_mpi_files

# assemble MPI flow field files and reorder
def assemble_mpi_flow_field_files_and_reorder(
    file_path_and_prefix,
    file_extension,
    num_procs,
    nValues_per_row,
    nElements_per_direction,
    poly_degree):
    
    #-----------------------------------------------------
    ''' Description:
        This function writes the file "{file_path_and_prefix}_reordered.{file_extension}"
    '''
    #-----------------------------------------------------
    
    #-----------------------------------------------------
    # Assemble the MPI files
    #-----------------------------------------------------
    assemble_mpi_files(file_path_and_prefix,file_extension,num_procs) # Note: This writes the file "{file_path_and_prefix}.{file_extension}"
    #-----------------------------------------------------

    #-----------------------------------------------------
    # Get DOF variables
    #-----------------------------------------------------
    nElements,nQuadPoints_per_element,nQuadPoints,nDOF,reduced_nQuadPoints,reduced_nDOF = get_DOF_vars(nElements_per_direction,poly_degree)
    #-----------------------------------------------------
    # Setup the loop bounds
    #-----------------------------------------------------
    nLoops = 4
    loop_bounds = np.ones(nLoops,dtype=np.int64)
    if(nElements_per_direction>=4):
        loop_bounds[0] = 2
    if(nElements_per_direction>=8):
        loop_bounds[1] = 2
    if(nElements_per_direction>=16):
        loop_bounds[2] = 2
    if(nElements_per_direction>=32):
        loop_bounds[3] = 2
    # if(nElements_per_direction>=64):
    #     loop_bounds[4] = 2
    # if(nElements_per_direction>=128):
    #     loop_bounds[5] = 2
    # if(nElements_per_direction>=256):
    #     loop_bounds[6] = 2
    # if(nElements_per_direction>=512):
    #     loop_bounds[7] = 2
    #-----------------------------------------------------

    #-------------------------------------------------------------
    # Store the flow field + coordinates
    #-------------------------------------------------------------
    stored_data = np.zeros((nElements_per_direction,nElements_per_direction,nElements_per_direction,nQuadPoints_per_element,nQuadPoints_per_element,nQuadPoints_per_element,1,nValues_per_row),dtype=np.float64)

    velocity_file_from_philip = file_path_and_prefix+"."+file_extension
    fin = open(velocity_file_from_philip,"r")

    # First line: Number of DOFs
    nDOF_expected = int(fin.readline())
    if(nDOF!=nDOF_expected):
        print("Error: nDOF does not match expected nDOF from file %s, check function inputs.",velocity_file_from_philip)
        print("Aborting...")
        exit()

    ''' must add more nested for loops for higher
        number of elements per direction
        currently can handle up to 32 (i.e. 2,4,8,16,32)
    '''

    iproc = 0
    iDOF_per_proc = 0
    start_new_file=True

    ez_L_base_base_base = 0
    for z_base_base_base in range(0,loop_bounds[3]):
        ey_L_base_base_base = 0
        for y_base_base_base in range(0,loop_bounds[3]):
            ex_L_base_base_base = 0
            for x_base_base_base in range(0,loop_bounds[3]):
                ez_L_base_base = ez_L_base_base_base
                for z_base_base in range(0,loop_bounds[2]):
                    ey_L_base_base = ey_L_base_base_base
                    for y_base_base in range(0,loop_bounds[2]):
                        ex_L_base_base = ex_L_base_base_base
                        for x_base_base in range(0,loop_bounds[2]):
                            ez_L_base = ez_L_base_base
                            for z_base in range(0,loop_bounds[1]):
                                ey_L_base = ey_L_base_base
                                for y_base in range(0,loop_bounds[1]):
                                    ex_L_base = ex_L_base_base
                                    for x_base in range(0,loop_bounds[1]):
                                        # algorithm for a cube with 64 (4^3) elements:
                                        ez_L = ez_L_base
                                        for cz in range(0,loop_bounds[0]):
                                            ez_R = ez_L + 1
                                            ey_L = ey_L_base
                                            for cy in range(0,loop_bounds[0]):
                                                ey_R = ey_L + 1
                                                ex_L = ex_L_base
                                                for cx in range(0,loop_bounds[0]):
                                                    ex_R = ex_L + 1
                                                    for ez in range(ez_L,ez_R+1):
                                                        for ey in range(ey_L,ey_R+1):
                                                            for ex in range(ex_L,ex_R+1):
                                                                for qz in range(0,nQuadPoints_per_element):
                                                                    for qy in range(0,nQuadPoints_per_element):
                                                                        for qx in range(0,nQuadPoints_per_element):
                                                                            row_string = fin.readline()
                                                                            row_data = np.fromstring(row_string, dtype=np.float64, sep=' ')
                                                                            for iValue in range(0,nValues_per_row):
                                                                                stored_data[ez,ey,ex,qz,qy,qx,0,iValue] = row_data[iValue] # modify to read in vorticity
                                                    ex_L += 2
                                                ey_L += 2
                                            ez_L += 2
                                        ex_L_base = ex_L
                                    ey_L_base = ey_L
                                ez_L_base = ez_L
                            ex_L_base_base = ex_L_base
                        ey_L_base_base = ey_L_base
                    ez_L_base_base = ez_L_base
                ex_L_base_base_base = ex_L_base_base
            ey_L_base_base_base = ey_L_base_base
        ez_L_base_base_base = ez_L_base_base
    #-------------------------------------------------------------
    # Write the reordered flow field file
    #-------------------------------------------------------------
    reordered_flow_field_file = file_path_and_prefix+"_reordered"+"."+file_extension
    file = open(reordered_flow_field_file,"w")

    # Write number of degrees of freedom
    wstr = "%i\n" % nDOF
    file.write(wstr)

    for ez in range(0,nElements_per_direction):
        for qz in range(0,nQuadPoints_per_element):
            for ey in range(0,nElements_per_direction):
                for qy in range(0,nQuadPoints_per_element):
                    for ex in range(0,nElements_per_direction):
                        for qx in range(0,nQuadPoints_per_element):
                            wstr = "%18.16e" % stored_data[ez,ey,ex,qz,qy,qx,0,0]
                            file.write(wstr)
                            for iValue in range(1,nValues_per_row):
                                wstr = " %18.16e" % stored_data[ez,ey,ex,qz,qy,qx,0,iValue]
                                file.write(wstr)
                            wstr = "\n"
                            file.write(wstr)
    file.close()
    return
