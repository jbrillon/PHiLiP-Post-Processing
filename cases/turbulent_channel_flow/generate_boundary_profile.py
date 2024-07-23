#-----------------------------------------------------
# Import public libraries
import numpy as np # NumPy: contains basic numerical routines
#-----------------------------------------------------
# Import personal libraries
# from finite_difference_library import first_derivative, fd_non_uniform_grid
import os;CURRENT_PATH = os.path.split(os.path.realpath(__file__))[0]+"/";
import sys
# load tools
sys.path.append(CURRENT_PATH+"../../src/tools");
from assemble_mpi_files import assemble_mpi_files
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


#-----------------------------------------------------
def generate_boundary_layer_profile_file_from_flow_field_file(
    file_without_extension,
    file_extension,
    number_of_elements_in_y_direction,
    poly_degree,
    n_skiprows=1):

    print("Generating boundary layer profile file from flow field file...")
    file = file_without_extension+"."+file_extension
    print(" - Loading file for just u+ vs y+ data: %s" % file)
    data = np.loadtxt(file,skiprows=n_skiprows,usecols=(1,3,7,8),dtype=np.float64)
    print(" - done.")
    # coordinates = [data[:,0],data[:,1],data[:,2]] # x,y,z
    print(" - Loading file for entire velocity field: %s" % file)
    vel_data = np.loadtxt(file,skiprows=n_skiprows,usecols=(3,4,5),dtype=np.float64)
    print(" - done.")
    velocity_field = [vel_data[:,0],vel_data[:,1],vel_data[:,2]] # u,v,w
    print(" - Computing the fluctuating velocity field...")
    velocity_field = get_fluctuating_velocity_field(velocity_field)
    print(" - done.")

    file_out_without_extension = file_without_extension+"_boundary_layer_profile"

    # NOTE FOR REYNOLDS STRESSES; call get_fluctuating_velocity_field(velocity_field)
    number_of_degrees_of_freedom = np.size(data[:,0])
    number_of_unique_points_per_direction = number_of_elements_in_y_direction*(poly_degree+1) - (number_of_elements_in_y_direction-1) # total DOFs minus number of interfaces
    # (1) compute the expected unique y-stations
    unique_y_stations = np.linspace(-1.0,1.0,number_of_unique_points_per_direction)
    number_of_y_stations = np.size(unique_y_stations) # same as number_of_unique_points_per_direction

    print(" - Averaging the x-velocity at each y-station...")
    # (2) sum all x-velocities and kinematic viscosity at each y-station
    # - x-velocity
    summation_of_xvelocity_at_each_y_station = np.zeros(number_of_y_stations)
    number_of_summations_of_xvelocity_at_each_y_station = np.zeros(number_of_y_stations)
    average_xvelocity_at_each_y_station = np.zeros(number_of_y_stations)
    # - kinematic viscosity
    summation_of_kinematic_viscosity_at_each_y_station = np.zeros(number_of_y_stations)
    number_of_summations_of_kinematic_viscosity_at_each_y_station = np.zeros(number_of_y_stations)
    average_kinematic_viscosity_at_each_y_station = np.zeros(number_of_y_stations)
    # - x-velocity fluctuation RMS
    summation_of_xvelocity_fluctuation_rms_at_each_y_station = np.zeros(number_of_y_stations)
    number_of_summations_of_xvelocity_fluctuation_rms_at_each_y_station = np.zeros(number_of_y_stations)
    average_xvelocity_fluctuation_rms_at_each_y_station = np.zeros(number_of_y_stations)
    # - y-velocity fluctuation RMS
    summation_of_yvelocity_fluctuation_rms_at_each_y_station = np.zeros(number_of_y_stations)
    number_of_summations_of_yvelocity_fluctuation_rms_at_each_y_station = np.zeros(number_of_y_stations)
    average_yvelocity_fluctuation_rms_at_each_y_station = np.zeros(number_of_y_stations)
    # - z-velocity fluctuation RMS
    summation_of_zvelocity_fluctuation_rms_at_each_y_station = np.zeros(number_of_y_stations)
    number_of_summations_of_zvelocity_fluctuation_rms_at_each_y_station = np.zeros(number_of_y_stations)
    average_zvelocity_fluctuation_rms_at_each_y_station = np.zeros(number_of_y_stations)
    for i in range(0,number_of_degrees_of_freedom):
        # check which y-station
        y_value = 1.0*data[i,0]
        index_of_y_station = np.where(np.abs(unique_y_stations-y_value)<1.0e-8)[0][0]
        # sum the x-velocity
        velocity_x_direction_value = 1.0*data[i,1]
        summation_of_xvelocity_at_each_y_station[index_of_y_station] += velocity_x_direction_value
        number_of_summations_of_xvelocity_at_each_y_station[index_of_y_station] += 1.0
        # sum the kinematic viscosity
        density_value = 1.0*data[i,2]
        dynamic_viscosity_value = 1.0*data[i,3]
        kinematic_viscosity_value = dynamic_viscosity_value/density_value
        summation_of_kinematic_viscosity_at_each_y_station[index_of_y_station] += kinematic_viscosity_value
        number_of_summations_of_kinematic_viscosity_at_each_y_station[index_of_y_station] += 1.0
        # sum the x-velocity fluctuation RMS
        velocity_x_direction_fluctuation_value = 1.0*velocity_field[0][i]
        summation_of_xvelocity_fluctuation_rms_at_each_y_station = velocity_x_direction_fluctuation_value*velocity_x_direction_fluctuation_value
        number_of_summations_of_xvelocity_fluctuation_rms_at_each_y_station[index_of_y_station] += 1.0
        # sum the y-velocity fluctuation RMS
        velocity_y_direction_fluctuation_value = 1.0*velocity_field[1][i]
        summation_of_yvelocity_fluctuation_rms_at_each_y_station = velocity_y_direction_fluctuation_value*velocity_y_direction_fluctuation_value
        number_of_summations_of_yvelocity_fluctuation_rms_at_each_y_station[index_of_y_station] += 1.0
        # sum the z-velocity fluctuation RMS
        velocity_z_direction_fluctuation_value = 1.0*velocity_field[2][i]
        summation_of_zvelocity_fluctuation_rms_at_each_y_station = velocity_z_direction_fluctuation_value*velocity_z_direction_fluctuation_value
        number_of_summations_of_zvelocity_fluctuation_rms_at_each_y_station[index_of_y_station] += 1.0
    # (3) compute the average
    for i in range(0,number_of_y_stations):
        average_xvelocity_at_each_y_station[i] = summation_of_xvelocity_at_each_y_station[i]/number_of_summations_of_xvelocity_at_each_y_station[i]
        average_kinematic_viscosity_at_each_y_station[i] = summation_of_kinematic_viscosity_at_each_y_station[i]/number_of_summations_of_kinematic_viscosity_at_each_y_station[i]
        average_xvelocity_fluctuation_rms_at_each_y_station[i] = np.sqrt(summation_of_xvelocity_fluctuation_rms_at_each_y_station[i]/number_of_summations_of_xvelocity_fluctuation_rms_at_each_y_station[i])
        average_yvelocity_fluctuation_rms_at_each_y_station[i] = np.sqrt(summation_of_yvelocity_fluctuation_rms_at_each_y_station[i]/number_of_summations_of_yvelocity_fluctuation_rms_at_each_y_station[i])
        average_zvelocity_fluctuation_rms_at_each_y_station[i] = np.sqrt(summation_of_zvelocity_fluctuation_rms_at_each_y_station[i]/number_of_summations_of_zvelocity_fluctuation_rms_at_each_y_station[i])
    print(" - done.")
    # (4) output file for the boundary layer profile    
    file_out = file_out_without_extension+"."+file_extension
    print(" - Writing file: %s" % file_out)
    np.savetxt(file_out,np.transpose(np.array([unique_y_stations, average_xvelocity_at_each_y_station, average_kinematic_viscosity_at_each_y_station, average_xvelocity_fluctuation_rms_at_each_y_station, average_yvelocity_fluctuation_rms_at_each_y_station, average_zvelocity_fluctuation_rms_at_each_y_station])))
    print(" - done.")

    print("done.")
    return
#-----------------------------------------------------
def batch_assemble_mpi_flow_field_files_generate_boundary_layer_profile(
    file_path=[],
    n_different_files_in_path=[],
    file_prefix=[],
    poly_degree=[],
    nElements_x_direction=[],
    nElements_y_direction=[],
    nElements_z_direction=[],
    nValues_per_row=[],
    num_procs=[],
    file_extension="dat"):
    
    #-----------------------------------------------------
    # Safeguard for when empty args are passed
    #-----------------------------------------------------
    if(file_path==[] or file_prefix==[] or poly_degree==[] or nElements_x_direction==[] or nElements_y_direction==[] or nElements_z_direction==[] or num_procs==[]):
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
            prefix = file_prefix[i][j]
            
            # get file path and prefix
            file_path_and_prefix = file_path[i] + prefix

            #-----------------------------------------------------
            # Assemble the MPI files
            #-----------------------------------------------------
            assemble_mpi_files(file_path_and_prefix,file_extension,num_procs[i]) # Note: This writes the file "{file_path_and_prefix}.{file_extension}"
            #-----------------------------------------------------
            
            # generate the boundary layer profile file
            velocity_file_for_boundary_layer_profile_without_extension = file_path_and_prefix
            generate_boundary_layer_profile_file_from_flow_field_file(
                velocity_file_for_boundary_layer_profile_without_extension,
                file_extension,
                nElements_y_direction[i],
                poly_degree[i],
                n_skiprows=1)
    return
#-----------------------------------------------------

# Inputs
# file_path=["/Users/Julien/NarvalFiles/viscous_TCF_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_Re5200_p4_20x10x10_turbulent_initialization/flow_field_files/"]
file_path=["/home/julien/Codes/dummy_dir_test_channel_flow/flow_field_files/"]
n_different_files_in_path=[1]

# generate the standard prefix files
file_prefix=[]
prefix_base = "velocity_vorticity-"
n_files_per_dir = 1
standard_prefix_files = []
for i in range(0,n_files_per_dir):
    standard_prefix_files.append(prefix_base+str(i))
file_prefix.append(standard_prefix_files)
# poly_degree_original = 4
# formula for with subdivisions: this->output_velocity_number_of_subvisions*(dg->max_degree+1)-1
poly_degree=[14]#4,9,14 
nElements_x_direction=[20]
nElements_y_direction=[10]
nElements_z_direction=[10]
nValues_per_row=[7]
num_procs=[8]
file_extension="dat"
# Call function
batch_assemble_mpi_flow_field_files_generate_boundary_layer_profile(
    file_path=file_path,
    n_different_files_in_path=n_different_files_in_path,
    file_prefix=file_prefix,
    poly_degree=poly_degree,
    nElements_x_direction=nElements_x_direction,
    nElements_y_direction=nElements_y_direction,
    nElements_z_direction=nElements_z_direction,
    nValues_per_row=nValues_per_row,
    num_procs=num_procs,
    file_extension="dat")