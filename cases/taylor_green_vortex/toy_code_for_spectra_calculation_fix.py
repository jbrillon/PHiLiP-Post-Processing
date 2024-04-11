#-----------------------------------------------------
# Import public libraries
import numpy as np # NumPy: contains basic numerical routines
#-----------------------------------------------------
import os;CURRENT_PATH = os.path.split(os.path.realpath(__file__))[0]+"/";
import sys
# load tools
sys.path.append(CURRENT_PATH+"../../src/tools");
from batch_assemble_mpi_flow_field_files_reorder_generate_spectra import batch_convert_velocity_field_at_equidistant_nodes_to_gLL_nodes_from_txt
from generate_spectra_files import generate_spectra_file_from_flow_field_file
# load submodules
sys.path.append(CURRENT_PATH+"../../submodules/quickplotlib/lib"); import quickplotlib as qp
# for testing quickplotlib changes; uncomment the lines below and comment the line above
# from sys import platform
# if platform == "linux" or platform == "linux2":
#     # linux
#     sys.path.append("/home/julien/Codes/quickplotlib/lib"); import quickplotlib as qp # uncomment if testing quickplotlib changes
# elif platform == "darwin":
#     # OS X
#     sys.path.append("/Users/Julien/Python/quickplotlib/lib"); import quickplotlib as qp # uncomment if testing quickplotlib changes

#-----------------------------------------------------
from sys import platform
if platform == "linux" or platform == "linux2":
    # linux
    filesystem="/media/julien/Samsung_T5/"
elif platform == "darwin":
    # OS X
    filesystem="/Volumes/Samsung_T5/"

# testing
file = sys.argv[1] # to pass file from terminal
print("Specified input file for paths: %s" % file)
batch_convert_velocity_field_at_equidistant_nodes_to_gLL_nodes_from_txt(file)
exit()
# filepath=filesystem+"NarvalFiles/2023_JCP/robustness/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs048_p5_procs64/flow_field_files/"
filepath=filesystem+"NarvalFiles/2023_JCP/robustness/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs024_p5_procs16/flow_field_files/"
input_vel_field_filename_="velocity_vorticity-0_reordered_gll_nodes.dat"
nElements_per_direction = 4
nQuadPoints_per_element = 6
nValues_per_row = 6
nDOF = 13824
num_procs = 1

import sys; sys.path.append(CURRENT_PATH+"../../../DHIT-Flow-Setup/v2");
from reorder_for_philip_prep_for_spectra_fix import generate_philip_input_files
generate_philip_input_files(
    nElements_per_direction,
    nQuadPoints_per_element,
    nValues_per_row,
    nDOF,
    num_procs,
    output_dir=filepath,
    input_vel_field_filename=input_vel_field_filename_)
