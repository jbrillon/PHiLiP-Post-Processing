# import matplotlib.tri as tri
import numpy as np
#-------------------------------------------------------------
import sys
from sys import platform
if platform == "linux" or platform == "linux2":
    # linux
    filesystem="/media/julien/Samsung_T5/"
elif platform == "darwin":
    # OS X
    filesystem="/Volumes/Samsung_T5/"
#-------------------------------------------------------------
# load tools
import os;CURRENT_PATH = os.path.split(os.path.realpath(__file__))[0]+"/";
sys.path.append(CURRENT_PATH+"../../src/tools");
import extract_quadrant_from_paraview_slice as eqfps
#-------------------------------------------------------------
paths=(
# filesystem+"NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/flow_field_files/",\
filesystem+"NarvalFiles/2023_JCP/spectra_fix/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/flow_field_files/",\
)
files_per_path=[2,2]#[10]
file_extension="dat"
scalar_field_name = "vorticity_magnitude"
n_paths=len(paths)

# x-coordinate of slice plane 
x_plane = np.float64(0.0000000000000000e+00)

# Quadrant bounding coordinates
y_min = -3.1415926535897931e+00
y_max = 3.1415926535897931e+00
z_min = -3.1415926535897931e+00
z_max = 3.1415926535897931e+00

# Generate quadrant files
for i in range(0,n_paths):
    for j in range(1,files_per_path[i]):
        prefix = "velocity_vorticity-%i_reordered" % j
        filename_without_extension = paths[i]+prefix
        eqfps.extract_section_from_flow_field_file(True,x_plane,y_min,y_max,z_min,z_max,scalar_field_name,filename_without_extension,file_extension)

