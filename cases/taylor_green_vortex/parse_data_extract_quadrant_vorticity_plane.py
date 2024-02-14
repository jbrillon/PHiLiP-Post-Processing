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
# paths=(
# "/home/julien/NarvalFiles/2023_JCP/filtered_dns_viscous_tgv/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p7_procs1024/solution_files/",\
# "/home/julien/NarvalFiles/2023_JCP/verification/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p3_procs1024/solution_files/",\
# "/home/julien/NarvalFiles/2023_JCP/verification/viscous_TGV_ILES_NSFR_cDG_IR_2PF-Roe_GL_OI-0_dofs0256_p3_procs1024/solution_files/",\
# "/home/julien/NarvalFiles/2023_JCP/verification/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-4_dofs0256_p3_CFL-0.15_procs1024/solution_files/",\
# )
# files_per_path=[10,10,7,7]
paths=(
# filesystem+"NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/solution_files/",\
filesystem+"NarvalFiles/2023_JCP/verification/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GLL_OI-0_dofs0256_p7_procs1024/solution_files/",\
)
files_per_path=[10]
prefix="vorticity_mag_slice_x0_plane_t_"
file_extension="txt"
scalar_field_name = "vorticity_magnitude"
n_paths=len(paths)

# x-coordinate of slice plane 
x_plane = 0.0

# Quadrant bounding coordinates
y_min = 0.0
y_max = np.pi+0.1
z_min = -(np.pi+0.1)
z_max = 0.0

# Generate quadrant files
for i in range(0,n_paths):
    for j in range(0,files_per_path[i]):
        filename_without_extension = paths[i]+prefix+str(j)
        eqfps.extract_section_from_paraview_slice(x_plane,y_min,y_max,z_min,z_max,scalar_field_name,filename_without_extension,file_extension)

