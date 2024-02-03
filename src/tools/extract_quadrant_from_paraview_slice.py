import numpy as np
import pandas as pd
# import matplotlib.tri as tri
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
path = filesystem+"NarvalFiles/2023_JCP/verification/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p3_procs1024/solution_files/"
input_filename = path+"vorticity_plane_from_paraview_x_zero.txt"

# Load values
df = pd.read_csv(input_filename,sep=",",header=0)
scalar = df["vorticity_magnitude"].to_numpy()
x = df["Points:0"].to_numpy()
y = df["Points:1"].to_numpy()
z = df["Points:2"].to_numpy()

# x-coordinate of slice plane 
x_plane = 0.0

# Quadrant bounding coordinates
y_min = 0.0
y_max = np.pi+0.1
z_min = -(np.pi+0.1)
z_max = 0.0

#-------------------------------------------------------------
# Write the quadrant file
#-------------------------------------------------------------
# vorticity_magnitude_plane_file = file_path_and_prefix+"_reordered"+"."+file_extension
quadrant_filename = "extracted_quadrant_from_slice.dat"
quadrant_filename_full = path+quadrant_filename
file = open(quadrant_filename_full,"w")
print("Writting to file: %s ..." % quadrant_filename_full)

nDOFs = len(x)

for i in range(0,nDOFs):
    if(x[i] == x_plane):
        if(((y[i]>=y_min) and (y[i]<=y_max)) and ((z[i]>=z_min) and (z[i]<=z_max))):
            wstr = "%1.16e %1.16e %1.16e\n" % (y[i],z[i],scalar[i])
            file.write(wstr)
file.close()
print("done.")
