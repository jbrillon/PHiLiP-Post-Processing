import numpy as np;

from sys import platform
if platform == "linux" or platform == "linux2":
    # linux
    filesystem="/media/julien/Samsung_T5/"
elif platform == "darwin":
    # OS X
    filesystem="/Volumes/Samsung_T5/"

# flow_field_filename="/Users/Julien/NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/flow_field_files/velocity_vorticity-0-00000.dat"
# nDOFs=1728

flow_field_filename=filesystem+"NarvalFiles/2023_JCP/verification/viscous_TGV_ILES_std_strong_DG_Roe_GL_OI-4_dofs0256_p3_CFL-0.15_procs1024/flow_field_files/velocity_vorticity-0_reordered.dat"
# flow_field_filename=filesystem+"NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/flow_field_files/velocity_vorticity-0_reordered.dat"
# flow_field_filename=filesystem+"NarvalFiles/2023_JCP/verification/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p3_procs1024/flow_field_files/velocity_vorticity-1_reordered.dat"
# nDOFs=np.loadtxt(filename,max_rows=1,dtype=np.int64)
# data=np.loadtxt(filename,skiprows=1,usecols=(0,1,2,6),dtype=np.float64)

# x-coordinate of plane at which we extract the values; reference from the HOPW-5 viscous TGV case
x_plane=-3.1415926535897931e+00
y_min=0.0000000000000000e+00
y_max=1.5707963267948966e+00
z_min=1.5707963267948966e+00
z_max=3.1415926535897931e+00

# y_plane_store=[]
# z_plane_store=[]
# vorticity_magnitude_plane_store=[]

#-------------------------------------------------------------
# Write the reordered flow field file
#-------------------------------------------------------------
# vorticity_magnitude_plane_file = file_path_and_prefix+"_reordered"+"."+file_extension
vorticity_magnitude_plane_file = "vorticity_x_plane_strong_DG_DNS-0.dat"
file = open(vorticity_magnitude_plane_file,"w")
print("Writting to file: %s ..." % vorticity_magnitude_plane_file)


fin = open(flow_field_filename,"r")
# First line: Number of DOFs
nDOFs = int(fin.readline())

# # Write number of degrees of freedom
# wstr = "%i\n" % nDOF
# file.write(wstr)
for i in range(0,nDOFs):
    row_string = fin.readline()
    data = np.fromstring(row_string, dtype=np.float64, sep=' ')
    x_coord = data[0]
    if(x_coord==x_plane):
        y_coord=data[1]
        z_coord=data[2]
        vorticity_magnitude = data[6]
        # if(((y_coord>=y_min) and (y_coord<=y_max)) and ((z_coord>=z_min) and (z_coord<=z_max)))
        # 
        wstr = "%18.16e %18.16e %18.16e" % (y_coord,z_coord,vorticity_magnitude)
        file.write(wstr)
        wstr = "\n"
        file.write(wstr)
fin.close()
file.close()
print("done.")
