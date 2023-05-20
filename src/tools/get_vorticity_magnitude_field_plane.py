import numpy as np;

# filename="/Users/Julien/NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/flow_field_files/velocity_vorticity-0-00000.dat"
# nDOFs=1728

# filename="/Users/Julien/NarvalFiles/2023_JCP/flux_nodes/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs096_p5_procs512/flow_field_files/velocity_vorticity-0_reordered.dat"
filename="/Users/Julien/NarvalFiles/2023_JCP/verification/viscous_TGV_ILES_NSFR_cDG_IR_2PF_GL_OI-0_dofs0256_p3_procs1024/flow_field_files/velocity_vorticity-0_reordered.dat"
nDOFs=np.loadtxt(filename,max_rows=1,dtype=np.int64)
data=np.loadtxt(filename,skiprows=1,usecols=(0,1,2,6),dtype=np.float64)

# x-coordinate of plane at which we extract the values; reference from the HOPW-5 viscous TGV case
x_plane=-3.1415926535897931e+00

y_plane_store=[]
z_plane_store=[]
vorticity_magnitude_plane_store=[]

#-------------------------------------------------------------
# Write the reordered flow field file
#-------------------------------------------------------------
# vorticity_magnitude_plane_file = file_path_and_prefix+"_reordered"+"."+file_extension
vorticity_magnitude_plane_file = "vorticity_x_plane_DNS.dat"
file = open(vorticity_magnitude_plane_file,"w")
print("Writting to file: %s ..." % vorticity_magnitude_plane_file)

# # Write number of degrees of freedom
# wstr = "%i\n" % nDOF
# file.write(wstr)
for i in range(0,nDOFs):
    x_coord = data[i,0]
    if(x_coord==x_plane):
        wstr = "%18.16e %18.16e %18.16e" % (data[i,1],data[i,2],data[i,3])
        file.write(wstr)
        wstr = "\n"
        file.write(wstr)
file.close()
print("done.")
