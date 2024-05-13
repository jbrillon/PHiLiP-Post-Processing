import numpy as np
import pandas as pd
#-------------------------------------------------------------
def extract_section_from_paraview_slice(
    x_plane,
    y_min,
    y_max,
    z_min,
    z_max,
    scalar_field_name,
    filename_without_extension,file_extension):
    input_filename = filename_without_extension+"."+file_extension
    # Load values
    df = pd.read_csv(input_filename,sep=",",header=0)
    scalar = df[scalar_field_name].to_numpy()
    x = df["Points:0"].to_numpy()
    y = df["Points:1"].to_numpy()
    z = df["Points:2"].to_numpy()

    #-------------------------------------------------------------
    # Write the quadrant file
    #-------------------------------------------------------------
    quadrant_filename = filename_without_extension+"_quadrant."+file_extension
    file = open(quadrant_filename,"w")
    print("Writting to file: %s ..." % quadrant_filename)

    nDOFs = len(x)

    for i in range(0,nDOFs):
        if(x[i] == x_plane):
            if(((y[i]>=y_min) and (y[i]<=y_max)) and ((z[i]>=z_min) and (z[i]<=z_max))):
                wstr = "%1.16e %1.16e %1.16e\n" % (y[i],z[i],scalar[i])
                file.write(wstr)
    file.close()
    print("done.")
    return
#-------------------------------------------------------------
def extract_section_from_flow_field_file(
    is_interface,
    x_plane,
    y_min,
    y_max,
    z_min,
    z_max,
    scalar_field_name,
    filename_without_extension,
    file_extension):
    input_filename = filename_without_extension+"."+file_extension
    # Load values
    column_index_of_scalar = 6 # default
    if(scalar_field_name=="vorticity_magnitude"):
        column_index_of_scalar = 6
    elif(scalar_field_name=="x-velocity"):
        column_index_of_scalar = 3
    elif(scalar_field_name=="x-velocity"):
        column_index_of_scalar = 4
    elif(scalar_field_name=="x-velocity"):
        column_index_of_scalar = 5
    else:
        print("ERROR: Invalid scalar_field_name passed to extract_section_from_flow_field_file().")
        print("Aborting...")
        exit()
    # load data
    x,y,z,scalar = np.loadtxt(input_filename,skiprows=1,usecols=(0,1,2,column_index_of_scalar),unpack=True)
    nDOFs = len(x)
    #-------------------------------------------------------------
    # Write the quadrant file
    #-------------------------------------------------------------
    quadrant_filename = filename_without_extension+"_slice."+file_extension
    file = open(quadrant_filename,"w")
    print("Scanning %i DOFs and writting to file: %s ..." % (nDOFs,quadrant_filename))

    for i in range(0,nDOFs):
        if(x[i] == x_plane):
            # if(((y[i]>=y_min) and (y[i]<=y_max)) and ((z[i]>=z_min) and (z[i]<=z_max))):
            wstr = "%1.16e %1.16e %1.16e\n" % (y[i],z[i],scalar[i])
            file.write(wstr)
    file.close()
    print("done.")

    if(is_interface):
        # average value at face
        quadrant_filename = filename_without_extension+"_slice."+file_extension
        y,z,scalar = np.loadtxt(quadrant_filename,skiprows=0,unpack=True)
        nDOFs = len(y)
        
        nDOFs_unique = int(nDOFs/2)
        output_quadrant_filename = filename_without_extension+"_slice_interface_averaged."+file_extension
        file = open(output_quadrant_filename,"w")
        
        print("x_plane is an interface! Scanning %i DOFs, averaging the interface values, and writting to file: %s ..." % (nDOFs,output_quadrant_filename))

        index_left = 0
        index_right = index_left+1
        for dummy in range(0,nDOFs_unique):
            # check if ok
            if(y[index_left]!=y[index_right]):
                print("ERROR: y[index_left]!=y[index_right] for index_left=%i and index_right=%i" % (index_left,index_right))
                print("Aborting...")
                exit()
            # do stuff
            scalar_average = 0.5*(scalar[index_left]+scalar[index_right])
            wstr = "%1.16e %1.16e %1.16e\n" % (y[index_left],z[index_left],scalar_average)
            file.write(wstr)

            # update index
            index_left += 2 # go to next pair
            index_right = index_left+1
            
        file.close()
        print("done.")
        print(index_left)
        print(index_right)

    return
#-------------------------------------------------------------