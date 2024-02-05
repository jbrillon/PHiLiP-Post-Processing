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