import numpy as np

def write_vtk_file_uniform_cube(
    unique_coordinates,
    velocities,
    filename="solution.vtk",
    additional_scalars_arrays=[],
    additional_scalars_names=[]):

    # =============================
    #       WRITE THE VTK FILE
    # =============================
    # examples: https://visit-sphinx-github-user-manual.readthedocs.io/en/task-allen-vtk9_master_ospray/data_into_visit/VTKFormat.html
    # official documentation: https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
    nDOF_total = int(unique_coordinates.shape[0])
    nDOF_per_dim = int(round(nDOF_total**(1.0/3.0)))
    print("Writing vtk file: %s ... " % filename)
    dataType = "double"
    file = open(filename,"w") # change your vtk file name
    file.write("# vtk DataFile Version 2.0\n")
    file.write("Cube example\n")
    file.write("ASCII\n")
    file.write("DATASET STRUCTURED_GRID\n")
    file.write("\n")
    file.write("DIMENSIONS %i %i %i\n" % (nDOF_per_dim,nDOF_per_dim,nDOF_per_dim))
    file.write("POINTS %i %s\n" % (nDOF_total,dataType))
    for i in range(0,nDOF_total):
        wstr = "%1.15f %1.15f %1.15f\n" % (unique_coordinates[i,0],unique_coordinates[i,1],unique_coordinates[i,2])
        file.write(wstr)
    file.write("\n")
    file.write("POINT_DATA %i\n" % (nDOF_total))
    # write 3 velocity components
    file.write("SCALARS u %s 1\n" % dataType)
    file.write("LOOKUP_TABLE default\n")
    for i in range(0,nDOF_total):
        wstr = "%1.15f\n" % (velocities[i,0])
        file.write(wstr)
    file.write("SCALARS v %s 1\n" % dataType)
    file.write("LOOKUP_TABLE default\n")
    for i in range(0,nDOF_total):
        wstr = "%1.15f\n" % (velocities[i,1])
        file.write(wstr)
    file.write("SCALARS w %s 1\n" % dataType)
    file.write("LOOKUP_TABLE default\n")
    for i in range(0,nDOF_total):
        wstr = "%1.15f\n" % (velocities[i,2])
        file.write(wstr)
    # write additional scalars
    if(additional_scalars_arrays!=[]):
        if(hasattr(additional_scalars_arrays[0],'__len__')):
            n_additional_scalars = int(len(additional_scalars_arrays))
            for j in range(0,n_additional_scalars):
                values = additional_scalars_arrays[j]
                name = additional_scalars_names[j]
                file.write("SCALARS %s %s 1\n" % (name,dataType))
                file.write("LOOKUP_TABLE default\n")
                for i in range(0,nDOF_total):
                    wstr = "%1.15f\n" % (values[i])
                    file.write(wstr)
    # write velocity vector
    file.write("VECTORS velocity %s\n" % dataType)
    for i in range(0,nDOF_total):
        wstr = "%1.15f %1.15f %1.15f\n" % (velocities[i,0],velocities[i,1],velocities[i,2])
        file.write(wstr)
    file.close()
    print("done.")
