# get padded mpi rank string
def get_padded_mpi_rank_string(mpi_rank):
    padding_length = 5
    mpi_rank_string = '%i' % mpi_rank
    mpi_rank_string = mpi_rank_string.zfill(padding_length)
    return mpi_rank_string

# assemble mpi files
def assemble_mpi_files(file_path_and_prefix,file_extension,num_procs):
    
    print("Assembling %i files with prefix: %s" % (num_procs,file_path_and_prefix))
    
    ext = "."+file_extension
    filename=file_path_and_prefix+ext

    print("Writting to file: %s ..." % filename)
    
    fout = open(filename, "w")
    for i in range(0,num_procs):
        mpi_rank_string = get_padded_mpi_rank_string(i)
        tempfile = file_path_and_prefix+"-"+mpi_rank_string+ext
        fin = open(tempfile,"r")
        for line in fin:
            fout.write(line)
        fin.close()
    fout.close()
    print("done.")
    return