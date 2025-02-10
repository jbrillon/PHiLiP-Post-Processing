restartIndices=({30..49})
# restartIndices=("49")
# INPUTS:
oversampling_factor="1"
restart_file_directory="restart_files_cPlus_C1/"
target_flow_field_subdir="cPlus_C1/"
# constants:
flow_field_directory="/home/julien/Codes/dummy_dir_test_channel_flow/flow_field_files/"
file_we_copy="${flow_field_directory}velocity_vorticity-0_boundary_layer_profile.dat"
for i in ${!restartIndices[@]}; do
    restart_index=${restartIndices[${i}]}
    filename="restart-000${restart_index}.prm"
    echo "Replacing text in file ${filename}"
    # do the math
    let output_time="10*${restart_index}"
    output_time_string="${output_time}.0"
    echo "    Output time is ${output_time_string}"
    /usr/bin/mpirun --use-hwthread-cpus "-np" "16" "/home/julien/Codes/2023-11-14/PHiLiP/build_release/bin/PHiLiP_3D" "-i" "${restart_file_directory}${filename}"
    python3 /home/julien/Codes/PHiLiP-Post-Processing/cases/turbulent_channel_flow/generate_boundary_profile.py
    target_filename="${flow_field_directory}${target_flow_field_subdir}velocity_vorticity-0_boundary_layer_profile_t0${output_time}_OS-${oversampling_factor}.dat"
    cp ${file_we_copy} ${target_filename}
done


