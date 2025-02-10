restartIndices=({48..67})
# restartIndices=("49")

for i in ${!restartIndices[@]}; do
    restart_index=${restartIndices[${i}]}
    filename="restart-000${restart_index}.prm"
    echo "Replacing text in file ${filename}"
    # do the math
    let output_time="10*${restart_index}"
    output_time_string="${output_time}.001"
    end_time_string="${output_time}.002"
    echo "    Output time is ${output_time_string}"

    replaceStrings=(\
        "  set final_time                                             = " \
        # "  set restart_files_directory_name                           = restart_files                      # default: ." \
        # "    set output_density_field_in_addition_to_velocity             = false" \
        # "    set output_flow_field_files_directory_name                   = ." \
        # "    set output_velocity_field_at_fixed_times                     = false" \
        # "    set output_velocity_field_times_string                       =       # default:" \
        # "    set output_velocity_number_of_subvisions                     = 2" \
        # "    set output_viscosity_field_in_addition_to_velocity           = false" \
        # "    set output_vorticity_magnitude_field_in_addition_to_velocity = false" \
    )

    newStrings=(\
        "  set final_time                                             = ${end_time_string}" \
        # "  set restart_files_directory_name                           = restart_files_cPlus_C4" \
        # "    set output_density_field_in_addition_to_velocity             = true" \
        # "    set output_flow_field_files_directory_name                   = flow_field_files" \
        # "    set output_velocity_field_at_fixed_times                     = true" \
        # "    set output_velocity_field_times_string                       = ${output_time_string}" \
        # "    set output_velocity_number_of_subvisions                     = 1" \
        # "    set output_viscosity_field_in_addition_to_velocity           = true" \
        # "    set output_vorticity_magnitude_field_in_addition_to_velocity = true" \
    )

    # now I need a way to loop through the indices and do a multiplication by 10 and change the filename
    for j in ${!replaceStrings[@]}; do
        replaceString=${replaceStrings[$j]}
        newString=${newStrings[$j]}
        sed -i -e "s/.*${replaceString}.*/${newString}/" ${filename}
    done
    
    # for changing the oversampling, uncomment and comment the loop above
    # replaceString="    set output_velocity_number_of_subvisions                     = 1"
    # newString="    set output_velocity_number_of_subvisions                     = 2"
    # sed -i -e "s/.*${replaceString}.*/${newString}/" ${filename}

    echo "    Done replacing text in file ${filename}."
    rm "${filename}-e"
    echo " "
done
