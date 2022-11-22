#===================================================
# INPUTS:
#---------------------------------------------------
narval_directory=${1}
target_directory=${2}
#---------------------------------------------------
# gets only the .txt files and avoids creating empty directories
rsync -zarvm --include="*/" --include="*.txt" --exclude="*" "brillon@narval.computecanada.ca:${narval_directory}" "${target_directory}"
rsync -zarvm --include="*/" --include="*.dat" --exclude="*" "brillon@narval.computecanada.ca:${narval_directory}" "${target_directory}"
# rsync -zarv  --prune-empty-dirs --include "*/"  --include="*.txt" --exclude="*" brillon@narval.computecanada.ca:scratch/2023_AIAA/2022-09-06/ .
# Example usage:
# ./get_files.sh "scratch/2023_AIAA/2022-09-06/" "/Users/Julien/julien_phd/post_processing/data/taylor_green_vortex/2022-09-06/"
# ./data/taylor_green_vortex/get_files.sh "scratch/2023_AIAA/2022-09-20/viscous_ILES_cDG_KG_two_point_flux_with_l2roe_dissipation_dofs0256_p3_procs4096" "/Users/Julien/julien_phd/post_processing/data/taylor_green_vortex/2022-09-06/"
# ./get_files.sh "scratch/2023_AIAA/2022-10-04/" "/Users/Julien/julien_phd/post_processing/data/taylor_green_vortex/2022-10-04/"
