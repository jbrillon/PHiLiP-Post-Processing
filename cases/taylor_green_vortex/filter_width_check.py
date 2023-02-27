import numpy as np

# inputs
nquad = 6.0 # P+1, P=5
n_cells_per_direction = 16.0 # 8 for 48^3 DOFs, 16 for 96^3 DOFs

# calculation
domain_length = 2.0*np.pi
domain_volume = domain_length*domain_length*domain_length # 3D
cell_length = domain_length/n_cells_per_direction
n_cells_per_volume = n_cells_per_direction*n_cells_per_direction*n_cells_per_direction # 3D
cell_volume = domain_volume/n_cells_per_volume

# quadrature points; reference: https://keisan.casio.com/exec/system/1280801905
if(nquad == 6.0):
    nodes = np.array([-1, -0.765055323929464692851, -0.2852315164806450963142, 0.2852315164806450963142, 0.765055323929464692851, 1],dtype=np.float64)

spacing_between_nodes = np.diff(nodes)
min_spacing_between_nodes = spacing_between_nodes[0] # first two points are always closest
jacobian_mapping_function = 0.5*cell_length
filter_width_based_on_smallest_spacing = min_spacing_between_nodes*jacobian_mapping_function

old_filter_width = cell_volume/(nquad*nquad*nquad)
print("Old filter width: %1.6e" % old_filter_width)

new_filter_width = cell_length/nquad
print("New filter width: %1.6e" % new_filter_width)

print("Filter witdth based on minimum node spacing: %1.6e" % filter_width_based_on_smallest_spacing)

print("Input to PHiLiP code such that it cancels the current filter width: %1.15f" % (filter_width_based_on_smallest_spacing/new_filter_width))
