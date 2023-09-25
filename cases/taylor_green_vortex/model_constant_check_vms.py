import numpy as np

num_el_per_dir = 16.0

poly_degree_large_scales = 3.0
poly_degree = 5.0
p_ratio = poly_degree_large_scales/poly_degree
curve_fit_constant = 1.174
mesh_size = 2.0*np.pi/(num_el_per_dir)
# print(mesh_size)
cDG = 0.172*mesh_size/poly_degree
# print(cDG)
cSSVMS = 1.0 - (curve_fit_constant * p_ratio)**(4.0/3.0)
cSSVMS = (cSSVMS)**(-3.0/4.0)
cSSVMS *= cDG 
cSSVMS_sqr = cSSVMS * cSSVMS
delta = mesh_size / (poly_degree+1.0)
cSmag = 0.1
cSmagDelta_sqr = (cSmag * delta)**2.0
cDG_sqr = cDG*cDG
print("cDG_sqr: %e" % cDG_sqr)
print("cSmagDelta_sqr: %e" % cSmagDelta_sqr)
print("cSSVMS_sqr: %e" % cSSVMS_sqr)
cFladCUSP_sqr = (0.33 * delta)**2.0
if(poly_degree==7.0):
    print("cFladCUSP_sqr: %e" % cFladCUSP_sqr)

# flad_Cs = 9.21*np.exp(-3.03*poly_degree_large_scales/poly_degree)
# print(flad_Cs)