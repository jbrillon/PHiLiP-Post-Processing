#-----------------------------------------------------
def get_DOF_vars(nElements_per_direction,poly_degree):
    nElements = nElements_per_direction*nElements_per_direction*nElements_per_direction
    nQuadPoints_per_element = poly_degree + 1
    nQuadPoints = nQuadPoints_per_element*nElements_per_direction
    nDOF = nQuadPoints*nQuadPoints*nQuadPoints
    reduced_nQuadPoints = nElements_per_direction*nQuadPoints_per_element - (nElements_per_direction-1)
    reduced_nDOF = reduced_nQuadPoints*reduced_nQuadPoints*reduced_nQuadPoints
    return [nElements,nQuadPoints_per_element,nQuadPoints,nDOF,reduced_nQuadPoints,reduced_nDOF]
#-----------------------------------------------------
def get_reduced_nDOF_and_nQuadPoints(nElements_per_direction,poly_degree):
    reduced_nQuadPoints = get_DOF_vars(nElements_per_direction,poly_degree)[-2]
    reduced_nDOF = get_DOF_vars(nElements_per_direction,poly_degree)[-1]
    return [reduced_nQuadPoints,reduced_nDOF]
#-----------------------------------------------------