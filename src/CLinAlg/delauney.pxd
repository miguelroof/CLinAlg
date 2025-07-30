
cdef double * normalize_3d_to_2d(double * vect3darray,  int num_point, double * zero_3d, double * axe_3d_u, double * axe_3d_v)

cdef bint is_inside(double * vA, double * vB, double * vC, double * vP)

cdef bint should_swap(double * vCom, double * vA, double * vB, double * vP)