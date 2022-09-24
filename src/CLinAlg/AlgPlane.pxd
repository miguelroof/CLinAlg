
cdef double offsetToPoint(Plane plane, double * p, double * direction)

cdef void getLocalCoords(double * toPoint, Plane plane, double * point)

cdef void getGlobalCoords(double * toPoint, Plane plane, double * point)

cdef void projectPoint(double * toPoint, Plane plane, double * point, double * direction)

cdef bint isInside(Plane plane, double * point)

cdef void system(double * mat3x3, Plane plane)

cdef class Plane():
    cdef double * _u
    cdef double * _v
    cdef double * _axis
    cdef double * _pos
