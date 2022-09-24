from .AlgMatrix cimport Matrix

cdef list tessellate(double * surfPoint, unsigned int numpoints)

cdef double getArea(Surface surf)

cdef void getCDG(double * cdg, Surface surf)

cdef Matrix getInertia(Surface surf)

cdef double staticMoment(Surface surf, double * point, double * dir)

cdef bint isInside(Surface surf, double * point, bint incEdge)

cdef tuple dictIntersectByLine(Surface surf, double * point, double * dir, bint incEdge)

cdef (double, double, double, double) boundRectangle(Surface surf, int v1, int v2)

cdef Matrix stenierInGlobalAxis_p(Matrix I, double area, double * cdg, double * to_point)

cpdef tuple mainAxisInertia(Matrix tensor)


cdef class MatrixTriangle():
    cdef int * _m
    cdef unsigned int _rows
    cpdef void setList(MatrixTriangle self, list data)

cdef class Surface():
    cdef Matrix vectMat
    cdef MatrixTriangle triMat