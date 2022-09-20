from AlgVector cimport Vector, PhantomVector
from AlgTool cimport presition

cdef double * projectedPoint(double * pnt, double * dir, double * fromPoint)

cdef double distanceToPoint(double * pnt, double * dir, double * point)

cdef bint isInside(double * pnt, double * dir, double * point)

cdef bint isPerpendicular(double * pnt1, double * dir1, double * pnt2, double * dir2)

cdef bint isParallel(double * pnt1, double * dir1, double * pnt2, double * dir2)

cdef double * getIntersectionPoint(double * pnt1, double * dir1, double * pnt2, double * dir2)

cdef class Line():
    cdef double * _pnt
    cdef double * _dir