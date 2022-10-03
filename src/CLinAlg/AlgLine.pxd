
cdef void projectedPoint(double * toPointer, double * pnt, double * dir, double * fromPoint)

cdef double distanceToPoint(double * pnt, double * dir, double * point)

cdef bint isInside(double * pnt, double * dir, double * point)

cdef bint isPerpendicular(double * pnt1, double * dir1, double * pnt2, double * dir2)

cdef bint isParallel(double * pnt1, double * dir1, double * pnt2, double * dir2)

cdef bint getIntersectionPoint(double * toPoint, double * pnt1, double * dir1, double * pnt2, double * dir2)

cdef void getParametricPoint(double * toPoint, double * pnt, double * dir, double s)

cdef (bint, double) getParameter(double * pnt, double * dir, double * point)

cdef class PhantomLine():
    cdef double * _pnt
    cdef double * _dir

cdef class Line(PhantomLine):
    pass