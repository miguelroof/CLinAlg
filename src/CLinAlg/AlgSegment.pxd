
cdef bint isInside(double * pini, double * pfin, double * point, bint incEdge)

cdef void  getNearestPoint(double * toPoint, double * pini, double * pfin, double * point)

cdef bint getIntersectionPointWithLine(double * toPoint, double * pini, double * pfin,
                                           double * lpnt, double * ldir, bint incEdge)

cdef bint getIntersectionPointWithSegment(double * toPoint, double * pini1, double * pfin1,
                                              double * pini2, double * pfin2, bint incEdge)

cdef class PhantomSegment():
    cdef double * _vini
    cdef double * _vfin


cdef class Segment(PhantomSegment):
    pass
