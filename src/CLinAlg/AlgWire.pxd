from . cimport AlgMatrix

cdef double getLength(double * wirepoint, int * indexes, unsigned int numvect)

cdef bint isInside(double * wirepoint, int * indexes, unsigned int numvect, double * punto, bint incEdges)

cdef bint isClosed(double * wirePoint, int * indexes, int numpoints)

cdef void getRotate(double * toPointer, double * wirePoint, unsigned int numpoints, double * center, double * axis,
                            double angle)

cdef void getCDG(double * cdg, double * wirePoint, int * indexes, unsigned int numpoints)

cdef bint getPointFromOrigin(double * toPoint, double * wirePoint, int * indexes, unsigned int numpoints, double length)

cdef double getLengthFromOrigin(double * wirePoint, int * indexes, unsigned int numpoints, double * point)

cdef bint getDirectionAtPoint(double * toDir, double * wirePoint, int * indexes, unsigned int numpoints, double * point)

cdef void getNearestPoint(double * toPoint, double * wirePoint, int * indexes, unsigned int numpoints, double * point)

cdef bint isClockWise(double * wirePoint, int * indexes, unsigned int numpoints, double * obsPoint)

cdef list getShifted(double * wirePoint, int * indexes, unsigned int numpoints, double value)

cdef list getIntersectionPointWithWire(double * wirePoint, int * indexes, unsigned int numpoints,
                                       double * otherWirePoint,  int * otherIndexes, unsigned int otherNumpoints,
                                       bint incEdge)

cdef list getIntersectionPointWithLine(double * wirePoint, int * indexes, unsigned int numpoints, double * linePoint, double * lineDir,
                                           bint incEdge)

cdef list getIntersectionPointWithSegment(double * wirePoint, int * indexes, unsigned int numpoints,
                                          double * vini, double * vfin,
                                          bint incEdge)

cdef Wire getSubWire(double * wirePoint, int * indexes, unsigned int numpoints, double length1, double length2)

cdef void getNormal(double * normal, double * wirePoint, int * indexes, unsigned int numpoints)




cdef class IndexPerimeter():
    cdef int * _m
    cdef unsigned int _npoint
    cdef unsigned int niterator
    cpdef void setList(IndexPerimeter self, list cdata)

cdef class PointIterator():
    cdef AlgMatrix.Matrix vectMat
    cdef IndexPerimeter indPer
    cdef unsigned int niterator

cdef class SegmentIterator():
    cdef AlgMatrix.Matrix vectMat
    cdef IndexPerimeter indPer
    cdef unsigned int niterator


cdef class Wire():
    cdef AlgMatrix.Matrix vectMat
    cdef IndexPerimeter indPer
    cpdef double length(self)
    cpdef list getBounds(self)


