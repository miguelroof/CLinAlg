
from . cimport AlgMatrix

cdef double getLength(double * wirepoint, unsigned int numvect)

cdef double getLengthByIndex(double * wirepoint, list wireIndex)

cdef bint isInside(double * wirepoint, unsigned int numvect, double * punto, bint incEdges)

cdef bint isClosed(double * wirePoint, int numpoints)

cdef bint isInsideByIndex(double * wirepoint, list wireIndex, double * punto, bint incEdges)

cdef void getRotate(double * toPointer, double * wirePoint, unsigned int numpoints, double * center, double * axis,
                            double angle)

cdef void getCDG(double * cdg, double * wirePoint, unsigned int numpoints)

cdef void getCDGByIndex(double * cdg, double * wirePoint, list wireInd)

cdef bint getPointFromOrigin(double * toPoint, double * wirePoint, unsigned int numpoints, double length)

cdef double getLengthFromOrigin(double * wirePoint, unsigned int numpoints, double * point)

cdef bint getDirectionAtPoint(double * toDir, double * wirePoint, unsigned int numpoints, double * point)

cdef void getNearestPoint(double * toPoint, double * wirePoint, unsigned int numpoints, double * point)

cdef bint isClockWise(double * wirePoint, unsigned int numpoints, double * obsPoint)

cdef list getShifted(double * wirePoint, unsigned int numpoints, double value)

cdef list getIntersectionPointWithWire(double * wirePoint, unsigned int numpoints,
                                       double * otherWirePoint, unsigned int otherNumpoints,
                                       bint incEdge)

cdef list getIntersectionPointWithLine(double * wirePoint, unsigned int numpoints, double * linePoint, double * lineDir, bint incEdge)

cdef list getIntersectionPointWithLineByIndex(double * basePoints, list wireInd, double * linePoint, double * lineDir,
                                           bint incEdge)

cdef list getIntersectionPointWithSegment(double * wirePoint, unsigned int numpoints,
                                          double * vini, double * vfin,
                                          bint incEdge)

cdef Wire getSubWire(double * wirePoint, unsigned int numpoints, double length1, double length2)

cdef void getNormal(double * normal, double * wirePoint, unsigned int numpoints)

cdef void getNormalByIndex(double * normal, double * wirePoint, list wireInd)



cdef class PointIterator():
    cdef AlgMatrix.Matrix mat
    cdef int n

cdef class SegmentIterator():
    cdef AlgMatrix.Matrix mat
    cdef int n

cdef class Wire():
    cdef AlgMatrix.Matrix vectMat
    cpdef double length(self)
