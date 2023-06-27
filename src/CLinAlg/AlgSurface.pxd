from .AlgMatrix cimport Matrix
from .AlgWire cimport IndexPerimeter
from .AlgTransformation cimport Transformation

cdef list tessellate(double * surfPoint, int * indPer, unsigned int numPer, double * axis)

cdef void sortTrianglesByNormal(double * surfPoint, int * indTri, int numTri, double * direction)

cdef double getArea(Surface surf)

cdef void getCDG(double * cdg, Surface surf)

cdef Matrix getInertia(Surface surf)

cdef double staticMoment(Surface surf, double * point, double * dir)

cdef bint isInside(Surface surf, double * point, bint incEdge)

cdef tuple dictIntersectByLine(Surface surf, double * point, double * dir, bint incEdge)

# cdef int _splitTriangleByLine(double * v1, double * v2, double * wirePoint, int * indexTriangle, double * linepnt, double * linedir)
#
# cdef list _getPerimeterByTriangleIndex(list indTri)
#
# cdef list _getSurfaceByTriangleIndex(list triIndex)

cdef list splitByLine(Surface surf, double * linepnt, double * linedir)


cdef (double, double, double, double) boundRectangle(Surface surf, int v1, int v2)

cdef Matrix stenierInGlobalAxis_p(Matrix I, double area, double * cdg, double * to_point)

cpdef tuple mainAxisInertia(Matrix tensor)

cpdef list trapezeToSurfaceYZ(list trapeces)

cdef class IndexTriangle():
    cdef int * _m
    cdef unsigned int _ntriangle
    cpdef void setList(IndexTriangle self, list data)

cdef class EdgeIterator():
    cdef Matrix vectMat
    cdef IndexPerimeter indPer
    cdef unsigned int niterator

cdef class Surface():
    cdef Matrix vectMat
    cdef IndexTriangle indTri
    cdef IndexPerimeter indPer
    cpdef Surface transform(Surface self, Transformation transf)