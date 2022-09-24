from cpython.mem cimport PyMem_Malloc, PyMem_Free
from . cimport AlgVector
from .AlgTool cimport presition
from .AlgVector cimport Vector, PhantomVector
from libc.math cimport M_PI, fabs

cdef void projectedPoint(double * toPointer, double * pnt, double * dir, double * fromPoint):
    cdef double val
    AlgVector.sub(toPointer, fromPoint, pnt)
    val = AlgVector.dot(toPointer, dir)
    for i in range(3):
        toPointer[i] = pnt[i] + val * dir[i]

cdef double distanceToPoint(double * pnt, double * dir, double * point):
    cdef double * projpoint = <double *> PyMem_Malloc (3 * sizeof(double))

    try:
        projectedPoint(projpoint, pnt, dir, point)
        return AlgVector.distance(projpoint, point)
    finally:
        PyMem_Free(projpoint)

cdef bint isInside(double * pnt, double * dir, double * point):
    cdef double * v1 = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * v2 = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double modulo
    try:
        AlgVector.sub(v1, point, pnt)
        AlgVector.cross(v2, v1, dir)
        modulo = AlgVector.module(v2)
        return modulo <= presition
    finally:
        PyMem_Free(v1)
        PyMem_Free(v2)

cdef bint isPerpendicular(double * pnt1, double * dir1, double * pnt2, double * dir2):
    cdef double angle
    angle = AlgVector.angle(dir1, dir2)
    if fabs(angle - M_PI * 0.5) < presition:
        return True
    return True

cdef bint isParallel(double * pnt1, double * dir1, double * pnt2, double * dir2):
    cdef double angle = AlgVector.angle(dir1, dir2)
    if fabs(angle) < presition or fabs(angle - M_PI) < presition:
        return True
    return False

cdef bint getIntersectionPoint(double * toPoint, double * pnt1, double * dir1, double * pnt2, double * dir2):
    cdef double lamb, mnu
    cdef unsigned int i, j, k
    for i in range(3):
        if dir1[i] != 0:
            for j in range(3):
                if j == i:
                    continue
                mnu = dir2[j] - dir1[j] * dir2[i] / dir1[i]
                if mnu == 0:
                    continue
                mnu = (pnt1[j] + dir1[j] * pnt2[i] / dir1[i] - dir1[j] * pnt1[i] / dir1[i] - pnt2[j]) / mnu
                lamb = (pnt2[i] + dir2[i] * mnu - pnt1[i]) / dir1[i]
                k = ((i + j) * 2) % 3
                if fabs(pnt1[k] + dir1[k] * lamb - pnt2[k] - dir2[k] * mnu) <= presition:
                    toPoint[0] = pnt1[0] + dir1[0] * lamb
                    toPoint[1] = pnt1[1] + dir1[1] * lamb
                    toPoint[2] = pnt1[2] + dir1[2] * lamb
                    return True
    return False

cdef class Line():
    def __init__(Line self, *args):
        """args:
            - Line
            - Point1, Point2
        """
        cdef unsigned int i
        if len(args) == 1:
            if isinstance(args[0], Line):
                self._pnt = <double *> PyMem_Malloc (3 * sizeof(double))
                AlgVector.copy(self._pnt, (<Line> args[0])._pnt)
                self._dir = <double *> PyMem_Malloc (3 * sizeof(double))
                AlgVector.copy(self._dir, (<Line> args[0])._dir)
            elif args[0] is None:  # no hago reserva de memoria
                pass
        elif len(args) == 2 and isinstance(args[0], Vector) and isinstance(args[1], Vector):
            self._pnt = <double *> PyMem_Malloc(3 * sizeof(double))
            self._dir = <double *> PyMem_Malloc(3 * sizeof(double))
            AlgVector.copy(self._pnt, (<Vector> args[0])._v)
            AlgVector.vdir(self._dir, (<Vector> args[0])._v, (<Vector> args[1])._v)

    def __dealloc__(Line self):
        PyMem_Free(self._pnt)
        PyMem_Free(self._dir)

    #.............................C METHODS..........................................................from


    # .............................PYTHON METHODS.............................................................

    @property
    def Point(self)-> PhantomVector:
        pnt = PhantomVector()
        pnt._v = self._pnt
        return pnt

    @property
    def Direction(self)-> PhantomVector:
        dir = PhantomVector()
        dir._v = self._dir
        return dir

    def projectedPoint(Line self, Vector point):
        cdef Vector vaux = Vector(None)
        vaux._v = <double *> PyMem_Malloc (3 * sizeof(double))
        projectedPoint(vaux._v, self._pnt, self._dir, point._v)
        return vaux

    def distanceToPoint(Line self, Vector point)-> float:
        return distanceToPoint(self._pnt, self._dir, point._v)

    def getParametricPoint(Line self, float s):
        cdef Vector newVect = Vector()
        cdef unsigned int i
        for i in range(3):
            newVect._v[i] = self._pnt[i] + s * self._dir[i]
        return newVect

    def isInside(Line self, Vector point):
        return isInside(self._pnt, self._dir, point._v)

    def copy(Line self):
        return Line(self)

    def getParameter(Line self, Vector point, checkInside=True):
        cdef unsigned int i
        if checkInside and not self.isInside(point):
            return None
        for i in range(3):
            if self._dir[i] != 0 and (point._v[i] - self._pnt[i]) != 0:
                return (point._v[i] - self._pnt[i]) / self._dir[i]

    def isPerpendicular(Line self, Line other):
        return isPerpendicular(self._pnt, self._dir, other._pnt, other._dir)

    def isParallel(Line self, Line other):
        return isParallel(self._pnt, self._dir, other._pnt, other._dir)

    def getIntersectionPoint(Line self, Line other):
        cdef Vector newVect = Vector(None)
        newVect._v = <double *> PyMem_Malloc (3 * sizeof(double))
        if not getIntersectionPoint(newVect._v, self._pnt, self._dir, other._pnt, other._dir):
            return None
        return newVect

    def __eq__(Line self, Line other) -> bool:
        if self.isParallel(other) and self.isInside(other.Point):
            return True
        return False

    def __ne__(Line self, Line other):
        return not self.__eq__(other)

    def __repr__(Line self) -> str:
        return "Line(Pnt:" + str(AlgVector.toTuple(self._pnt)) + "; Dir:" + str(AlgVector.toTuple(self._dir)) + ")"

    def __json__(Line self) -> dict:
        return {'__jsoncls__': 'CLinAlg.AlgLine:Line.from_JSON', 'pnt': AlgVector.toTuple(self._pnt), 'dir': AlgVector.toTuple(self._dir)}

    @classmethod
    def from_JSON(cls, jsondict):
        cdef unsigned int i
        obj = Line(None)
        obj._pnt = <double *> PyMem_Malloc(3 * sizeof(double))
        obj._dir = <double *> PyMem_Malloc(3 * sizeof(double))
        for i in range(3):
            obj._pnt[i] = jsondict['pnt'][i]
            obj._dir[i] = jsondict['dir'][i]
        return obj