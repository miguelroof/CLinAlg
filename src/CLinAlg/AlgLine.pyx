from cpython.mem cimport PyMem_Malloc, PyMem_Free
from .cimport AlgVector
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
    cdef double * projpoint = <double *> PyMem_Malloc(3 * sizeof(double))

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

cdef void getParametricPoint(double * toPoint, double * pnt, double * dir, double s):
    """
    Returns parametric point for a line
    :param toPoint: return poitner to point
    :param pnt: pointer with point within line
    :param dir: pointer with direction of line
    :param s: parameter
    :return: void
    """
    cdef unsigned int i
    for i in range(3):
        toPoint[i] = pnt[i] + s * dir[i]

cdef (bint, double) getParameter(double * pnt, double * dir, double * point):
    """
    Return parameter for a given point. Return True and the parameter if the point is within the line, otherwise
    return False and 0
    :param pnt: pointer with point of line
    :param dir: pointer with direction of line
    :param point: pointer to get the parameter
    :return: double of parametric line
    """
    cdef unsigned int i
    if isInside(pnt, dir, point):
        for i in range(3):
            if dir[i] != 0:
                return True, (point[i] - pnt[i]) / dir[i]
        return False, 0
    else:
        return False, 0

    #.............................C METHODS..........................................................from

cdef class PhantomLine():
    # .............................PYTHON METHODS.............................................................

    @property
    def point(self):
        cdef PhantomVector pnt = PhantomVector()
        pnt._v = self._pnt
        return pnt

    @property
    def direction(self):
        cdef PhantomVector dir = PhantomVector()
        dir._v = self._dir
        return dir

    def projectedPoint(self, point: PhantomVector):
        cdef Vector vaux = Vector()
        projectedPoint(vaux._v, self._pnt, self._dir, point._v)
        return vaux

    def distanceToPoint(self, point: PhantomVector):
        return distanceToPoint(self._pnt, self._dir, point._v)

    def getParametricPoint(self, s: float):
        cdef Vector newVect = Vector()
        getParametricPoint(newVect._v, self._pnt, self._dir, s)
        return newVect

    def isInside(self, point: PhantomVector):
        return isInside(self._pnt, self._dir, point._v)

    def copy(self):
        return Line(self)

    def getParameter(self, point: PhantomVector):
        param = getParameter(self._pnt, self._dir, point._v)
        if param[0]:
            return param[1]
        else:
            return None

    def isPerpendicular(self, other: PhantomLine):
        return isPerpendicular(self._pnt, self._dir, other._pnt, other._dir)

    def isParallel(self, other: PhantomLine):
        return isParallel(self._pnt, self._dir, other._pnt, other._dir)

    def getIntersectionPoint(self, other: PhantomLine):
        cdef Vector newVect = Vector()
        if not getIntersectionPoint(newVect._v, self._pnt, self._dir, other._pnt, other._dir):
            return None
        return newVect

    def __eq__(self, other: PhantomLine) -> bool:
        if self.isParallel(other) and self.isInside(other.point):
            return True
        return False

    def __ne__(self, other: PhantomLine):
        return not self.__eq__(other)

    def __repr__(self) -> str:
        return "Line(Pnt:" + str(AlgVector.toTuple(self._pnt)) + "; Dir:" + str(AlgVector.toTuple(self._dir)) + ")"

    def __hash__(self):
        cdef unsigned int i
        cdef tuple tup
        tup = tuple(self._pnt[i] for i in range(3)) + tuple(self._dir[i] for i in range(3))
        return hash(tup)

cdef class Line(PhantomLine):
    def __cinit__(self, *args, **kwargs):
        if args != (None):
            self._pnt = <double *> PyMem_Malloc(3 * sizeof(double))
            self._dir = <double *> PyMem_Malloc(3 * sizeof(double))

    def __init__(Line self, *args, **kwargs):
        """args:
            - Line
            - Point1, Point2
        """
        cdef unsigned int i
        if len(args) == 1:
            if isinstance(args[0], PhantomLine):
                AlgVector.copy(self._pnt, (<Line> args[0])._pnt)
                AlgVector.copy(self._dir, (<Line> args[0])._dir)
            elif args[0] is None:  # no hago reserva de memoria
                pass
        elif len(args) == 2 and isinstance(args[0], PhantomVector) and isinstance(args[1], PhantomVector):
            AlgVector.copy(self._pnt, (<PhantomVector> args[0])._v)
            AlgVector.vdir(self._dir, (<PhantomVector> args[0])._v, (<PhantomVector> args[1])._v)

    def __dealloc__(Line self):
        PyMem_Free(self._pnt)
        PyMem_Free(self._dir)

    def __json__(Line self) -> dict:
        return {'__jsoncls__': 'CLinAlg.AlgLine:Line.from_JSON', 'pnt': AlgVector.toTuple(self._pnt),
                'dir': AlgVector.toTuple(self._dir)}

    @classmethod
    def from_JSON(cls, jsondict):
        cdef unsigned int i
        obj = Line()
        for i in range(3):
            obj._pnt[i] = jsondict['pnt'][i]
            obj._dir[i] = jsondict['dir'][i]
        return obj
