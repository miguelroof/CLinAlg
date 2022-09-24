from cpython.mem cimport PyMem_Malloc, PyMem_Free
from libc.math cimport fabs
from . cimport AlgVector, AlgLine, AlgMatrix
from .AlgTool cimport presition

cdef bint isInside(double * pini, double * pfin, double * point, bint incEdge):
    cdef unsigned int i, j, k
    cdef double mnu
    for i in range(3):
        if pini[i] == pfin[i]:
            continue
        mnu = (point[i] - pini[i]) / (pfin[i] - pini[i])
        if incEdge:
            if mnu < -presition or mnu > 1 + presition:
                return False
        else:
            if mnu < presition or mnu > 1 - presition:
                return False
        j = (i + 1) % 3
        k = ((i + j) * 2) % 3
        if (fabs(pini[j] + mnu * (pfin[j] - pini[j]) - point[j]) < presition and
                fabs(pini[k] + mnu * (pfin[k] - pini[k]) - point[k]) < presition):
            return True
        else:
            return False
    return False

cdef void  getNearestPoint(double * toPoint, double * pini, double * pfin, double * point):
    cdef double angle
    cdef unsigned int i
    cdef double * cdir = <double *> PyMem_Malloc(3 * sizeof(double))
    AlgLine.projectedPoint(toPoint, pini, cdir, point)
    AlgVector.sub(cdir, pfin, pini)
    try:
        if isInside(pini, pfin, toPoint, (<bint>True)):
            return
        else:
            if AlgVector.distance(pini, point) < AlgVector.distance(pfin, point):
                for i in range(3):
                    toPoint[i] = pini[i]
            else:
                for i in range(3):
                    toPoint[i] = pfin[i]
    finally:
        PyMem_Free(cdir)

cdef bint getIntersectionPointWithLine(double * toPoint, double * pini, double * pfin,
                                           double * lpnt, double * ldir, bint incEdge):
    # las lineas no pueden ser paralelas
    cdef double * cdir = <double *> PyMem_Malloc(3 * sizeof(double))
    AlgVector.sub(cdir, pfin, pini)
    try:
        if AlgLine.getIntersectionPoint(toPoint, pini, cdir, lpnt, ldir) and isInside(pini, pfin, toPoint, <bint>True):
            return True
        else:
            return False
    finally:
        PyMem_Free(cdir)

cdef bint getIntersectionPointWithSegment(double * toPoint, double * pini1, double * pfin1,
                                              double * pini2, double * pfin2, bint incEdge):
    cdef double * cdir1 = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * cdir2 = <double *> PyMem_Malloc(3 * sizeof(double))
    AlgVector.vdir(cdir1, pini1, pfin1)
    AlgVector.vdir(cdir2, pini2, pfin2)
    try:
        if AlgLine.getIntersectionPoint(toPoint, pini1, cdir1, pini2, cdir2) and isInside(pini1, pfin1, toPoint, <bint>True) and isInside(pini2, pfin2,
                                                                        toPoint,
                                                                        <bint>True):
            return True
        else:
            return False
    finally:
        PyMem_Free(cdir1)
        PyMem_Free(cdir2)



cdef class Segment():

    def __cinit__(self, p_ini, p_fin):
        cdef unsigned int i
        self._vini = <double *> PyMem_Malloc(3 * sizeof(double))
        self._vfin = <double *> PyMem_Malloc(3 * sizeof(double))
        for i in range(3):
            self._vini[i] = (<AlgVector.Vector>p_ini)._v[i]
            self._vfin[i] = (<AlgVector.Vector>p_fin)._v[i]

    def __dealloc__(Segment self):
        PyMem_Free(self._vini)
        PyMem_Free(self._vfin)

    def __json__(self):
        cdef unsigned int i
        pini = [self._vini[i] for i in range(3)]
        pfin = [self._vfin[i] for i in range(3)]
        return {'__jsoncls__': 'CLinAlg.AlgSegment:Segment.from_JSON', 'v1': pini, 'v2': pfin}

    @classmethod
    def from_JSON(cls, jsondict):
        obj = cls(jsondict['v1'], jsondict['v2'])
        return obj

    @property
    def pini(Segment self):
        cdef AlgVector.PhantomVector v = AlgVector.PhantomVector()
        v._v = self._vini
        return v

    @pini.setter
    def pini(Segment self, p_ini:AlgVector.Vector):
        cdef unsigned int i
        for i in range(3):
            self._vini[i] = p_ini._v[i]

    @property
    def pfin(Segment self):
        cdef AlgVector.PhantomVector v = AlgVector.PhantomVector()
        v._v = self._vfin
        return v

    @pfin.setter
    def pfin(Segment self, p_fin:AlgVector.Vector):
        cdef unsigned int i
        for i in range(3):
            self._vfin[i] = p_fin._v[i]

    @property
    def line(Segment self):
        cdef AlgLine.Line newLine
        cdef unsigned int i
        cdef double * cdir = <double *> PyMem_Malloc (3 * sizeof(double))
        AlgVector.vdir(cdir, self._vini, self._vfin)
        newLine = AlgLine.Line(None)
        newLine._pnt = self._vini
        newLine._dir = cdir
        return newLine

    def isInside(Segment self, point:AlgVector.Vector, incEdge=True):
        return isInside(self._vini, self._vfin, point._v, <bint>incEdge)

    def getNearestPoint(Segment self, point:AlgVector.Vector):
        cdef AlgVector.Vector newVect = AlgVector.Vector(None)
        newVect._v = <double *> PyMem_Malloc (3 * sizeof(double))
        getNearestPoint(newVect._v, self._vini, self._vfin, point._v)
        return newVect

    def getIntersectionPoint(Segment self, other, incEdge=True):
        cdef AlgVector.Vector newVect = AlgVector.Vector(None)
        newVect._v = <double *> PyMem_Malloc(3 * sizeof(double))
        if isinstance(other, AlgLine.Line):
            if getIntersectionPointWithLine(newVect._v, self._vini, self._vfin, (<AlgLine.Line> other)._pnt,
                                                                      (<AlgLine.Line> other)._dir,
                                                                      <bint>incEdge):
                return newVect
            else:
                return None
        elif isinstance(other, Segment):
            if getIntersectionPointWithSegment(newVect._v, self._vini, self._vfin,
                                                                         (<Segment> other)._vini,
                                                                         (<Segment> other)._vfin, <bint>incEdge):
                return newVect
            else:
                return None

    def isIntersect(Segment self, Segment other, incEdge=True):
        cdef double * ipoint = <double *> PyMem_Malloc (3* sizeof(double))
        try:
            return getIntersectionPointWithSegment(ipoint, self._vini, self._vfin, (<Segment> other)._vini,
                                                                     (<Segment> other)._vfin, <bint>incEdge)
        finally:
            PyMem_Free(ipoint)

    def length(Segment self):
        return AlgVector.distance(self._vini, self._vfin)

    def copy(Segment self):
        return Segment(self.pini, self.pfin)

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memodict={}):
        return self.copy()

    def __repr__(self):
        return 'Segment(' + str(AlgVector.toTuple(self._vini)) + ',' + str(AlgVector.toTuple(self._vfin)) + ')'

    def __str__(self):
        return 'Segment(' + str(AlgVector.toTuple(self._vini)) + ',' + str(AlgVector.toTuple(self._vfin)) + ')'

    def __eq__(self, other):
        return ((self.pini == other.pini and self.pfin == other.pfin) or
                (self.pini == other.pfin and self.pfin == other.pini))

    def __ne__(self, other):
        return not self == other


cdef class PhantomSegment(Segment):
    def __cinit__(self, *args):
        if args:
            self.setPointer(*args)

    def setPointer(self, *args):
        # 2 BaseVectores
        # 2 matrices y dos posiciones
        if len(args) == 2 and isinstance(args[0], AlgVector.Vector) and isinstance(args[1], AlgVector.Vector):
            self._vini = (<AlgVector.Vector> args[0])._v
            self._vfin = (<AlgVector.Vector> args[1])._v
        elif len(args) == 4 and isinstance(args[0], AlgMatrix.Matrix) and isinstance(args[2], AlgMatrix.Matrix):
            self._vini = &(<AlgMatrix.Matrix> args[0])._m[args[1]]
            self._vfin = &(<AlgMatrix.Matrix> args[2])._m[args[3]]

    def __dealloc__(self):
        pass  # no borro lo referenciado