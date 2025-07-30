from cpython.mem cimport PyMem_Malloc, PyMem_Free
from libc.math cimport cos
from .AlgTool cimport presition, M_PI
from .AlgMatrix cimport Matrix
from . cimport AlgVector
from .AlgVector cimport Vector


cdef double offsetToPoint(Plane plane, double * p, double * direction):
    cdef double * punto = <double *> PyMem_Malloc(3 * sizeof(double))
    try:
        if punto == NULL:
            raise MemoryError()
        AlgVector.sub(punto, plane._pos, p)
        return AlgVector.dot(direction, punto)
    finally:
        if punto != NULL:
            PyMem_Free(punto)
            punto = NULL
        
cdef void  getLocalCoords(double * toPoint, Plane plane, double * point):
    cdef double x, y, z, angle
    cdef double * pt = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * xv = <double *> PyMem_Malloc(3 * sizeof(double))
    try:
        if pt == NULL or xv == NULL:
            raise MemoryError()
        AlgVector.sub(pt, point, plane._pos)
        AlgVector.project(xv, pt, plane._u)
        x = AlgVector.module(xv)
        angle = AlgVector.angle(xv, plane._u)
        if angle > M_PI * 0.5:
            x = -x
        AlgVector.project(xv, pt, plane._v)
        y = AlgVector.module(xv)
        angle = AlgVector.angle(xv, plane._v)
        if angle > M_PI * 0.5:
            y = -y
        AlgVector.project(xv, pt, plane._axis)
        z = AlgVector.module(xv)
        angle = AlgVector.angle(xv, plane._axis)
        if angle > M_PI * 0.5:
            z = -z
        toPoint[0] = x
        toPoint[1] = y
        toPoint[2] = z
    finally:
        if pt != NULL:
            PyMem_Free(pt)
            pt = NULL
        if xv != NULL:
            PyMem_Free(xv)
            xv = NULL
        
cdef void getGlobalCoords(double * toPoint, Plane plane, double * point):
    cdef unsigned int i
    for i in range(3):
        toPoint[i] = plane._u[i] * point[0] + plane._v[i] * point[1] + plane._axis[i] * point[
            2] + plane._pos[i]

cdef void projectPoint(double * toPoint, Plane plane, double * point, double * direction):
    cdef double a, hyp, dirmod
    cdef double * vproj = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * lp = <double *> PyMem_Malloc (3 * sizeof(double))
    cdef double * gp = <double *> PyMem_Malloc (3 * sizeof(double))
    cdef unsigned int i
    # """ Alternativa: v = point-orig
    # dist = vx*nx + vy*ny + vz*nz
    # projected_point = point - dist*n  donde n es la normal. Esta proyeccion da el punto en 3d"""
    try:
        if vproj == NULL or lp == NULL or gp == NULL:
            raise MemoryError()
        dirmod = AlgVector.module(direction)
        getLocalCoords(lp, plane, point)
        lp[2] = 0.0
        getGlobalCoords(gp, plane, lp)
        AlgVector.sub(vproj, gp, point)
        a = AlgVector.angle(direction, vproj)
        if a > M_PI * 0.5:
            a = M_PI - a
            hyp = AlgVector.module(vproj) / cos(a)
            for i in range(3):
                toPoint[i] = point[i] - direction[i] * hyp / dirmod
        else:
            hyp = AlgVector.module(vproj) / cos(a)
            for i in range(3):
                toPoint[i] = point[i] + direction[i] * hyp / dirmod
    finally:
        # noinspection PyTypeChecker
        if vproj != NULL:
            PyMem_Free(vproj)
            vproj = NULL
        if lp != NULL:
            PyMem_Free(lp)
            lp = NULL
        if gp != NULL:
            PyMem_Free(gp)
            gp = NULL

cdef bint isInside(Plane plane, double * point):
    cdef unsigned int i
    cdef double suma = 0
    for i in range(3):
        suma += (point[i] - plane._pos[i]) * plane._axis[i]
    return suma < presition

cdef void system(double * mat3x3, Plane plane):
    cdef unsigned int i
    for i in range(3):
        mat3x3[i] = plane._u[i]
        mat3x3[i+3] = plane._v[i]
        mat3x3[i+6] = plane._axis[i]

cdef class Plane:
    '''A plane with local axis'''

    def __cinit__(self, u:Vector, v:Vector, axis:Vector, pos:Vector):
        # print("he creado un Plano")
        self._u = <double *> PyMem_Malloc(3 * sizeof(double))
        self._v = <double *> PyMem_Malloc(3 * sizeof(double))
        self._axis = <double *> PyMem_Malloc(3 * sizeof(double))
        self._pos = <double *> PyMem_Malloc(3 * sizeof(double))

    def __init__(self, u:Vector, v:Vector, axis:Vector, pos:Vector):
        cdef unsigned int i
        cdef double vmod
        vmod = AlgVector.module(u._v)
        if vmod < presition:
            raise RuntimeError("U vector for plane cannot be zero")
        for i in range(3):
            self._u[i] = u._v[i] / vmod
        vmod = AlgVector.module(v._v)
        if vmod < presition:
            raise RuntimeError("V vector for plane cannot be zero")
        for i in range(3):
            self._v[i] = v._v[i] / vmod
        vmod = AlgVector.module(axis._v)
        if vmod < presition:
            raise RuntimeError("Axis vector for plane cannot be zero")
        for i in range(3):
            self._axis[i] = axis._v[i] / vmod
        for i in range(3):
            self._pos[i] = pos._v[i]

            
    def __dealloc__(self):
        if self._u:
            PyMem_Free(self._u)
            self._u = NULL
        if self._v:
            PyMem_Free(self._v)
            self._v = NULL
        if self._axis:
            PyMem_Free(self._axis)
            self._axis = NULL
        if self._pos:
            PyMem_Free(self._pos)
            self._pos = NULL

    #..................C METHODS................................
    
    def offsetToPoint(self, p:Vector, direction= None):
        if direction is None:
            direction = self.axis
        else:
            direction = direction.normalize()
        return offsetToPoint(self, p._v, (<Vector>direction)._v)
        
    def getLocalCoords(self, point:Vector):
        cdef Vector v = Vector()
        getLocalCoords(v._v, self, point._v)
        return v
    
    def getGlobalCoords(self, point:Vector):
        cdef Vector v = Vector()
        getGlobalCoords(v._v, self, point._v)
        return v
            
    def projectPoint(self, p:Vector, direction=None):
        cdef Vector v = Vector()
        if direction is None:
            direction = self.axis
        projectPoint(v._v, self, p._v, (<Vector>direction)._v)
        return v

    #....................PYTHON METHODS......................
    @property
    def u(self):
        cdef AlgVector.PhantomVector pv = AlgVector.PhantomVector()
        pv._v = self._u
        pv._ref_object = self # evitar que lo borren
        return pv

    @u.setter
    def u(self, Vector u):
        cdef unsigned int i
        cdef double vmod
        vmod = AlgVector.module(u._v)
        if vmod < presition:
            raise RuntimeError("u vector cannot be zero")
        for i in range(3):
            self._u[i] = u._v[i] / vmod

    @property
    def v(self):
        cdef AlgVector.PhantomVector pv = AlgVector.PhantomVector()
        pv._v = self._v
        pv._ref_object = self  # evitar que lo borren
        return pv

    @v.setter
    def v(self, Vector v):
        cdef unsigned int i
        cdef double vmod
        vmod = AlgVector.module(v._v)
        if vmod < presition:
            raise RuntimeError("v vector cannot be zero")
        for i in range(3):
            self._v[i] = v._v[i] / vmod

    @property
    def axis(self):
        cdef AlgVector.PhantomVector pv = AlgVector.PhantomVector()
        pv._v = self._axis
        pv._ref_object = self  # evitar que lo borren
        return pv

    @axis.setter
    def axis(self, Vector axis):
        cdef unsigned int i
        cdef double vmod
        vmod = AlgVector.module(axis._v)
        if vmod < presition:
            raise RuntimeError("axis vector cannot be zero")
        for i in range(3):
            self._axis[i] = axis._v[i] / vmod

    @property
    def position(self):
        cdef AlgVector.PhantomVector pv = AlgVector.PhantomVector()
        pv._v = self._pos
        pv._ref_object = self # evitar que lo borren
        return pv

    @position.setter
    def position(self, pos:Vector):
        cdef unsigned int i
        for i in range(3):
            self._pos[i] = pos._v[i]

    def __repr__(self):
        return "Workplane u=" + str(AlgVector.toTuple(self._u)) + " v=" + str(AlgVector.toTuple(self._v)) + " axis=" + str(
            AlgVector.toTuple(self._axis)) + " orig=" + str(AlgVector.toTuple(self._pos))

    def alignToPointAndAxis(self, point:Vector, axis:Vector, offset=0.0):
        cdef double ax, bx, ay, pby, az, bz, dmod
        cdef double * pvect = <double *> PyMem_Malloc(3 * sizeof(double))
        cdef double * pvect2 = <double *> PyMem_Malloc(3 * sizeof(double))
        cdef unsigned int i
        try:
            if pvect == NULL or pvect2 == NULL:
                raise MemoryError()
            dmod = AlgVector.module(axis._v)
            if dmod < presition:
                raise RuntimeError("axis vector cannot be zero")
            for i in range(3):
                self._axis[i] = axis._v[i] / dmod
            pvect[0] = 1
            pvect[1] = 0
            pvect[2] = 0
            ax = AlgVector.angle(self._axis, pvect)
            bx = M_PI - ax
            pvect[0] = 0
            pvect[1] = 1
            ay = AlgVector.angle(self._axis, pvect)
            pby = M_PI - ay
            pvect[1] = 0
            pvect[2] = 1
            az = AlgVector.angle(self._axis, pvect)
            bz = M_PI - az
            b = min(ax, ay, az, bx, pby, bz)
            if b in [ax, bx]:
                pvect2[0] = 0
                pvect2[1] = 1
                pvect2[0] = 0
            else:
                pvect2[0] = 1
                pvect2[1] = 0
                pvect2[0] = 0
            AlgVector.cross(self._v, pvect2, self._axis)
            AlgVector.cross(self._u, self._axis, self._v)
            for i in range(3):
                self._pos[i] = point._v[i] + offset * self._axis[i]
        finally:
            if pvect != NULL:
                PyMem_Free(pvect)
                pvect = NULL
            if pvect2 != NULL:
                PyMem_Free(pvect2)
                pvect2 = NULL

    def isInside(self, point:Vector):
        return isInside(self, point._v)

    @property
    def system(self):
        cdef Matrix newMat = Matrix()
        newMat._m = <double *> PyMem_Malloc (9 * sizeof(double))
        if not newMat._m:
            raise MemoryError()
        newMat._rows = 3
        newMat._cols = 3
        system(newMat._m, self)
        return newMat

    def from3Points(self, v1:Vector, v2:Vector, v3:Vector, isCounterClockwise=True):
        if isCounterClockwise:
            axis = (v3 - v2) % (v2 - v1)
            self.alignToPointAndAxis(v1, axis, 0)
        else:
            axis = (v2 - v1) % (v3 - v2)
            self.alignToPointAndAxis(v1, axis, 0)

    def getClosestAxis(self, point:Vector):
        cdef double ax, ay, az, bx, pby, bz, b
        ax = AlgVector.angle(point._v, self._u)
        ay = AlgVector.angle(point._v, self._v)
        az = AlgVector.angle(point._v, self._axis)
        bx = M_PI - ax
        pby = M_PI - ay
        bz = M_PI - az
        b = min(ax, ay, az, bx, pby, bz)
        if b in [ax, bx]:
            return self.u
        elif b in [ay, pby]:
            return self.v
        elif b in [az, bz]:
            return self.axis
        else:
            return None

    def getRotation(self):
        return self.system.transpose()

    def inverse(self):
        cdef unsigned int i
        for i in range(3):
            self._u[i] = -self._u[i]
            self._axis[i] = -self._axis[i]

    def copy(self):
        return Plane(self.u, self.v, self.axis, self.position)

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memodict={}):
        return self.copy()

    def __json__(self):
        return {'__jsoncls__': "CLinAlg.AlgPlane:Plane.from_JSON", 'u': self.u, 'v': self.v,
                'axis': self.axis, 'position': self.position}

    @classmethod
    def from_JSON(cls, jsondict):
        obj = cls(**jsondict.fromkeys(['u', 'v', 'axis', 'position']))
        return obj