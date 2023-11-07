from cpython.mem cimport PyMem_Malloc, PyMem_Free
from libc.math cimport sqrt, M_PI, acos, sin, cos, fabs
from .AlgTool cimport presition

cdef void vdir(double * todouble, double * pini, double * pfin):
    cdef double modulo
    cdef unsigned int i
    sub(todouble, pfin, pini)
    modulo = module(todouble)
    if modulo:
        for i in range(3):
            todouble[i] = todouble[i] / modulo

cdef void add(double * todouble, double * v1, double * v2):
    todouble[0] = v1[0] + v2[0]
    todouble[1] = v1[1] + v2[1]
    todouble[2] = v1[2] + v2[2]

cdef void sub(double * todouble, double * v1, double * v2):
    todouble[0] = v1[0] - v2[0]
    todouble[1] = v1[1] - v2[1]
    todouble[2] = v1[2] - v2[2]

cdef double dot(double * v1, double * v2):
    cdef double val = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]
    return val

cdef double module(double * v):
    return sqrt(dot(v, v))

cdef double angle(double * v1, double * v2):
    cdef double x
    x = dot(v1, v2)
    if fabs(x) > presition:
        x = x / (module(v1) * module(v2))
    else:
        return M_PI * 0.5
    if x < -1: x = -2 - x
    if x > 1: x = 2 - x
    return acos(x)

cdef double distance(double * v1, double * v2):
    cdef double val = 0
    cdef unsigned int i
    for i in range(3):
        val += (v1[i] - v2[i]) * (v1[i] - v2[i])
    return sqrt(val)

cdef bint isEqual(double * v1, double * v2):
    cdef unsigned int i
    for i in range(3):
        if fabs(v1[i] - v2[i]) > presition:
            return False
    return True

cdef void copy(double * tto, double * ffrom):
    cdef unsigned int i
    for i in range(3):
        tto[i] = ffrom[i]

cdef void cross(double * todouble, double * v1, double * v2):
    """Producto BaseVectorial"""
    todouble[0] = v1[1] * v2[2] - v1[2] * v2[1]
    todouble[1] = v1[2] * v2[0] - v1[0] * v2[2]
    todouble[2] = v1[0] * v2[1] - v1[1] * v2[0]

cdef tuple toTuple(double * v):
    return (v[0], v[1], v[2])

cdef void project(double * todouble, double * v1, double * v2):
    cdef double v1v2, v2v2
    v2v2 = dot(v2, v2)
    if fabs(v2v2) < presition:
        todouble[0] = 0
        todouble[1] = 0
        todouble[2] = 0
    v1v2 = dot(v1, v2) / v2v2
    todouble[0] = v2[0] * v1v2
    todouble[1] = v2[1] * v1v2
    todouble[2] = v2[2] * v1v2

cdef class PhantomVector():
    @property
    def module(self) -> float:
        return module(self._v)

    @module.setter
    def module(PhantomVector self, float value):
        cdef double m = module(self._v)
        if m:
            self._v[0] = self._v[0] * value / m
            self._v[1] = self._v[1] * value / m
            self._v[2] = self._v[2] * value / m

    def angle(PhantomVector v1, PhantomVector v2)-> float:
        return angle(v1._v, v2._v)

    @property
    def x(PhantomVector self) -> float:
        return self._v[0]

    @x.setter
    def x(PhantomVector self, value):
        self._v[0] = value

    @property
    def y(PhantomVector self) -> float:
        return self._v[1]

    @y.setter
    def y(PhantomVector self, value):
        self._v[1] = value

    @property
    def z(PhantomVector self) -> float:
        return self._v[2]

    @z.setter
    def z(PhantomVector self, value):
        self._v[2] = value

    def normalize(self):
        cdef double m = module(self._v)
        if m:
            self._v[0] = self._v[0] / m
            self._v[1] = self._v[1] / m
            self._v[2] = self._v[2] / m
        return self

    def copy(PhantomVector self):
        return Vector(self)

    def distance(PhantomVector v1, PhantomVector v2) -> float:
        return distance(v1._v, v2._v)

    def rotate(PhantomVector v, PhantomVector axis, angle):
        cdef double * uvectv
        cdef unsigned int i
        cdef Vector u = Vector()
        cdef double sinangle, cosangle, uv, uu
        cdef double axismodule
        uvectv = <double *> PyMem_Malloc(3 * sizeof(double))
        try:
            axismodule = module(axis._v)
            sinangle = sin(angle * 0.5)
            cosangle = cos(angle * 0.5)
            for i in range(3):
                u._v[i] = sinangle * axis[i] / axismodule
            cross(uvectv, u._v, v._v)
            uu = dot(u._v, u._v)
            uv = dot(u._v, v._v)
            for i in range(3):
                u._v[i] = 2 * uv * u._v[i] + (cosangle * cosangle - uu) * v._v[i] + 2 * cosangle * uvectv[i]
            return u
        finally:
            PyMem_Free(uvectv)

    def project(PhantomVector v1, PhantomVector v2):
        "project one Vector onto the second"
        cdef Vector newVect = Vector()
        project(newVect._v, v1._v, v2._v)
        return newVect

    def toList(PhantomVector self) -> list:
        cdef unsigned int i
        return [self._v[i] for i in range(3)]

    def toTuple(PhantomVector self) -> tuple:
        return (self._v[0], self._v[1], self._v[2])

    def __hash__(PhantomVector self):
        return hash(self.toTuple())

    def __repr__(PhantomVector self) -> str:
        return 'Vector(%.2f,%.2f,%.2f)' % (self[0], self[1], self[2])

    def __str__(Vector self) -> str:
        return 'Vector(%.2f,%.2f,%.2f)' % (self[0], self[1], self[2])

    def __getitem__(PhantomVector self, key) -> float:
        return self._v[key]

    def __setitem__(PhantomVector self, key, value):
        self._v[key] = value

    def __copy__(PhantomVector self):
        return self.copy()

    def __deepcopy__(Vector self, memo=None):
        return self.copy()

    #....................OPERADORES......................................
    def __add__(self, PhantomVector other):
        cdef Vector newVector = Vector()
        if isinstance(self, PhantomVector):
            add(newVector._v, (<PhantomVector> self)._v, other._v)
        else:
            copy(newVector._v, other._v)
        return newVector

    def __iadd__(PhantomVector self, PhantomVector other):  #cuando se produce sumas dentro en un bucle tipo sum
        add(self._v, self._v, other._v)

    def __radd__(self, PhantomVector other):  #cuando se produce sumas dentro en un bucle tipo sum
        add(self._v, self._v, other._v)

    def __sub__(self, PhantomVector other):
        cdef Vector  newVector = Vector()
        if isinstance(self, PhantomVector):
            sub(newVector._v, (<PhantomVector> self)._v, other._v)
        else:
            copy(newVector._v, other._v)
        return newVector

    def __isub__(PhantomVector self, PhantomVector other):
        sub(self._v, self._v, other._v)

    def __rsub__(PhantomVector self, PhantomVector other):
        sub(self._v, self._v, other._v)

    def __eq__(PhantomVector self, PhantomVector other) -> bool:
        return bool(isEqual(self._v, other._v))

    def __ne__(PhantomVector self, PhantomVector other):
        return not self.__eq__(other)

    def __neg__(PhantomVector self):
        cdef Vector newVector = Vector()
        newVector._v[0] = -self._v[0]
        newVector._v[1] = -self._v[1]
        newVector._v[2] = -self._v[2]
        return newVector

    def __mul__(PhantomVector self, other):
        "dot product"
        cdef Vector newVector
        if isinstance(other, PhantomVector):
            return dot(self._v, (<PhantomVector> other)._v)
        elif isinstance(other, (int, float)):
            newVector = Vector()
            newVector._v[0] = self._v[0] * other
            newVector._v[1] = self._v[1] * other
            newVector._v[2] = self._v[2] * other
            return newVector
        else:
            raise RuntimeError("Not allowed this argument for vector multiplication")

    def __rmul__(PhantomVector self, other):
        "dot product"
        cdef Vector newVector
        if isinstance(other, PhantomVector):
            return dot(self._v, (<PhantomVector> other)._v)
        elif isinstance(other, (int, float)):
            newVector = Vector()
            newVector._v[0] = self._v[0] * other
            newVector._v[1] = self._v[1] * other
            newVector._v[2] = self._v[2] * other
            return newVector
        else:
            raise RuntimeError("Not allowed this argument for vector multiplication")

    def __bool__(PhantomVector self)-> bool:
        cdef unsigned int i
        for i in range(3):
            if fabs(self._v[i]) > presition:
                return True
        return False

    def __truediv__(PhantomVector self, other):
        cdef Vector newVector = Vector()
        if isinstance(other, (int, float)) and other != 0:
            newVector._v[0] = self._v[0] / other
            newVector._v[1] = self._v[1] / other
            newVector._v[2] = self._v[2] / other
            return newVector
        else:
            del newVector
            raise ArithmeticError(self, other)

    def __mod__(PhantomVector self, PhantomVector other):
        cdef Vector newVector = Vector()
        cross(newVector._v, self._v, other._v)
        return newVector

cdef class Vector(PhantomVector):
    def __cinit__(self, *coord, **kwargs):
        if coord == (None):
            pass
        else:
            self._v = <double *> PyMem_Malloc(3 * sizeof(double))

    def __init__(self, *coord, **kwargs):
        if len(coord) == 1:
            if isinstance(coord[0], PhantomVector):
                self._v[0] = (<PhantomVector> coord[0])._v[0]
                self._v[1] = (<PhantomVector> coord[0])._v[1]
                self._v[2] = (<PhantomVector> coord[0])._v[2]
                return
            elif isinstance(coord[0], (list, tuple)) and len(coord[0]) == 3:
                self._v[0] = coord[0][0]
                self._v[1] = coord[0][1]
                self._v[2] = coord[0][2]
            elif coord[0] is None:  # esto lo usare cuando quiera reservar memoria pero no inicializarla
                return
                # raise NotImplementedError("You cannot use None as argument")
        elif len(coord) == 3:
            self._v[0] = coord[0]
            self._v[1] = coord[1]
            self._v[2] = coord[2]
        elif len(coord) == 0:
            return
        else:
            raise NotImplementedError("You cannot initialize a vector with %d arguments" % len(coord))

    def __dealloc__(self):
        PyMem_Free(self._v)

    def __json__(self):
        return {'__jsoncls__': 'CLinAlg.AlgVector:Vector.from_JSON', 'vector': (self._v[0], self._v[1], self._v[2])}

    @classmethod
    def from_JSON(cls, jsondict):
        cdef Vector newVector = Vector(jsondict['vector'])
        return newVector
