from cpython.mem cimport PyMem_Malloc, PyMem_Free
from libc.math cimport sqrt, M_PI, acos,  sin, cos, fabs
from AlgTool cimport presition
from AlgMatrix cimport Matrix


cdef void dir(double * todouble, double * pini, double * pfin):
    cdef double modulo
    cdef unsigned int i
    sub(todouble, pfin, pini)
    modulo = module(todouble)
    if modulo:
        for i in range(3):
            todouble[i] = todouble[i] / modulo

cdef void add(double * todouble, double * v1, double * v2):
    cdef unsigned int i
    for i in range(3):
        todouble[i] = v1[i] + v2[i]

cdef void sub(double * todouble, double * v1, double * v2):
    cdef unsigned int i
    for i in range(3):
        todouble[i] = v1[i] - v2[i]

cdef double dot(double * v1, double * v2):
    cdef unsigned int i
    cdef double val = 0
    for i in range(3):
        val += v1[i] * v2[i]
    return val

cdef double module(double * v):
    return sqrt(dot(v, v))

cdef double angle(double * v1, double * v2):
    cdef double X
    X = dot(v1, v2)
    if fabs(X) > presition:
        X = X / (module(v1) * module(v2))
    else:
        return M_PI * 0.5
    if X < -1: X = -2 - X
    if X > 1: X = 2 - X
    return acos(X)

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
    "Producto BaseVectorial"
    todouble[0] = v1[1] * v2[2] - v1[2] * v2[1]
    todouble[1] = v1[2] * v2[0] - v1[0] * v2[2]
    todouble[2] = v1[0] * v2[1] - v1[1] * v2[0]

cdef double * clone(double * ffrom):
    cdef double * todouble = <double *> PyMem_Malloc(3 * sizeof(double))
    copy(todouble, ffrom)
    return todouble

cdef tuple totuple(double * v):
    return (v[0], v[1], v[2])

cdef void project(double * todouble, double * v1, double * v2):
    cdef double v1v2, v2v2
    v2v2 = dot(v2, v2)
    if v2v2 == 0:
        todouble[0] = 0
        todouble[1] = 0
        todouble[2] = 0
    v1v2 = dot(v1, v2) / v2v2
    todouble[0] = v2[0] * v1v2
    todouble[1] = v2[1] * v1v2
    todouble[2] = v2[2] * v1v2

cdef class Vector():

    def __cinit__(self, *coord):
        if coord != (None,):
            self._v = <double *> PyMem_Malloc(3 * sizeof(double))
        if len(coord) == 1:
            if isinstance(coord[0], Vector):
                self._v[0] = (<Vector>coord[0])._v[0]
                self._v[1] = (<Vector>coord[0])._v[1]
                self._v[2] = (<Vector>coord[0])._v[2]
                return
            elif isinstance(coord[0], (list, tuple)) and len(coord[0]) == 3:
                self._v[0] = coord[0][0]
                self._v[1] = coord[0][1]
                self._v[2] = coord[0][2]
            elif coord[0] is None:  # esto lo usare cuando quiera reservar memoria pero no llenarlo
                return
                # raise NotImplementedError("You cannot use None as argument")
        elif len(coord) == 3:
            self._v[0] = coord[0]
            self._v[1] = coord[1]
            self._v[2] = coord[2]
        elif len(coord) == 0:
            self._v[0] = 0
            self._v[1] = 0
            self._v[2] = 0
            return
        else:
            raise NotImplementedError("You cannot initialize a vector with %d arguments" % len(coord))
        print("Called cinit")

    def __dealloc__(self):
        PyMem_Free(self._v)
        print("Called dealloc")

    @property
    def module(self)-> float:
        return module(self._v)

    @module.setter
    def module(Vector self, value):
        cdef double m = module(self._v)
        if m:
            self._v[0] = self._v[0] * value / m
            self._v[1] = self._v[1] * value / m
            self._v[2] = self._v[2] * value / m

    def angle(Vector v1, Vector v2)-> float:
        return angle(v1._v, v2._v)

    @property
    def x(Vector self) -> float:
        return self._v[0]

    @x.setter
    def x(Vector self, value):
        self._v[0] = value

    @property
    def y(Vector self) -> float:
        return self._v[1]

    @y.setter
    def y(Vector self, value):
        self._v[1] = value

    @property
    def z(Vector self) -> float:
        return self._v[2]

    @z.setter
    def z(Vector self, value):
        self._v[2] = value

    def normalize(self)-> Vector:
        cdef double m = module(self._v)
        if m:
            self._v[0] = self._v[0] / m
            self._v[1] = self._v[1] / m
            self._v[2] = self._v[2] / m
        return self

    def copy(Vector self) -> Vector:
        return Vector(self)

    def distance(Vector v1, Vector v2) -> float:
        return distance(v1._v, v2._v)

    def rotate(Vector v, Vector axis, angle) -> Vector:
        cdef double * uvectv = <double *> PyMem_Malloc(3 * sizeof(double))
        cdef unsigned int i
        cdef Vector u = Vector()
        cdef double sinangle, cosangle, uv, uu
        try:
            sinangle = sin(angle * 0.5)
            cosangle = cos(angle * 0.5)
            for i in range(3):
                u._v[i] = sinangle * axis[i] / axis.module
            cross(uvectv, u._v, v._v)
            uu = dot(u._v, u._v)
            uv = dot(u._v, v._v)
            for i in range(3):
                u._v[i] = 2 * uv * u._v[i] + (cosangle * cosangle - uu) * v._v[i] + 2 * cosangle * uvectv[i]
            return u
        finally:
            PyMem_Free(uvectv)

    def project(Vector v1, Vector v2) -> Vector:
        "project one Vector onto the second"
        cdef Vector newVect = Vector()
        project(newVect._v, v1._v, v2._v)
        return newVect

    def toList(Vector self) -> list:
        cdef unsigned int i
        return [self._v[i] for i in range(3)]

    def __hash__(Vector self):
        return hash(tuple(self.toList()))

    def __repr__(Vector self) -> str:
        return 'Vector(%.2f,%.2f,%.2f)' % (self[0], self[1], self[2])

    def __str__(Vector self) -> str:
        return 'Vector(%.2f,%.2f,%.2f)' % (self[0], self[1], self[2])

    def __getitem__(Vector self, key) -> float:
        return self._v[key]

    def __setitem__(Vector self, key, value):
        self._v[key] = value

    def __copy__(Vector self) -> Vector:
        return self.copy()

    def __deepcopy__(Vector self, memo=None)-> Vector:
        return self.copy()

    #....................OPERADORES......................................
    def __add__(Vector self, Vector other) -> Vector:
        cdef Vector newVector = Vector()
        add(newVector._v, self._v, other._v)
        return newVector

    def __iadd__(Vector self, Vector other):  #cuando se produce sumas dentro en un bucle tipo sum
        self._v[0] = self._v[0] + other._v[0]
        self._v[1] = self._v[1] + other._v[1]
        self._v[2] = self._v[2] + other._v[2]

    def __sub__(Vector self, Vector other) -> Vector:
        cdef Vector  newVector = Vector()
        newVector._v[0] = self._v[0] - other._v[0]
        newVector._v[1] = self._v[1] - other._v[1]
        newVector._v[2] = self._v[2] - other._v[2]
        return newVector

    def __isub__(Vector self, Vector other):
        self._v[0] = self._v[0] - other._v[0]
        self._v[1] = self._v[1] + - other._v[1]
        self._v[2] = self._v[2] - other._v[2]

    def __eq__(Vector self, Vector other) -> bool:
        return isEqual(self._v, other._v)

    def __ne__(Vector self, Vector other):
        return not self.__eq__(other)

    def __neg__(Vector self) -> Vector:
        cdef Vector newVector = Vector()
        newVector._v[0] = -self._v[0]
        newVector._v[1] = -self._v[1]
        newVector._v[2] = -self._v[2]
        return newVector

    def __mul__(Vector self, other):
        "dot product"
        cdef Vector newVector
        # print(type(other))
        if isinstance(other, Vector) and isinstance(self, Vector):
            return dot(self._v, (<Vector> other)._v)
        elif isinstance(other, (int, float)):
            newVector = Vector()
            newVector._v[0] = self._v[0] * other
            newVector._v[1] = self._v[1] * other
            newVector._v[2] = self._v[2] * other
            return newVector
        else:
            raise IOError(other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __bool__(Vector self)-> bool:
        return True if (self._v[0] or self._v[1] or self._v[2]) else False

    def __truediv__(Vector self, other):
        cdef Vector newVector = Vector()
        if isinstance(other, (int, float)) and other != 0:
            newVector._v[0] = self._v[0] / other
            newVector._v[1] = self._v[1] / other
            newVector._v[2] = self._v[2] / other
            return newVector
        else:
            del newVector
            raise ArithmeticError(self, other)

    def __mod__(Vector self, Vector other) -> Vector:
        cdef Vector newVector = Vector()
        cross(newVector._v, self._v, other._v)
        return newVector

    def __json__(self):
        return {'__jsoncls__': 'CLinAlg.AlgVector:Vector.from_JSON', 'vector': (self._v[0], self._v[1], self._v[2])}

    @classmethod
    def from_JSON(cls, jsondict):
        cdef Vector newVector = Vector(jsondict['vector'])
        return newVector


cdef class PhantomVector(Vector):
    def __cinit__(self,*args):
        pass

    def setPointer(PhantomVector self, *args):
        if len(args) == 1 and isinstance(args[0], Vector):
            self._v = (<Vector> args[0])._v
        elif len(args) == 2 and isinstance(args[0], Matrix) and isinstance(args[1], int):
            self._v = &(<Matrix> args[0])._m[args[1]]

    def __dealloc__(PhantomVector self):
        pass  # no quiero que lo elimine dado que solo es una referencia
