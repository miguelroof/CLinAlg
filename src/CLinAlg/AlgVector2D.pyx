from cpython.mem cimport PyMem_Malloc, PyMem_Free
from libc.math cimport sqrt, acos, sin, cos, fabs
from .AlgTool cimport presition, M_PI



cdef void add(double * todouble, double * v1, double * v2):
    todouble[0] = v1[0] + v2[0]
    todouble[1] = v1[1] + v2[1]

cdef void sub(double * todouble, double * v1, double * v2):
    todouble[0] = v1[0] - v2[0]
    todouble[1] = v1[1] - v2[1]

cdef double dot(double * v1, double * v2):
    return v1[0] * v2[0] + v1[1] * v2[1]

cdef double cross(double * v1, double * v2):
    return v1[0] * v2[1] - v1[1] * v2[0]

cdef double module(double * v):
    if v == NULL:
        raise RuntimeError("Pointer to vector cannot be NULL")
    return sqrt(dot(v, v))

cdef void vdir(double * todouble, double * pini, double * pfin):
    cdef double modul
    sub(todouble, pfin, pini)
    modul = module(todouble)
    if modul:
        todouble[0] = todouble[0] / modul
        todouble[1] = todouble[1] / modul

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
    return sqrt((v1[0] - v2[0]) ** 2 + (v1[1] - v2[1]) ** 2)

cdef bint isEqual(double * v1, double * v2):
    return fabs(v1[0] - v2[0]) <= presition or fabs(v1[1] - v2[1]) <= presition

cdef void copy(double * tto, double * ffrom):
    tto[0] = ffrom[0]
    tto[1] = ffrom[1]

cdef tuple toTuple(double * v):
    return (v[0], v[1])

cdef void project(double * todouble, double * v1, double * v2):
    cdef double v1v2, v2v2
    v2v2 = dot(v2, v2)
    if fabs(v2v2) < presition:
        todouble[0] = 0.0
        todouble[1] = 0.0
    v1v2 = dot(v1, v2) / v2v2
    todouble[0] = v2[0] * v1v2
    todouble[1] = v2[1] * v1v2

cdef class PhantomVector2D():

    def __cinit__(self):
        self._v = NULL
        self._ref_object = None

    def __dealloc__(self):
        if self._ref_object:
            # print("Elimino referencia en vector2D")
            self._ref_object = None
        if self._v != NULL:
            self._v = NULL

    def set_ref_object(self, object):
        self._ref_object = object

    @property
    def module(self) -> float:
        if self._v != NULL:
            return module(self._v)
        else:
            raise RuntimeError("PhantomVector2D not initialized")

    @module.setter
    def module(PhantomVector2D self, float value):
        cdef double m = module(self._v)
        if m:
            self._v[0] = self._v[0] * value / m
            self._v[1] = self._v[1] * value / m

    def angle(PhantomVector2D v1, PhantomVector2D v2)-> float:
        return angle(v1._v, v2._v)

    @property
    def x(PhantomVector2D self) -> float:
        return self._v[0]

    @x.setter
    def x(PhantomVector2D self, value):
        self._v[0] = value

    @property
    def y(PhantomVector2D self) -> float:
        return self._v[1]

    @y.setter
    def y(PhantomVector2D self, value):
        self._v[1] = value

    def normalize(self):
        cdef double m = module(self._v)
        if m:
            self._v[0] = self._v[0] / m
            self._v[1] = self._v[1] / m
        return self

    def copy(PhantomVector2D self):
        return Vector2D(self)

    def distance(PhantomVector2D v1, PhantomVector2D v2) -> float:
        return distance(v1._v, v2._v)

    def rotate(PhantomVector2D v, PhantomVector2D axis_point, angle) -> Vector2D:
        "Rotate 2D Vector around axis_point defined by other 2D Vector and vertical direction"
        cdef double uvectv
        cdef unsigned int i
        cdef Vector2D u = Vector2D()
        cdef double sinangle, cosangle, uv, uu
        cdef double axismodule
        axismodule = module(axis_point._v)
        sinangle = sin(angle * 0.5)
        cosangle = cos(angle * 0.5)
        for i in range(2):
            u._v[i] = sinangle * axis_point[i] / axismodule
        uvectv = cross(u._v, v._v)
        uu = dot(u._v, u._v)
        uv = dot(u._v, v._v)
        for i in range(2):
            u._v[i] = 2 * uv * u._v[i] + (cosangle * cosangle - uu) * v._v[i]
        return u

    def project(PhantomVector2D v1, PhantomVector2D v2):
        "project one Vector2D onto the second"
        cdef Vector2D newVect = Vector2D()
        project(newVect._v, v1._v, v2._v)
        return newVect

    def __serialize__(PhantomVector2D self) -> list:
        cdef unsigned int i
        return [self._v[i] for i in range(2)]

    def toList(PhantomVector2D self) -> list:
        return self.__serialize__()

    def pythonized(PhantomVector2D self) -> list:
        return self.__serialize__()

    def toTuple(PhantomVector2D self) -> tuple:
        return (self._v[0], self._v[1])

    def __hash__(PhantomVector2D self):
        return hash(self.toTuple())

    def __repr__(PhantomVector2D self) -> str:
        return 'Vector2D(%.2f,%.2f)' % (self[0], self[1])

    def __str__(Vector2D self) -> str:
        return 'Vector2D(%.2f,%.2f)' % (self[0], self[1])

    def __getitem__(PhantomVector2D self, key) -> float:
        return self._v[key]

    def __setitem__(PhantomVector2D self, key, value):
        self._v[key] = value

    def __copy__(PhantomVector2D self):
        return self.copy()

    def __deepcopy__(Vector2D self, memo=None):
        return self.copy()

    #....................OPERADORES......................................
    def __add__(PhantomVector2D self, PhantomVector2D other):
        cdef Vector2D newVector = Vector2D()
        add(newVector._v, (<PhantomVector2D> self)._v, other._v)
        return newVector

    def __iadd__(PhantomVector2D self, PhantomVector2D other):  #cuando se produce sumas dentro en un bucle tipo sum
        add(self._v, self._v, other._v)

    def __radd__(PhantomVector2D self, PhantomVector2D other):  #cuando se produce sumas dentro en un bucle tipo sum
        add(self._v, self._v, other._v)

    def __sub__(self, PhantomVector2D other):
        cdef Vector2D  newVector = Vector2D()
        sub(newVector._v, (<PhantomVector2D> self)._v, other._v)
        return newVector

    def __isub__(PhantomVector2D self, PhantomVector2D other):
        sub(self._v, self._v, other._v)

    def __rsub__(PhantomVector2D self, PhantomVector2D other):
        sub(self._v, self._v, other._v)

    def __eq__(PhantomVector2D self, PhantomVector2D other) -> bool:
        return bool(isEqual(self._v, other._v))

    def __ne__(PhantomVector2D self, PhantomVector2D other):
        return not self.__eq__(other)

    def __neg__(PhantomVector2D self):
        cdef Vector2D newVector = Vector2D()
        newVector._v[0] = -self._v[0]
        newVector._v[1] = -self._v[1]
        return newVector

    def __mul__(PhantomVector2D self, other):
        "dot product"
        cdef Vector2D newVector
        if isinstance(other, PhantomVector2D):
            return dot(self._v, (<PhantomVector2D> other)._v)
        elif isinstance(other, (int, float)):
            newVector = Vector2D()
            newVector._v[0] = self._v[0] * other
            newVector._v[1] = self._v[1] * other
            return newVector
        else:
            raise RuntimeError("Not allowed this argument for vector multiplication")

    def __rmul__(PhantomVector2D self, other):
        "dot product"
        cdef Vector2D newVector
        if isinstance(other, PhantomVector2D):
            return dot(self._v, (<PhantomVector2D> other)._v)
        elif isinstance(other, (int, float)):
            newVector = Vector2D()
            newVector._v[0] = self._v[0] * other
            newVector._v[1] = self._v[1] * other
            return newVector
        else:
            raise RuntimeError("Not allowed this argument for vector multiplication")

    def __bool__(PhantomVector2D self)-> bool:
        return (abs(self._v[0])>presition or fabs(self._v[1])>presition)

    def __truediv__(PhantomVector2D self, other):
        cdef Vector2D newVector = Vector2D()
        if isinstance(other, (int, float)) and other != 0:
            newVector._v[0] = self._v[0] / other
            newVector._v[1] = self._v[1] / other
            return newVector
        else:
            del newVector
            raise ArithmeticError(self, other)

    def __mod__(PhantomVector2D self, PhantomVector2D other) -> float:
        "cross product"
        return cross(self._v, other._v)

cdef class Vector2D(PhantomVector2D):
    def __cinit__(self, *coord, **kwargs):
        if len(coord) == 1 and coord[0] is None:
            pass
        else:
            self._v = <double *> PyMem_Malloc(2 * sizeof(double))

    def __init__(self, *coord, **kwargs):
        if len(coord) == 1:
            if isinstance(coord[0], PhantomVector2D):
                self._v[0] = (<PhantomVector2D> coord[0])._v[0]
                self._v[1] = (<PhantomVector2D> coord[0])._v[1]
                return
            elif isinstance(coord[0], (list, tuple)) and len(coord[0]) == 2:
                self._v[0] = coord[0][0]
                self._v[1] = coord[0][1]
            elif coord[0] is None:  # esto lo usare cuando quiera reservar memoria pero no inicializarla
                return
                # raise NotImplementedError("You cannot use None as argument")
        elif len(coord) == 2:
            self._v[0] = coord[0]
            self._v[1] = coord[1]
        elif len(coord) == 0:
            return
        else:
            raise NotImplementedError("You cannot initialize a vector with %d arguments" % len(coord))

    def __dealloc__(self):
        if self._v:
            PyMem_Free(self._v)
            self._v = NULL

    def __json__(self):
        return {'__jsoncls__': 'CLinAlg.AlgVector:Vector2D.from_JSON', 'vector': (self._v[0], self._v[1])}

    @classmethod
    def from_JSON(cls, jsondict):
        cdef Vector2D newVector = Vector2D(jsondict['vector'])
        return newVector
    
