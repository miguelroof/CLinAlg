from libc.math cimport sqrt, M_PI, acos, sin, cos, fabs
from .AlgTool cimport presition

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






