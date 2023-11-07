cdef void vdir(double * todouble, double * pini, double * pfin)

cdef void add(double * todouble, double * v1, double * v2)

cdef void sub(double * todouble, double * v1, double * v2)

cdef double dot(double * v1, double * v2)

cdef double module(double * v)

cdef double angle(double * v1, double * v2)

cdef double distance(double * v1, double * v2)

cdef bint isEqual(double * v1, double * v2)

cdef void copy(double * tto, double * ffrom)

cdef double cross(double * v1, double * v2)

cdef tuple toTuple(double * v)

cdef void project(double * todouble, double * v1, double * v2)

cdef class PhantomVector2D():
    cdef double * _v

cdef class Vector2D(PhantomVector2D):
    pass