from .AlgMatrix cimport Matrix
from .AlgVector cimport Vector

cdef void quatFromAxisAngle(double * toquat, double * axis, double angle)

cdef void cross(double * todouble, double * q1, double * q2)

cdef double dot(double * q1, double * q2)

cdef void div(double * todouble, double * q1, double * q2)

cdef void add(double * todouble, double * q1, double * q2)

cdef void sub(double * todouble, double * q1, double * q2)

cdef double module(double * q)

cdef void rotateVectorMatrix(double * toPointer, double * q, double * center, double * mat, unsigned int numvect)

cdef Quaternion fromRotationMatrix(Matrix m)

cdef Quaternion fromEulerAngles(double R, double P, double Y)

cdef class Quaternion():
    cdef double * _q
    cdef Vector c_getAxis(Quaternion self)
    cdef Matrix c_getRotationMatrix(Quaternion self)
    cdef void c_rotateVector_p(Quaternion self, double * topointer, double * vpointer)
    cdef (double, double, double) c_eulerAngles(Quaternion self)