cimport AlgMatrix
from AlgMatrix cimport Matrix
from AlgVector cimport Vector
from AlgQuaternion cimport Quaternion, fromRotationMatrix

cpdef Transformation fromRotTrans(Quaternion quat, Vector vect)

cdef class Transformation():
    cdef Matrix _mat
    cdef Transformation c_inverse(Transformation self)
    cdef Vector c_transformVector(Transformation self, Vector v)
    cdef list c_transformVectorList(Transformation self, list vList)
