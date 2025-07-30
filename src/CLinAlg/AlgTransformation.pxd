from .AlgMatrix cimport Matrix
from .AlgVector cimport Vector
from .AlgQuaternion cimport Quaternion




cdef class Transformation():
    cdef Matrix _mat
    cdef Transformation c_inverse(Transformation self)
    cdef Vector c_transformVector(Transformation self, Vector v)
    cdef list c_transformVectorList(Transformation self, list vList)
    cpdef Transformation copy(Transformation self, bint withRotation, bint withTranslation)

cdef Transformation fromRotTrans(Quaternion quat, Vector vect)
cdef Transformation fromPointAndAxis(Vector pos, Vector axis_x, Vector axis_y)