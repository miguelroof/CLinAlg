from cpython.mem cimport PyMem_Malloc, PyMem_Free
from libc.math cimport sqrt, sin, cos, acos, atan2, asin
from . cimport AlgVector


cdef void quatFromAxisAngle(double * toquat, double * axis, double angle):
    cdef double vmod = sin(angle / 2) / AlgVector.module(axis)
    toquat[0] = cos(angle / 2)
    toquat[1] = axis[0] * vmod
    toquat[2] = axis[1] * vmod
    toquat[3] = axis[2] * vmod

cdef void cross(double * todouble, double * q1, double * q2):
    todouble[0] = q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3]
    todouble[1] = q1[0] * q2[1] + q1[1] * q2[0] + q1[2] * q2[3] - q1[3] * q2[2]
    todouble[2] = q1[0] * q2[2] - q1[1] * q2[3] + q1[2] * q2[0] + q1[3] * q2[1]
    todouble[3] = q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1] + q1[3] * q2[0]

cdef double dot(double * q1, double * q2):
    cdef double suma = 0
    cdef unsigned int i
    for i in range(4):
        suma += q1[i] * q2[i]
    return suma

cdef void div(double * todouble, double * q1, double * q2):
    cdef double qmod
    qmod = dot(q2, q2)
    todouble[0] = (q2[0] * q1[0] + q2[1] * q1[1] + q2[2] * q1[2] + q2[3] * q1[3]) / qmod
    todouble[0] = (q2[0] * q1[1] - q2[1] * q1[0] - q2[2] * q1[3] + q2[3] * q1[2]) / qmod
    todouble[0] = (q2[0] * q1[2] + q2[1] * q1[3] - q2[2] * q1[0] - q2[3] * q1[1]) / qmod
    todouble[0] = (q2[0] * q1[3] - q2[1] * q1[2] + q2[2] * q1[1] - q2[3] * q1[0]) / qmod

cdef void add(double * todouble, double * q1, double * q2):
    cdef unsigned int i
    for i in range(4):
        todouble[i] = q1[i] + q2[i]

cdef void sub(double * todouble, double * q1, double * q2):
    cdef unsigned int i
    for i in range(4):
        todouble[i] = q1[i] - q2[i]

cdef double module(double * q):
    cdef double mag = 0
    cdef unsigned int i
    for i in range(4):
        mag += q[i] * q[i]
    return sqrt(mag)

cdef void rotateVectorMatrix(double * tomat, double * q, double * center, double * mat, unsigned int numvect):
    """
    Rotate vector matrix 
    :param tomat: double pointer with the same size as the mat to be rotated
    :param q: double pointer with the quaternion
    :param center: double pointer with the center of rotation
    :param mat: double pointer with the points to be rotated (in rows)
    :param numvect: number of points to rotate
    """
    cdef double * uvv
    cdef double * uuu
    cdef unsigned int i, j
    cdef double dotuv, dotuu
    try:
        uuu = <double *> PyMem_Malloc(3 * sizeof(double))
        uvv = <double *> PyMem_Malloc(3 * sizeof(double))
        dotuu = q[0] * q[0] - AlgVector.dot(&q[1], &q[1])
        for i in range(numvect):
            dotuv = 0
            for j in range(3):
                uuu[j] = mat[j + i * 3] - center[j]
                dotuv += q[j + 1] * uuu[j]
            AlgVector.cross(uvv, &q[1], uuu)
            for j in range(3):
                tomat[j + i * 3] = 2 * dotuv * q[j + 1] + dotuu * uuu[j] + 2 * q[0] * uvv[j] + center[j]
    finally:
        PyMem_Free(uvv)
        PyMem_Free(uuu)

cpdef Quaternion fromRotationMatrix(Matrix m):
    cdef Quaternion nq = Quaternion()
    cdef double tr, S, qmod
    cdef unsigned int i
    if not m._rows == 3 and m._cols == 3:
        del nq
        raise ValueError("Matrix should be 3x3")
    tr = m._m[0] + m._m[4] + m._m[8]
    if tr > 0:
        S = sqrt(1 + tr) * 2
        nq._q[0] = 0.25 * S
        nq._q[1] = (m._m[7] - m._m[5]) / S
        nq._q[2] = (m._m[2] - m._m[6]) / S
        nq._q[3] = (m._m[3] - m._m[1]) / S
    elif m._m[0] > m._m[4] and m._m[0] > m._m[8]:
        S = sqrt(1 + m._m[0] - m._m[4] - m._m[8]) * 2
        nq._q[0] = (m._m[7] - m._m[5]) / S
        nq._q[1] = 0.25 * S
        nq._q[2] = (m._m[1] + m._m[3]) / S
        nq._q[3] = (m._m[2] + m._m[6]) / S
    elif m._m[4] > m._m[8]:
        S = sqrt(1 + m._m[4] - m._m[0] - m._m[8]) * 2
        nq._q[0] = (m._m[2] - m._m[6]) / S
        nq._q[1] = (m._m[1] + m._m[3]) / S
        nq._q[2] = 0.25 * S
        nq._q[3] = (m._m[5] + m._m[7]) / S
    else:
        S = sqrt(1 + m._m[8] - m._m[0] - m._m[4]) * 2
        nq._q[0] = (m._m[3] - m._m[1]) / S
        nq._q[1] = (m._m[2] + m._m[6]) / S
        nq._q[2] = (m._m[5] + m._m[7]) / S
        nq._q[3] = 0.25 * S
    qmod = module(nq._q)
    for i in range(4):
        nq._q[i] = nq._q[i] / qmod
    return nq

cpdef Quaternion fromEulerAngles(double R, double P, double Y):
    cdef double t0, t1, t2, t3, t4, t5
    cdef Quaternion quat = Quaternion()
    t0 = cos(Y * 0.5)
    t1 = sin(Y * 0.5)
    t2 = cos(R * 0.5)
    t3 = sin(R * 0.5)
    t4 = cos(P * 0.5)
    t5 = sin(P * 0.5)
    quat._q[0] = t0 * t2 * t4 + t1 * t3 * t5
    quat._q[1] = t0 * t3 * t4 - t1 * t2 * t5
    quat._q[2] = t0 * t2 * t5 + t1 * t3 * t4
    quat._q[3] = t1 * t2 * t4 - t0 * t3 * t5
    return quat



cdef class Quaternion():

    def __cinit__(Quaternion self, *args):
        cdef double vmod
        cdef unsigned int i
        self._q = <double *> PyMem_Malloc (4 * sizeof(double))
        if len(args) == 2:  # axis and angle
            quatFromAxisAngle(self._q, (<Vector>args[0])._v, args[1])
        elif len(args) == 4:
            for i in range(4):
                self._q[i] = args[i]
        elif len(args) == 1 and isinstance(args[0], Quaternion):
            for i in range(4):
                self._q[i] = (<Quaternion>args[0])._q[i]
        else:
            self._q[0] = 1
            self._q[1] = 0
            self._q[2] = 0
            self._q[3] = 0

    def __dealloc__(Quaternion self):
        PyMem_Free(self._q)

    def __json__(Quaternion self):
        return {'__jsoncls__': 'CLinAlg.AlgQuaternion:Quaternion.from_JSON', 'q': self.toList()}


    @classmethod
    def from_JSON(cls, jsondict):
        obj = cls(jsondict['q'])
        return obj

    #........................cmethods.............................

    cdef Vector c_getAxis(Quaternion self):
        cdef double factor
        cdef Vector axis = Vector()
        if self._q[0] == 1:
            axis._v[0] = 1
            axis._v[1] = 0
            axis._v[2] = 0
        else:
            factor = sqrt(1 - self._q[0] * self._q[0])
            axis._v[0] = self._q[1] / factor
            axis._v[1] = self._q[2] / factor
            axis._v[2] = self._q[3] / factor
        return axis

    def getAxis(Quaternion self):
        return self.c_getAxis()

    def getAngle(Quaternion self) -> double:
        return 2 * acos(self._q[0])

    cdef Matrix c_getRotationMatrix(Quaternion self):
        cdef Matrix mm = Matrix()
        mm._rows = 3
        mm._cols = 3
        mm._m = <double *> PyMem_Malloc(9 * sizeof(double))
        mm._m[0] = 1 - 2 * self._q[2] * self._q[2] - 2 * self._q[3] * self._q[3]
        mm._m[1] = 2 * self._q[1] * self._q[2] - 2 * self._q[0] * self._q[3]
        mm._m[2] = 2 * self._q[1] * self._q[3] + 2 * self._q[0] * self._q[2]
        mm._m[3] = 2 * self._q[1] * self._q[2] + 2 * self._q[0] * self._q[3]
        mm._m[4] = 1 - 2 * self._q[1] * self._q[1] - 2 * self._q[3] * self._q[3]
        mm._m[5] = 2 * self._q[2] * self._q[3] - 2 * self._q[0] * self._q[1]
        mm._m[6] = 2 * self._q[1] * self._q[3] - 2 * self._q[0] * self._q[2]
        mm._m[7] = 2 * self._q[2] * self._q[3] + 2 * self._q[0] * self._q[1]
        mm._m[8] = 1 - 2 * self._q[1] * self._q[1] - 2 * self._q[2] * self._q[2]
        return mm

    def getRotationMatrix(Quaternion self):
        return self.c_getRotationMatrix()

    cdef void c_rotateVector_p(Quaternion self, double * topointer, double * vpointer):
        "Funcion auxiliar para rotar punteros"
        cdef double * u
        cdef double * uvv
        cdef unsigned int i
        cdef double dotuv, dotuu
        try:
            u = <double *> PyMem_Malloc(3 * sizeof(double))
            uvv = <double *> PyMem_Malloc(3 * sizeof(double))
            u[0] = self._q[1]
            u[1] = self._q[2]
            u[2] = self._q[3]
            AlgVector.cross(uvv, u, vpointer)
            dotuv = AlgVector.dot(u, vpointer)
            dotuu = AlgVector.dot(u, u)
            for i in range(3):
                topointer[i] = 2 * dotuv * u[i] + (self._q[0] * self._q[0] - dotuu) * vpointer[i] + 2 * self._q[0] * \
                               uvv[i]
        finally:
            PyMem_Free(u)
            PyMem_Free(uvv)

    def rotateVector(Quaternion self, v: Vector):
        cdef Vector newVector = Vector()
        self.c_rotateVector_p(newVector._v, v._v)
        return newVector

    cdef (double, double, double) c_eulerAngles(Quaternion self):
        cdef double t0, t1, roll, pitch, yaw
        t0 = 2 * (self._q[0] * self._q[1] + self._q[2] * self._q[3])
        t1 = 1 - 2 * (self._q[1] * self._q[1] + self._q[2] * self._q[2])
        roll = atan2(t0, t1)
        t0 = 2 * (self._q[0] * self._q[2] - self._q[3] * self._q[1])
        if t0 > 1:
            t0 = 1
        elif t0 < -1:
            t0 = -1
        pitch = asin(t0)
        t0 = 2 * (self._q[0] * self._q[3] + self._q[1] * self._q[2])
        t1 = 1 - 2 * (self._q[2] * self._q[2] + self._q[3] * self._q[3])
        yaw = atan2(t0, t1)
        return roll, pitch, yaw

    def eulerAngles(Quaternion self):
        return self.c_eulerAngles()

    def module(Quaternion self):
        return module(self._q)

    def normalize(Quaternion self):
        cdef double qmod
        cdef Quaternion newQuat = Quaternion()
        cdef unsigned int i
        qmod = module(self._q)
        for i in range(4):
            newQuat._q[i] = self._q[i] / qmod
        return newQuat

    def conjugate(Quaternion self):
        cdef Quaternion newQuat = Quaternion(self._q[0], -self._q[1], -self._q[2], -self._q[3])
        return newQuat

    def inverse(Quaternion self):
        cdef Quaternion newQuat
        cdef double qmod
        qmod = module(self._q)
        qmod = qmod * qmod
        newQuat = Quaternion(self._q[0]/qmod, -self._q[1]/qmod, -self._q[2]/qmod, -self._q[3]/qmod)
        return newQuat

    def __mod__(Quaternion self, other:Quaternion):
        cdef Quaternion newQuat = Quaternion()
        cross(newQuat._q, self._q, other._q)
        return newQuat

    def __mul__(Quaternion self, other):
        cdef Quaternion newQuat = Quaternion()
        cdef unsigned int i
        if isinstance(other, Quaternion):
            del newQuat
            return dot(self._q, (<Quaternion>other)._q)
        elif isinstance(other, (int, float)):
            for i in range(4):
                newQuat._q[i] = self._q[i] * other
            return newQuat
        else:
            raise ValueError(other)

    def __truediv__(Quaternion self, other):
        cdef Quaternion newQuat = Quaternion()
        cdef unsigned int i
        if isinstance(other, Quaternion):
            div(newQuat._q, self._q, (<Quaternion>other)._q)
            return newQuat
        elif isinstance(other, (int, float)):
            for i in range(4):
                newQuat._q[i] = self._q[i] / other
            return newQuat
        else:
            del newQuat
            raise ValueError(other)

    def __add__(Quaternion self, other: Quaternion):
        cdef Quaternion newQuat = Quaternion()
        add(newQuat._q, self._q, other._q)
        return newQuat

    def __sub__(Quaternion self, Quaternion other):
        cdef Quaternion newQuat = Quaternion()
        sub(newQuat._q, self._q, other._q)
        return newQuat

    def __repr__(Quaternion self):
        return 'Quaternion' + str(self.toList())

    def __str__(Quaternion self):
        return 'Quaternion' + str(self.toList())

    def toList(Quaternion self):
        cdef unsigned int i
        return [self._q[i] for i in range(4)]

    @property
    def w(Quaternion self):
        return self._q[0]

    @property
    def x(Quaternion self):
        return self._q[1]

    @property
    def y(Quaternion self):
        return self._q[2]

    @property
    def z(Quaternion self):
        return self._q[3]

    def __getitem__(Quaternion self, key):
        return self._q[key]

    def __setitem__(Quaternion self, key, value):
        self._q[key] = value


