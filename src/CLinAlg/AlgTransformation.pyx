from cpython.mem cimport PyMem_Malloc, PyMem_Free
from .AlgQuaternion cimport Quaternion, fromRotationMatrix
from .AlgTool cimport presition
from . cimport AlgMatrix



cdef class Transformation():
    def __cinit__(self, matrix4x4=None):
        cdef unsigned int i, j
        self._mat = Matrix()
        self._mat._rows = 4
        self._mat._cols = 4
        self._mat._m = <double *> PyMem_Malloc(16 * sizeof(double))
        if isinstance(matrix4x4, list):
            for i in range(4):
                for j in range(4):
                    self._mat.c_set(i, j, matrix4x4[i][j])
        elif isinstance(matrix4x4, Transformation):
            for i in range(16):
                self._mat._m[i] = (<Transformation> matrix4x4)._mat._m[i]
        elif isinstance(matrix4x4, Matrix):
            if (<Matrix> matrix4x4)._rows != 4 or (<Matrix> matrix4x4)._cols != 4:
                raise ValueError("The matrix size should be 4x4")
            for i in range(16):
                self._mat._m[i] = (<Matrix> matrix4x4)._m[i]
        else:
            self.setIdentity()

    def __dealloc__(self):
        pass

    #...........................C.....................................
    cdef Transformation c_inverse(Transformation self):
        # create inverse of self transforamtion
        cdef Transformation newTransform
        newTransform = Transformation()
        try:
            AlgMatrix.inverse(newTransform._mat._m, self._mat._m, 4, 4)
            return newTransform
        except Exception as err:
            del(newTransform)
            raise err

    def inverse(Transformation self):
        return self.c_inverse()

    cdef Vector c_transformVector(Transformation self, Vector v):
        cdef Vector newVect = Vector()
        cdef unsigned int i, j
        cdef double val
        for i in range(3):
            val = 0
            for j in range(3):
                val += self._mat._m[j + 4 * i] * v._v[j]
            val += self._mat._m[3 + 4 * i]
            newVect._v[i] = val
        return newVect

    def transformVector(Transformation self, Vector v):
        return self.c_transformVector(v)

    cdef list c_transformVectorList(Transformation self, list vList):
        cdef Matrix matVect
        cdef unsigned int vind, i, j
        cdef size_t vlen
        cdef Vector curVect
        cdef Vector newVector
        cdef double val
        Vectores = []
        vlen = len(vList)
        for vind in range(vlen):
            curVect = vList[vind]
            newVector = Vector()
            for i in range(3):
                val = 0
                for j in range(3):
                    val += self._mat.c_get(i, j) * curVect._v[j]
                val += self._mat.c_get(i, 3)
                newVector._v[i] = val
            Vectores.append(newVector)
        return Vectores

    def transformVectorList(Transformation self, list vList):
        return self.c_transformVectorList(vList)

    #........................PYTHON...................................

    def setIdentity(self):
        cdef unsigned int i, j, k
        k = 0
        for i in range(4):
            for j in range(4):
                if i == j:
                    self._mat._m[k] = 1
                else:
                    self._mat._m[k] = 0
                k += 1

    def toRotTrans(self) -> tuple:
        cdef Quaternion quat
        cdef Vector v
        quat = fromRotationMatrix(self._mat[:3, :3])
        v = Vector()
        v._v[0] = self._mat.c_get(0, 3)
        v._v[1] = self._mat.c_get(1, 3)
        v._v[2] = self._mat.c_get(2, 3)
        return quat, v


    cpdef Transformation copy(Transformation self, bint withRotation, bint withTranslation):
        cdef Transformation new_transform
        if withRotation and withTranslation:
            new_transform = Transformation(self)
        elif withRotation:
            m = self._mat.copy()
            m.c_set(0, 3, 0)
            m.c_set(1, 3, 0)
            m.c_set(2, 3, 0)
            new_transform = Transformation(m)
            del(m)
        elif withTranslation:
            m = AlgMatrix.Matrix.identity(4)
            for i in range(3):
                m.c_set(i, 3, self._mat.c_get(i, 3))
            new_transform = Transformation(m)
            del(m)
        else:
            new_transform = Transformation()
            new_transform.setIdentity()
        return new_transform

    @property
    def mat(Transformation self):
        return self._mat

    @staticmethod
    def fromRotTrans(Quaternion quat, Vector v):
        return fromRotTrans(quat, v)

    @staticmethod
    def fromPointAndAxis(Vector point, Vector axis_x, Vector axis_y):
        return fromPointAndAxis(point, axis_x, axis_y)

    def __mul__(Transformation self, other):
        cdef Transformation newTransform
        cdef unsigned int i, j
        newTransform = Transformation()

        if isinstance(other, Transformation):
            AlgMatrix.mult(newTransform._mat._m, self._mat._m, 4, 4, (<Transformation> other)._mat._m, 4, 4)
        elif isinstance(other, Vector):
            return self.c_transformVector((<Vector> other))
        else:
            del newTransform
            raise ValueError("Multiplication for transformation not implemented")

    def __truediv__(Transformation self, Transformation other):
        cdef Transformation newTransform
        cdef double * inversed = <double *> PyMem_Malloc(16 * sizeof(double))
        try:
            if inversed == NULL:
                raise MemoryError()
            newTransform = Transformation()
            AlgMatrix.inverse(inversed, other._mat._m, 4, 4)
            AlgMatrix.mult(newTransform._mat._m, self._mat._m, 4, 4, inversed, 4, 4)
            return newTransform
        except Exception as err:
            del newTransform
            raise err
        finally:
            if inversed != NULL:
                PyMem_Free(inversed)
                inversed = NULL

    def __add__(Transformation self, Transformation other):
        return self * other

    def __sub__(Transformation self, Transformation other):
        return self / other

    def __repr__(self):
        return "Transformation: \n" + self._mat.__repr__()

    def __json__(self):
        return {'__jsoncls__': 'CLinAlg.AlgTransformation:Transformation.from_JSON', 'mat': self._mat}

    @classmethod
    def from_JSON(cls, jsondict):
        obj = cls(jsondict['mat'])
        return obj

cdef Transformation fromRotTrans(Quaternion quat, Vector vect):
    cdef Matrix m
    cdef int i, j
    cdef Transformation newTransform
    newTransform = Transformation()
    m = quat.c_getRotationMatrix()
    for i in range(3):
        for j in range(3):
            newTransform._mat.c_set(i, j, m.c_get(i, j))
    for i in range(3):
        newTransform._mat.c_set(i, 3, vect._v[i])
        newTransform._mat.c_set(3, i, 0)
    newTransform._mat.c_set(3, 3, 1)
    del m
    return newTransform

cdef Transformation fromPointAndAxis(Vector pos, Vector axis_x, Vector axis_y):
    cdef Vector axis_z = axis_x % axis_y
    if axis_z.module < presition:
        raise RuntimeError("Axis X and Axis Y should be different")
    cdef unsigned int i
    cdef Transformation newTransform = Transformation()
    for i in range(3):
        newTransform._mat._m[i] = axis_x._v[i]
        newTransform._mat._m[i+4] = axis_y._v[i]
        newTransform._mat._m[i+8] = axis_z._v[i]
        newTransform._mat._m[3+i*4] = pos._v[i]
        newTransform._mat._m[12+i] = 0
    newTransform._mat._m[15] = 1
    return newTransform