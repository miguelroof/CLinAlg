

cdef double determinant(double * mat, unsigned int rows)

cdef tuple eig(double * _mat, unsigned int nsize, bint tosorted)

cdef void transpose(double * tomat, double * fmat, unsigned int rows, unsigned int cols)

cdef void adjugate(double * tomat, double * fmat, unsigned int rows, unsigned int cols)

cdef void inverse(double * tomat, double * fmat, unsigned int rows, unsigned int cols)

cdef double threshold(unsigned int n,double * a)

cdef void protate(double * a,unsigned int nsize, double * p,unsigned int k,unsigned int l)

cdef void matAdj(double * toMat, double* fromMat, unsigned int rows, unsigned int cols, unsigned int prow, unsigned int pcol)

cdef void mult(double * tomat, double * mat1, unsigned int rows1, unsigned int cols1, double * mat2, unsigned int rows2, unsigned int cols2)

cdef void multByScalar(double * tomat, double * mat, unsigned int rows, unsigned int cols, double scalar)

cdef void add(double * tomat, double * mat1, double * mat2, unsigned int rows, unsigned int cols)

cdef void sub(double * tomat, double * mat1, double * mat2, unsigned int rows, unsigned int cols)


cdef bint hasLine(double * mat, unsigned int rows, unsigned int cols, double * vector, bint byRows)

cdef class Matrix:
    cdef unsigned int _rows, _cols
    cdef double * _m
    cdef void c_set(self, int row, int col, double data)
    cdef double c_get(self,int row, int col)
    cdef void pushdata(Matrix self, unsigned int rows,unsigned int cols, double * datalist)
    cpdef void deleteRow(Matrix self, unsigned int row)
    cpdef void deleteColumn(Matrix self, unsigned int col)
    cpdef Matrix copy(Matrix self)
    cpdef Matrix transpose(Matrix self)
    @staticmethod
    cdef c_zeros(int row, int col)
    @staticmethod
    cdef c_identity(int n)


