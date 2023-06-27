from cpython.mem cimport PyMem_Malloc, PyMem_Free
from libc.math cimport sqrt, pow, fabs
from .AlgTool cimport presition


cdef double determinant(double * mat, unsigned int rows):
    cdef double data, compose
    cdef unsigned int i,j,k, curind
    cdef double * subarray
    data = 0
    if rows == 3: # regla de sarus
        return mat[0]*mat[4]*mat[8]+mat[3]*mat[7]*mat[2]+mat[6]*mat[1]*mat[5]-mat[2]*mat[4]*mat[6]-mat[5]*mat[7]*mat[0]-mat[8]*mat[1]*mat[3]
    elif rows == 2:
        return mat[0]*mat[3]-mat[1]*mat[2]
    elif rows == 1:
        return mat[0]
    else:
        try:
            subarray = <double *> PyMem_Malloc((rows-1)*(rows-1)*sizeof(double))
            if not subarray:
                raise MemoryError()
            for i in range(rows):
                curind = 0
                for j in range(1, rows):
                    for k in range(rows):
                        if k == i:
                            continue
                        subarray[curind] = mat[k + rows*j]
                        curind += 1
                compose = determinant(subarray, rows-1)
                data += mat[i] * (pow(-1,i) * compose)
        finally:
            PyMem_Free(subarray)
    return data


cdef tuple eig(double * _mat, unsigned int nsize, bint tosorted):
    cdef unsigned int i,j,k, iterCounter, index
    cdef double * p
    cdef double * a
    cdef double * eigvalue
    cdef double mu, val
    cdef Matrix matVect
    try:
        p = <double *> PyMem_Malloc(nsize*nsize*sizeof(double)) # puntero a los autovectores
        a = <double *> PyMem_Malloc(nsize*nsize*sizeof(double))  # puntero copia de la matriz inicial, para hacer los giros
        eigvalue = <double *> PyMem_Malloc(nsize*sizeof(double)) # puntero a los elementos autovalores
        if not p or not a or not eigvalue:
            raise MemoryError()
        for i in range(nsize):
            for j in range(nsize):
                if i == j:
                    p[j+nsize*i] = 1
                else:
                    p[j+nsize*i] = 0
                a[j+nsize*i] = _mat[j+nsize*i]
        iterCounter = 0
        while iterCounter < 30:
            mu = threshold(nsize, a)
            for i in range(<unsigned int>(nsize-1)):
                for j in range(i+1, nsize):
                    if fabs(a[j+nsize*i]) >= mu:
                        protate(a,nsize, p,i,j)
            if mu <= presition*0.001:
                for i in range(nsize):
                    eigvalue[i] = a[i+nsize*i]
                break
            iterCounter += 1
        if iterCounter >= 30:
            raise ArithmeticError("No convergency for Jacobi method")

        if tosorted:
            for i in range(<unsigned int>(nsize-1)):
                index = i
                val = eigvalue[i]
                for j in range(i+1, nsize):
                    if eigvalue[j] > val:
                        index = j
                        val = eigvalue[j]
                if index != i:
                    for k in range(nsize):
                        val = p[k + i * nsize]
                        p[k + i * nsize] = p[k + index * nsize]
                        p[k + index * nsize] = val
                    val = eigvalue[i]
                    eigvalue[i] = eigvalue[index]
                    eigvalue[index] = val
        matVect = Matrix()
        matVect.pushdata(nsize, nsize, p)
        return [eigvalue[i] for i in range(nsize)], matVect
    finally:
        # print("borro")
        PyMem_Free(a)
        PyMem_Free(p)
        PyMem_Free(eigvalue)

cdef void transpose(double * tomat, double * fmat, unsigned int rows, unsigned int cols):
    cdef unsigned int i,j
    for i in range(rows):
        for j in range(cols):
            tomat[i+rows*j] = fmat[j+cols*i]

cdef void adjugate(double * tomat, double * fmat, unsigned int rows, unsigned int cols):
    # los tamanos tienen que estar equilibrados
    cdef unsigned int i,j,k, curind
    cdef double composedDet
    cdef double * subarray
    if rows == 1:
        tomat[0] = fmat[0]
        return
    subarray = <double *> PyMem_Malloc((rows-1)*(cols-1)*sizeof(double))
    if not subarray:
        raise MemoryError()
    try:
        if cols != rows:
            raise ArithmeticError("The matrix should be squared")
        curind = 0
        for i in range(rows):
            for j in range(cols):
                matAdj(subarray, fmat, rows, cols, i,j)
                composedDet = determinant(subarray, rows-1)
                tomat[curind] = composedDet*(pow(-1,(i+j)))
                curind += 1
    finally:
        PyMem_Free(subarray)

cdef void inverse(double * tomat, double * fmat, unsigned int rows, unsigned int cols):
    cdef unsigned int i
    cdef double determinante
    cdef double * tarray
    try:
        tarray = <double *> PyMem_Malloc(rows*cols*sizeof(double))
        if not tarray:
            raise MemoryError()
        if cols != rows:
            raise ArithmeticError("Is not possible to inverse not square matrix")
        determinante = determinant(fmat, rows)
        if fabs(determinante) < presition:
            raise ArithmeticError("The matrix determinant should be NOT 0")
        transpose(tarray, fmat, rows, cols)
        adjugate(tomat, tarray, rows, cols)
        for i in range(cols*rows):
            tomat[i] = tomat[i]/determinante
    finally:
        PyMem_Free(tarray)


cdef double threshold(unsigned int n,double * a):
    cdef double vsum
    cdef unsigned int i,j
    vsum = 0
    for i in range(<unsigned int>(n-1)):
        for j in range(i+1, n):
            vsum += fabs(a[j+i*n])
    return 0.5*vsum/n/(n-1)

cdef void protate(double * a,unsigned int nsize, double * p,unsigned int k,unsigned int l):
    cdef double aDiff, t, phi, c, s, tau, temp
    cdef unsigned int i

    aDiff = a[l+nsize*l] - a[k+nsize*k]
    if fabs(a[l+nsize*k]) < fabs(aDiff)*1e-20:
        t = a[l+nsize*k]/aDiff
    else:
        phi = aDiff/(2*a[l+nsize*k])
        t = 1.0/(fabs(phi) + sqrt(phi*phi+1))
        if phi < 0:
            t = -t
    c = 1/sqrt(t*t+1)
    s = t*c
    tau = s/(1+c)
    temp = a[l+nsize*k]
    a[l+nsize*k] = 0
    a[k+nsize*k] = a[k+nsize*k] - t*temp
    a[l+nsize*l] = a[l+nsize*l] + t*temp
    for i in range(k):
        temp = a[k+nsize*i]
        a[k+nsize*i] = temp - s*(a[l+nsize*i] + tau*temp)
        a[l+nsize*i] = a[l+nsize*i] + s*(temp - tau*a[l+nsize*i])
    for i in range(k+1,l):
        temp = a[i+nsize*k]
        a[i+nsize*k] = temp - s*(a[l+nsize*i]+tau*temp)
        a[l+nsize*i] = a[l+nsize*i] + s*(temp-tau*a[l+nsize*i])
    for i in range(l+1,nsize):
        temp = a[i+nsize*k]
        a[i+nsize*k] = temp - s*(a[i+nsize*l] + tau*temp)
        a[i+nsize*l] = a[i+nsize*l] + s*(temp-tau*a[i+nsize*l])
    for i in range(nsize):
        temp = p[k+nsize*i]
        p[k+nsize*i] = temp - s*(p[l+nsize*i]+ tau*p[k+nsize*i])
        p[l+nsize*i] = p[l+nsize*i] + s*(temp-tau*p[l+nsize*i])

cdef void matAdj(double * toMat, double * fromMat, unsigned int rows, unsigned int cols, unsigned int prow, unsigned int pcol):
    # newmat should be dimensioned acordly. It dont will make the verifiacation
    cdef unsigned int i,j, ipos
    ipos = 0
    for i in range(rows):
        if i == prow:
            continue
        for j in range(cols):
            if j == pcol:
                continue
            toMat[ipos] = fromMat[j+i*cols]
            ipos += 1

cdef void mult(double * tomat, double * mat1, unsigned int rows1, unsigned int cols1, double * mat2, unsigned int rows2, unsigned int cols2):
    # no hace la verificacion de tamanos !!!!!!!
    cdef unsigned int row, col, j
    cdef double value
    for row in range(rows1):
        for col in range(cols2):
            value = 0
            for j in range(cols1):
                value += mat1[j+cols1*row] * mat2[col+cols2*j]
            tomat[col + row*cols2] = value

cdef void multByScalar(double * tomat, double * mat, unsigned int rows, unsigned int cols, double scalar):
    cdef unsigned int i,j
    for i in range(rows):
        for j in range(cols):
            tomat[j+cols*i] = mat[j+cols*i] * scalar

cdef void add(double * tomat, double * mat1, double * mat2, unsigned int rows, unsigned int cols):
    # no hace la verificacion de tamanos !!!!!!!
    cdef unsigned int i, j, p
    for i in range(rows):
        for j in range(cols):
            p = j+i*cols
            tomat[p] = mat1[p] + mat2[p]

cdef void sub(double * tomat, double * mat1, double * mat2, unsigned int rows, unsigned int cols):
    # no hace la verificacion de tamanos !!!!!!!
    cdef unsigned int i, j, p
    for i in range(rows):
        for j in range(cols):
            p = j+i*cols
            tomat[p] = mat1[p] - mat2[p]

cpdef Matrix zeros(int row, int col):
    cdef int i
    cdef Matrix newMat
    newMat = Matrix()
    col = col or row
    newMat._m = <double *> PyMem_Malloc((row)*(col)*sizeof(double))
    if not newMat._m:
        raise MemoryError()
    for i in range(col*row):
        newMat._m[i] = 0.0
    newMat._rows = (<unsigned int>row)
    newMat._cols = (<unsigned int>col)
    return newMat

cpdef Matrix identity(int n):
    cdef int i,j,k
    cdef Matrix newMat
    newMat = Matrix()
    newMat._m = <double *> PyMem_Malloc(n*n*sizeof(double))
    if not newMat._m:
        raise MemoryError()
    k = 0
    for i in range(n):
        for j in range(n):
            if i == j:
                newMat._m[k] = 1.0
            else:
                newMat._m[k] = 0.0
            k += 1
    newMat._rows = <unsigned int>n
    newMat._cols = <unsigned int>n
    return newMat

cdef bint hasLine(double * mat, unsigned int rows, unsigned int cols, double * vector, bint byRows):
    """
    Function that returns if a matrix contains a vector, secuentially. dont spect great choses
    :param mat: double pointer with the matrix
    :param vector: double pointer with the values to compare
    :param numval: num of point of the vector
    :returns: True if founded, otherwhise false
    """
    cdef unsigned int i = 0
    cdef unsigned int r, c
    cdef bint isEqual
    if byRows:
        for r in range(rows):
            isEqual = True
            for c in range(cols):
                if fabs(mat[r*cols+c]-vector[c]) > presition:
                    isEqual = False
                    break
            if isEqual:
                return True
    else:
        for c in range(cols):
            isEqual = True
            for r in range(rows):
                if fabs(mat[r*cols+c]-vector[r]) > presition:
                    isEqual = False
                    break
            if isEqual:
                return True
    return False



cdef class Matrix():
    # cdef unsigned int _rows, _cols
    # cdef double * _m

    def __cinit__(Matrix self, mat=None):
        # cdef unsigned int csize
        cdef unsigned int i,j

        if isinstance(mat, (list, tuple)) and mat:
            self._rows = int(len(mat))
            self._cols = int(len(mat[0]))
            self._m = <double *> PyMem_Malloc(self._rows * self._cols * sizeof(double))
            if not self._m:
                raise MemoryError()
            for i in range(self._rows):
                for j in range(self._cols):
                    self._m[j+i*self._cols] = mat[i][j]
        elif isinstance(mat, Matrix):
            self.pushdata(mat.rows, mat.cols, (<Matrix>mat)._m)
        else:
            self._rows = 0
            self._cols = 0


    def __dealloc__(Matrix self):
        PyMem_Free(self._m)

    #.........................................C METHODS...........................................................

    cdef pushdata(Matrix self, unsigned int rows,unsigned int cols, double * datalist):
        cdef unsigned int i
        PyMem_Free(self._m)
        self._m = <double *> PyMem_Malloc(rows * cols * sizeof(double))
        if not self._m:
            raise MemoryError()
        self._rows = rows
        self._cols = cols
        for i in range(rows*cols):
            self._m[i] = datalist[i]

    cdef double c_get(self,int row, int col):
        return self._m[col + row*self._cols]

    cdef void c_set(self, int row, int col, double data):
        self._m[col + row*self._cols] = data


    cpdef void deleteColumn(Matrix self, unsigned int column):
        "funcion que elimina la columna indicada"
        cdef double * newMat
        cdef unsigned int i,j,k, numcols
        if column < 0 or column > <unsigned int>(self._cols-1):
            raise ValueError("Column should be included in interval")
        numcols = <unsigned int>(self._cols-1)
        newMat = <double *> PyMem_Malloc (self._rows*numcols*sizeof(double))
        if not newMat:
            raise MemoryError()
        k = -1
        for j in range(self._cols):
            k += 1
            if k == column:
                k -= 1
                continue
            for i in range(self._rows):
                newMat[i*numcols+k] = self.c_get(i, j)
        if self._m: #libero la memoria. borro los datos
            PyMem_Free(self._m)
        self._m = newMat
        self._cols = numcols

    cpdef void deleteRow(Matrix self, unsigned int row):
        "funcion que elimina la fila indicada"
        cdef double * newMat
        cdef unsigned int i,j,k, numrows
        if row < 0 or row > <unsigned int>(self._rows-1):
            raise ValueError("Row should be included in interval")
        numrows = <unsigned int>(self._rows-1)
        newMat = <double *> PyMem_Malloc (self._cols*numrows*sizeof(double))
        if not newMat:
            raise MemoryError()
        k = -1
        for i in range(self._rows):
            k += 1
            if i == row:
                k -= 1
                continue
            for j in range(self._cols):
                newMat[k*self._cols+j] = self.c_get(i,j)
        if self._m: #libero la memoria. borro los datos
            PyMem_Free(self._m)
        self._m = newMat
        self._rows = numrows

    #........................PYTHON METHODS.............................

    def pythonize(self):
        cdef unsigned int i,j
        mat = [[self.c_get(i,j) for j in range(self._cols)] for i in range(self._rows)]
        return mat

    def __repr__(self):
        mat = self.pythonize()
        return 'Matrix(\t' + ",\n\t\t".join([str(["%.3f" % x for x in row]) for row in mat]) + ")"

    def __str__(self):
        return self.__repr__()

    def __getitem__(Matrix self, pos):
        if isinstance(pos, tuple) and len(pos) == 2:
            if isinstance(pos[0],int) and isinstance(pos[1], int):
                return self._m[pos[1]+pos[0]*self._cols]
            elif isinstance(pos[0], slice) and isinstance(pos[1], slice):
                rangerows = range(pos[0].start or 0, pos[0].stop or self._rows, pos[0].step or 1)
                rangecols = range(pos[1].start or 0, pos[1].stop or self._cols, pos[1].step or 1)
                return Matrix([[self.c_get(i,j) for j in rangecols] for i in rangerows])
            elif isinstance(pos[0], int) and isinstance(pos[1], slice):
                rangerows = [pos[0]]
                rangecols = range(pos[1].start or 0, pos[1].stop or self._cols, pos[1].step or 1)
                return Matrix([[self.c_get(i,j) for j in rangecols] for i in rangerows])
            elif isinstance(pos[0], slice) and isinstance(pos[1], int):
                rangerows = range(pos[0].start or 0, pos[0].stop or self._rows, pos[0].step or 1)
                rangecols = [pos[1]]
                return Matrix([[self.c_get(i,j) for j in rangecols] for i in rangerows])
        elif isinstance(pos, int):
            if self._rows == 1:
                return self._m[pos]
            else:
                return Matrix([[self._m[pos*self._cols+j] for j in range(self._cols)]])

    def __setitem__(Matrix self, pos:tuple, value):
        self.c_set(pos[0],pos[1],value)


    def __json__(self):
        return {'__jsoncls__': 'CLinAlg.AlgMatrix:Matrix.from_JSON', 'm': self.pythonize()}

    @classmethod
    def from_JSON(cls, jsondict):
        obj = cls(jsondict['m'])
        return obj

    cpdef Matrix copy(Matrix self):
        return Matrix(self)

    def size(Matrix self) -> tuple:
        return self._rows, self._cols

    def eigMPmath(self):
        import mpmath
        # print("AlgMatrix.Matrix cannot solve eigenvalue. Use mpmath")
        n = self.rows
        eigvalue, rot = mpmath.eig(mpmath.matrix(self.pythonize()))
        eigvalue = mpmath.diag(eigvalue)
        eigvalue = [float(eigvalue[i, i]) for i in range(n)]
        eigvect = [[float(x) for x in row] for row in rot.tolist()]
        return eigvalue, Matrix(eigvect)

    def eig(self, tosorted=True):
        cdef unsigned int i,j,k, nsize, index
        cdef Matrix eigvect
        cdef double val

        if not self.isSimetric():
            eigvalue, eigvect = self.eigMPmath()
        else:
            try:
                eigvalue, eigvect = eig(self._m, self._rows, tosorted)
            except:
                eigvalue, eigvect = self.eigMPmath()

        if tosorted:
            nsize = self._rows
            for i in range(<unsigned int>(nsize-1)):
                index = i
                val = eigvalue[i]
                for j in range(i + 1, nsize):
                    if eigvalue[j] > val:
                        index = j
                        val = eigvalue[j]
                if index != i:
                    # tengo que permutar las lineas i, index
                    for k in range(nsize):
                        val = (<Matrix>eigvect)._m[k + i * nsize]
                        (<Matrix>eigvect)._m[k + i * nsize] = (<Matrix>eigvect)._m[k + index * nsize]
                        (<Matrix>eigvect)._m[k + index * nsize] = val
                    val = eigvalue[i]
                    eigvalue[i] = eigvalue[index]
                    eigvalue[index] = val
        return eigvalue, eigvect

    #..............................operadores...............................................
    def __add__(Matrix self, other):
        cdef unsigned int i
        cdef Matrix newMat
        newMat = Matrix()

        if isinstance(other, Matrix):
            if (<Matrix>other)._cols != self._cols or (<Matrix>other)._rows != self._rows:
                del newMat
                raise ValueError("Dimmensions shoud be compatible")
            newMat._m = <double *> PyMem_Malloc(self._rows * self._cols * sizeof(double))
            if not newMat._m:
                raise MemoryError()
            add(newMat._m, self._m, (<Matrix>other)._m, self._rows, self._cols)
        elif isinstance(other, (int, float)):
            newMat._m = <double *> PyMem_Malloc(self._rows * self._cols * sizeof(double))
            if not newMat._m:
                raise MemoryError()
            for i in range(self._cols * self._rows):
                newMat._m[i] = self._m[i] + other
        else:
            raise ValueError("Not implemented add")
        newMat._rows = self._rows
        newMat._cols = self._cols
        return newMat

    def __sub__(Matrix self, other):
        cdef unsigned int i
        cdef Matrix newMat
        newMat = Matrix()
        if isinstance(other, Matrix):
            if (<Matrix>other)._rows != self._rows or (<Matrix>other)._cols != self._cols:
                del newMat
                raise ValueError("Dimmensions shoud be compatible")
            newMat._m = <double *> PyMem_Malloc(self._rows * self._cols * sizeof(double))
            if not newMat._m:
                raise MemoryError()
            sub(newMat._m, self._m, (<Matrix>other)._m, self._rows, self._cols)
        elif isinstance(other, (int, float)):
            for i in range(self._cols * self._rows):
                newMat._m[i] = self._m[i] - other
        else:
            raise ValueError("Not implemented add")
        newMat._rows = self._rows
        newMat._cols = self._cols
        return newMat

    def __serialize__(self) -> list:
        cdef unsigned int i
        return [self._m[i] for i in range(self._rows*self._cols)]

    def __eq__(Matrix self, Matrix other) -> bool:
        cdef unsigned int i
        if self.size() != other.size():
            return False
        return all([(fabs(self._m[i]-other._m[i]) < presition) for i in range(self._cols*self._rows)])

    def __ne__(Matrix self, Matrix other) -> bool:
        return not self.__eq__(other)

    def __neg__(Matrix self):
        cdef unsigned int i
        cdef Matrix newMat
        newMat = Matrix()
        newMat._m = <double *> PyMem_Malloc(self._rows*self._cols*sizeof(double))
        if not newMat._m:
            raise MemoryError()
        newMat._rows = self._rows
        newMat._cols = self._cols
        for i in range(self._cols*self._rows):
            newMat._m[i] = self._m[i]*-1
        return newMat

    def __mul__(first, other):
        cdef unsigned int row, col, j
        cdef double value
        cdef Matrix newMat
        newMat = Matrix()
        if isinstance(other, Matrix) and isinstance(first, Matrix):
            if first.cols != other.rows:
                raise ArithmeticError("The matrix dimmensions should be compatible")
            newMat._m = <double *> PyMem_Malloc(other.cols * first.rows * sizeof(double))
            if not newMat._m:
                raise MemoryError()
            newMat._rows = (<Matrix>first)._rows
            newMat._cols = (<Matrix>other)._cols
            mult(newMat._m, (<Matrix>first)._m, first.rows, first.cols,
                                (<Matrix>other)._m, other.rows, other.cols)
        elif isinstance(other, (int, float)) and isinstance(first, Matrix):
            newMat._m = <double *> PyMem_Malloc(first.cols * first.rows * sizeof(double))
            if not newMat._m:
                raise MemoryError()
            newMat._rows = (<Matrix>first)._rows
            newMat._cols = (<Matrix>first)._cols
            multByScalar(newMat._m, (<Matrix>first)._m, first.rows, first.cols, other)
        elif isinstance(other,Matrix) and isinstance(first, (int, float)):
            newMat._m = <double *> PyMem_Malloc(other.cols * other.rows * sizeof(double))
            if not newMat._m:
                raise MemoryError()
            newMat._rows = (<Matrix>other)._rows
            newMat._cols = (<Matrix>other)._cols
            multByScalar(newMat._m, (<Matrix>other)._m, other.rows, other.cols, (<double>first))
        else:
            raise ValueError("Not implemented add")
        # newMat.pushdata(self._rows, self._cols, mm)
        return newMat


    def __bool__(Matrix self) -> bool:
        cdef unsigned int i
        return any([fabs(self._m[i]) > presition for i in range(self._rows*self._cols)])

    def __div__(Matrix self, other):
        cdef unsigned int i
        cdef Matrix newMat

        if isinstance(other, Matrix):
            newMat = (<Matrix>other).inverse()
            return self.__mul__(newMat)
        elif isinstance(other, (int,float)):
            if other == 0:
                raise ValueError("Cannot be divided by 0")
            newMat = Matrix()
            newMat._m = <double *>PyMem_Malloc(self._rows*self._cols*sizeof(double))
            if not newMat._m:
                raise MemoryError()
            for i in range(self._rows*self._cols):
                newMat._m[i] = self._m[i]/other
            return newMat

    def __truediv__(Matrix self, other):
        return self.__div__(other)

    def __copy__(Matrix self) -> Matrix:
        return self.copy()

    def __hash__(Matrix self):
        return hash(tuple(self.__serialize__()))


    @property
    def rows(Matrix self):
        return self._rows

    @rows.setter
    def rows(Matrix self, unsigned int rows):
        cdef double * mat
        cdef unsigned int i,j,k
        mat = <double *>PyMem_Malloc(rows*self._cols*sizeof(double))
        if not mat:
            raise MemoryError()
        k = 0
        for i in range(rows):
            for j in range(self._cols):
                if i < self._rows:
                    mat[k] = self.c_get(i,j)
                else:
                    mat[k] = 0
                k += 1
        self.pushdata(rows, self._cols, mat)
        PyMem_Free(mat)

    @property
    def cols(Matrix self):
        return self._cols

    @cols.setter
    def cols(Matrix self, unsigned int cols):
        cdef double * mat
        cdef unsigned int i, j, k
        mat = <double *> PyMem_Malloc(self._rows * cols * sizeof(double))
        if not mat:
            raise MemoryError()
        k = 0
        for i in range(self._rows):
            for j in range(cols):
                if j < self._cols:
                    mat[k] = self.c_get(i, j)
                else:
                    mat[k] = 0
                k += 1
        self.pushdata(self._rows, cols, mat)
        PyMem_Free(mat)

    def isSimetric(Matrix self):
        cdef unsigned int row, col, rowm
        if self._cols != self._rows:
            return False
        rowm = self._cols-1
        for row in range(rowm):
            for col in range(row+1, self._cols):
                if fabs(self.c_get(row,col)-self.c_get(col,row)) > presition:
                    return False
        return True


    cpdef Matrix transpose(Matrix self):
        cdef Matrix newMat = Matrix()
        newMat._rows = self._cols
        newMat._cols = self._rows
        newMat._m = <double *> PyMem_Malloc(self._rows*self._cols*sizeof(double))
        if not newMat._m:
            raise MemoryError()
        transpose(newMat._m, self._m, self._rows, self._cols)
        return newMat


    def adjugate(Matrix self):
        cdef Matrix newMat
        if self._rows != self._cols:
            raise ArithmeticError("The matrix should be squared for adjugate calculation")
        newMat = Matrix()
        newMat._rows = self._rows
        newMat._cols = self._cols
        newMat._m = <double *> PyMem_Malloc(self._rows*self._cols*sizeof(double))
        if not newMat._m:
            raise MemoryError()
        adjugate(newMat._m, self._m, self._rows, self._cols)
        return newMat

    def determinant(Matrix self) -> double:
        return determinant(self._m, self._rows)

    def inverse(Matrix self):
        cdef Matrix newMat
        if self._rows != self._cols:
            raise ArithmeticError("The matrix should be squared for adjugate calculation")
        newMat = Matrix()
        newMat._rows = self._rows
        newMat._cols = self._cols
        newMat._m = <double *> PyMem_Malloc(self._rows*self._cols*sizeof(double))
        if not newMat._m:
            raise MemoryError()
        inverse(newMat._m, self._m, self._rows, self._cols)
        return newMat






