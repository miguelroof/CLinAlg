from cpython.mem cimport PyMem_Malloc, PyMem_Free, PyMem_Realloc
from libc.math cimport sqrt, pow, fabs
from .AlgTool cimport presition



cdef double determinant(double * mat, unsigned int rows):
    """
    Calculate the determinant of a square matrix using recursive expansion.
    
    Uses Sarrus rule for 3x3 matrices and direct calculation for 2x2 and 1x1.
    For larger matrices, uses cofactor expansion along the first row.
    
    Args:
        mat: Pointer to matrix data in row-major order
        rows: Number of rows (and columns) in the square matrix
    
    Returns:
        The determinant value as a double
    
    Raises:
        MemoryError: If memory allocation fails for subarray
    """
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
        subarray = <double *> PyMem_Malloc((rows - 1) * (rows - 1) * sizeof(double))

        try:
            if not subarray:
                raise MemoryError("No subarray in matrix module")
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
            if subarray != NULL:
                PyMem_Free(subarray)
                subarray = NULL
    return data


cdef tuple eig(double * _mat, unsigned int nsize, bint tosorted):
    """
    Calculate eigenvalues and eigenvectors using the Jacobi method.
    
    This method is suitable for symmetric matrices and uses iterative rotations
    to diagonalize the matrix.
    
    Args:
        _mat: Pointer to matrix data in row-major order
        nsize: Size of the square matrix (number of rows/columns)
        tosorted: If True, sort eigenvalues in descending order
    
    Returns:
        A tuple containing:
        - List of eigenvalues
        - Matrix object containing eigenvectors as columns
    
    Raises:
        MemoryError: If memory allocation fails
        ArithmeticError: If Jacobi method doesn't converge within 30 iterations
    """
    cdef unsigned int i,j,k, iterCounter, index
    cdef double * p
    cdef double * a
    cdef double * eigvalue
    cdef double mu, val
    cdef Matrix matVect
    p = <double *> PyMem_Malloc(nsize * nsize * sizeof(double))  # puntero a los autovectores
    a = <double *> PyMem_Malloc(
        nsize * nsize * sizeof(double))  # puntero copia de la matriz inicial, para hacer los giros
    eigvalue = <double *> PyMem_Malloc(nsize * sizeof(double))  # puntero a los elementos autovalores

    try:
        if not p or not a or not eigvalue:
            raise MemoryError("Error on eigenvalues in matrix module")
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
        if a != NULL:
            PyMem_Free(a)
            a = NULL
        if p != NULL:
            PyMem_Free(p)
            p = NULL
        if eigvalue != NULL:
            PyMem_Free(eigvalue)
            eigvalue = NULL

cdef void transpose(double * tomat, double * fmat, unsigned int rows, unsigned int cols):
    """
    Transpose a matrix by swapping rows and columns.
    
    Args:
        tomat: Pointer to destination matrix (should be allocated as cols x rows)
        fmat: Pointer to source matrix (rows x cols)
        rows: Number of rows in source matrix
        cols: Number of columns in source matrix
    """
    cdef unsigned int i,j
    for i in range(rows):
        for j in range(cols):
            tomat[i+rows*j] = fmat[j+cols*i]

cdef void adjugate(double * tomat, double * fmat, unsigned int rows, unsigned int cols):
    """
    Calculate the adjugate (adjoint) matrix.
    
    The adjugate is the transpose of the cofactor matrix. Each element is the
    determinant of the minor matrix multiplied by the appropriate sign.
    
    Args:
        tomat: Pointer to destination matrix for adjugate result
        fmat: Pointer to source matrix
        rows: Number of rows in the matrix
        cols: Number of columns in the matrix
    
    Raises:
        ArithmeticError: If the matrix is not square
        MemoryError: If memory allocation fails for subarray
    """
    # los tamanos tienen que estar equilibrados
    cdef unsigned int i,j,k, curind
    cdef double composedDet
    cdef double * subarray
    if rows == 1:
        tomat[0] = fmat[0]
        subarray = NULL
        return
    try:
        if cols != rows:
            raise ArithmeticError("The matrix should be squared")
        subarray = <double *> PyMem_Malloc((rows - 1) * (cols - 1) * sizeof(double))
        if not subarray:
            raise MemoryError("No subarray in matrix module")

        curind = 0
        for i in range(rows):
            for j in range(cols):
                matAdj(subarray, fmat, rows, cols, i,j)
                composedDet = determinant(subarray, rows-1)
                tomat[curind] = composedDet*(pow(-1,(i+j)))
                curind += 1
    except Exception as err:
        raise err
    finally:
        if subarray != NULL:
            PyMem_Free(subarray)
            subarray = NULL

cdef void inverse(double * tomat, double * fmat, unsigned int rows, unsigned int cols):
    """
    Calculate the inverse of a square matrix.
    
    Uses the formula: A^(-1) = (1/det(A)) * adj(A)^T
    where adj(A) is the adjugate matrix.
    
    Args:
        tomat: Pointer to destination matrix for inverse result
        fmat: Pointer to source matrix
        rows: Number of rows in the matrix
        cols: Number of columns in the matrix
    
    Raises:
        ArithmeticError: If matrix is not square or determinant is zero
        MemoryError: If memory allocation fails
    """
    cdef unsigned int i
    cdef double determinante
    cdef double * tarray
    try:
        if cols != rows:
            raise ArithmeticError("Is not possible to inverse not square matrix")
        tarray = <double *> PyMem_Malloc(rows*cols*sizeof(double))
        if not tarray:
            raise MemoryError("No tarray in inverse funciton at Matrix module")
        determinante = determinant(fmat, rows)
        if fabs(determinante) < presition:
            raise ArithmeticError("The matrix determinant should be NOT 0")
        transpose(tarray, fmat, rows, cols)
        adjugate(tomat, tarray, rows, cols)
        for i in range(cols*rows):
            tomat[i] = tomat[i]/determinante
    finally:
        if tarray != NULL:
            PyMem_Free(tarray)
            tarray = NULL


cdef double threshold(unsigned int n,double * a):
    """
    Calculate the threshold value for Jacobi iteration convergence.
    
    Computes the average of absolute values of off-diagonal elements
    in the upper triangle of the matrix.
    
    Args:
        n: Size of the square matrix
        a: Pointer to matrix data
    
    Returns:
        Threshold value used to determine convergence
    """
    cdef double vsum
    cdef unsigned int i,j
    vsum = 0
    for i in range(<unsigned int>(n-1)):
        for j in range(i+1, n):
            vsum += fabs(a[j+i*n])
    return 0.5*vsum/n/(n-1)

cdef void protate(double * a,unsigned int nsize, double * p,unsigned int k,unsigned int l):
    """
    Perform a Jacobi rotation to zero out element a[l,k].
    
    This is a key step in the Jacobi eigenvalue algorithm. It applies a
    rotation to both the matrix and the eigenvector accumulator.
    
    Args:
        a: Pointer to the matrix being diagonalized (modified in place)
        nsize: Size of the square matrix
        p: Pointer to eigenvector matrix (modified in place)
        k: First index for rotation
        l: Second index for rotation (l > k)
    """
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
    """
    Extract a minor matrix by removing specified row and column.
    
    Creates a submatrix by excluding the specified row and column from the
    source matrix. Used in determinant and adjugate calculations.
    
    Args:
        toMat: Pointer to destination matrix (should be (rows-1) x (cols-1))
        fromMat: Pointer to source matrix
        rows: Number of rows in source matrix
        cols: Number of columns in source matrix
        prow: Row index to exclude
        pcol: Column index to exclude
    
    Note:
        Destination matrix should be pre-allocated with correct size.
    """
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
    """
    Multiply two matrices: tomat = mat1 * mat2.
    
    Performs standard matrix multiplication. Assumes dimensions are compatible
    (cols1 == rows2).
    
    Args:
        tomat: Pointer to destination matrix (rows1 x cols2)
        mat1: Pointer to first matrix (rows1 x cols1)
        rows1: Number of rows in first matrix
        cols1: Number of columns in first matrix
        mat2: Pointer to second matrix (rows2 x cols2)
        rows2: Number of rows in second matrix
        cols2: Number of columns in second matrix
    
    Note:
        Does not verify dimension compatibility - caller must ensure cols1 == rows2.
    """
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
    """
    Multiply a matrix by a scalar value.
    
    Args:
        tomat: Pointer to destination matrix
        mat: Pointer to source matrix
        rows: Number of rows in the matrix
        cols: Number of columns in the matrix
        scalar: Scalar value to multiply each element by
    """
    cdef unsigned int i,j
    for i in range(rows):
        for j in range(cols):
            tomat[j+cols*i] = mat[j+cols*i] * scalar

cdef void add(double * tomat, double * mat1, double * mat2, unsigned int rows, unsigned int cols):
    """
    Add two matrices element-wise: tomat = mat1 + mat2.
    
    Args:
        tomat: Pointer to destination matrix
        mat1: Pointer to first matrix
        mat2: Pointer to second matrix
        rows: Number of rows in the matrices
        cols: Number of columns in the matrices
    
    Note:
        Does not verify dimension compatibility - caller must ensure same dimensions.
    """
    # no hace la verificacion de tamanos !!!!!!!
    cdef unsigned int i, j, p
    for i in range(rows):
        for j in range(cols):
            p = j+i*cols
            tomat[p] = mat1[p] + mat2[p]

cdef void sub(double * tomat, double * mat1, double * mat2, unsigned int rows, unsigned int cols):
    """
    Subtract two matrices element-wise: tomat = mat1 - mat2.
    
    Args:
        tomat: Pointer to destination matrix
        mat1: Pointer to first matrix
        mat2: Pointer to second matrix
        rows: Number of rows in the matrices
        cols: Number of columns in the matrices
    
    Note:
        Does not verify dimension compatibility - caller must ensure same dimensions.
    """
    # no hace la verificacion de tamanos !!!!!!!
    cdef unsigned int i, j, p
    for i in range(rows):
        for j in range(cols):
            p = j+i*cols
            tomat[p] = mat1[p] - mat2[p]


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
    """
    A high-performance matrix class implemented in Cython.

    This class provides efficient matrix operations using C-level memory management
    and optimized algorithms. Supports basic arithmetic operations, linear algebra
    operations (determinant, inverse, eigenvalues), and matrix manipulations.

    Attributes:
        _rows (unsigned int): Number of rows in the matrix
        _cols (unsigned int): Number of columns in the matrix
        _m (double*): Pointer to matrix data stored in row-major order

    Examples:
        >>> m = Matrix([[1, 2], [3, 4]])
        >>> m.determinant()
        -2.0
        >>> m2 = Matrix.zeros(3, 3)
        >>> m3 = Matrix.identity(4)
    """
    # cdef unsigned int _rows, _cols
    # cdef double * _m

    def __cinit__(Matrix self, mat=None):
        """
        Initialize a Matrix object.

        Args:
            mat: Optional initialization data. Can be:
                - A list or tuple of lists (2D array)
                - Another Matrix object (creates a copy)
                - None (creates an empty matrix)

        Raises:
            MemoryError: If memory allocation fails

        Examples:
            >>> m1 = Matrix([[1, 2], [3, 4]])
            >>> m2 = Matrix(m1)  # Copy constructor
            >>> m3 = Matrix()    # Empty matrix
        """
        # cdef unsigned int csize
        cdef unsigned int i,j

        if isinstance(mat, (list, tuple)) and mat:
            self._rows = int(len(mat))
            self._cols = int(len(mat[0]))
            self._m = <double *> PyMem_Malloc(self._rows * self._cols * sizeof(double))
            if not self._m:
                raise MemoryError("No Matrix initialization")
            for i in range(self._rows):
                for j in range(self._cols):
                    self._m[j+i*self._cols] = mat[i][j]
        elif isinstance(mat, Matrix):
            self.pushdata(mat.rows, mat.cols, (<Matrix>mat)._m)
        else:
            self._rows = 0
            self._cols = 0


    def __dealloc__(Matrix self):
        """
        Deallocate matrix memory when the object is destroyed.

        This method is automatically called by Python's garbage collector.
        Frees the C-level memory allocated for matrix data.
        """
        if self._m:
            # print("BORRADO DE MATRIZ %i x %i" % (self.rows, self.cols))
            PyMem_Free(self._m)
            self._m = NULL

    #.........................................C METHODS...........................................................

    cdef void pushdata(Matrix self, unsigned int rows,unsigned int cols, double * datalist):
        """
        Replace matrix data with new data from a C array.
        
        This is a C-level method that frees existing matrix memory and allocates
        new memory to store the provided data.
        
        Args:
            rows: Number of rows in the new matrix
            cols: Number of columns in the new matrix
            datalist: Pointer to data array in row-major order
        
        Raises:
            MemoryError: If memory allocation fails
        """
        cdef unsigned int i
        if self._m != NULL:
            PyMem_Free(self._m)
        self._m = NULL
        self._m = <double *> PyMem_Malloc(rows * cols * sizeof(double))
        if not self._m:
            raise MemoryError("Error at alocating memory")
        self._rows = rows
        self._cols = cols
        for i in range(rows*cols):
            self._m[i] = datalist[i]

    cdef double c_get(self,int row, int col):
        """
        Get a matrix element at the specified position (C-level method).
        
        Args:
            row: Row index (0-based)
            col: Column index (0-based)
        
        Returns:
            The value at the specified position
        """
        return self._m[col + row*self._cols]

    cdef void c_set(self, int row, int col, double data):
        """
        Set a matrix element at the specified position (C-level method).
        
        Args:
            row: Row index (0-based)
            col: Column index (0-based)
            data: Value to set at the specified position
        """
        self._m[col + row*self._cols] = data


    cpdef void deleteColumn(Matrix self, unsigned int column):
        """
        Delete a column from the matrix.
        
        Removes the specified column and resizes the matrix. The matrix is modified
        in place with a new memory allocation.
        
        Args:
            column: Index of the column to delete (0-based)
        
        Raises:
            ValueError: If column index is out of range or matrix has only 1 column
            MemoryError: If memory allocation fails
        
        Examples:
            >>> m = Matrix([[1, 2, 3], [4, 5, 6]])
            >>> m.deleteColumn(1)
            >>> # m is now [[1, 3], [4, 6]]
        """
        cdef double * newMat
        cdef unsigned int i,j,k, numcols
        if column >= self._cols:
            raise ValueError("Column should be included in interval")
        if self._cols == 1:
            raise ValueError("Cannot delete column from 1-column matrix")
        numcols = <unsigned int>(self._cols-1)
        newMat = <double *> PyMem_Malloc (self._rows*numcols*sizeof(double))
        if not newMat:
            raise MemoryError()
        k = 0
        for j in range(self._cols):
            if j == column:
                continue
            for i in range(self._rows):
                newMat[i*numcols+k] = self.c_get(i, j)
            k += 1
        if self._m: #libero la memoria. borro los datos
            PyMem_Free(self._m)
            self._m = NULL
        self._m = newMat
        self._cols = numcols

    cpdef void deleteRow(Matrix self, unsigned int row):
        """
        Delete a row from the matrix.
        
        Removes the specified row and resizes the matrix. The matrix is modified
        in place with a new memory allocation.
        
        Args:
            row: Index of the row to delete (0-based)
        
        Raises:
            ValueError: If row index is out of range or matrix has only 1 row
            MemoryError: If memory allocation fails
        
        Examples:
            >>> m = Matrix([[1, 2], [3, 4], [5, 6]])
            >>> m.deleteRow(1)
            >>> # m is now [[1, 2], [5, 6]]
        """
        cdef double * newMat
        cdef unsigned int i,j,k, numrows
        if row >= self._rows:
            raise ValueError("Row should be included in interval")
        if self._rows == 1:
            raise ValueError("Cannot delete row from 1-row matrix")
        numrows = <unsigned int>(self._rows-1)
        newMat = <double *> PyMem_Malloc (self._cols*numrows*sizeof(double))
        if not newMat:
            raise MemoryError()
        k = 0
        for i in range(self._rows):
            if i == row:
                continue
            for j in range(self._cols):
                newMat[k*self._cols+j] = self.c_get(i,j)
            k += 1
        if self._m: #libero la memoria. borro los datos
            PyMem_Free(self._m)
            self._m = NULL
        self._m = newMat
        self._rows = numrows

    #........................PYTHON METHODS.............................


    def pythonized(self):
        """
        Convert the matrix to a Python list of lists.

        Returns:
            A 2D list representation of the matrix

        Examples:
            >>> m = Matrix([[1, 2], [3, 4]])
            >>> m.pythonized()
            [[1.0, 2.0], [3.0, 4.0]]
        """
        cdef unsigned int i,j
        mat = [[self.c_get(i,j) for j in range(self._cols)] for i in range(self._rows)]
        return mat

    def toList(self):
        """
        Convert the matrix to a Python list of lists.

        This is an alias for pythonized().

        Returns:
            A 2D list representation of the matrix
        """
        return self.pythonized()


    def __repr__(self):
        """
        Return a string representation of the matrix.

        Returns:
            A formatted string showing matrix contents
        """
        mat = self.pythonized()
        return 'Matrix(\t' + ",\n\t\t".join([str(["%.3f" % x for x in row]) for row in mat]) + ")"

    def __str__(self):
        """
        Return a string representation of the matrix.

        Returns:
            A formatted string showing matrix contents
        """
        return self.__repr__()

    def __getitem__(Matrix self, pos):
        """
        Get matrix element(s) using indexing and slicing notation.

        Supports various indexing modes:
        - Single element: m[i, j]
        - Row slicing: m[i, :] or m[i:j, :]
        - Column slicing: m[:, j] or m[:, i:j]
        - Submatrix: m[i:j, k:l]
        - Row access: m[i] (returns row as 1xN matrix)

        Args:
            pos: Index specification (int, tuple, or slice)

        Returns:
            Either a scalar value or a new Matrix object

        Raises:
            IndexError: If indices are out of range

        Examples:
            >>> m = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
            >>> m[0, 1]  # 2.0
            >>> m[0, :]  # Matrix([[1, 2, 3]])
            >>> m[:, 1]  # Matrix([[2], [5], [8]])
            >>> m[0:2, 1:3]  # Matrix([[2, 3], [5, 6]])
        """
        if isinstance(pos, tuple) and len(pos) == 2:
            if isinstance(pos[0],int) and isinstance(pos[1], int):
                if pos[0] >= self._rows or pos[1] >= self._cols:
                    raise IndexError("Matrix size are incompatible")
                return self._m[pos[1]+pos[0]*self._cols]
            elif isinstance(pos[0], slice) and isinstance(pos[1], slice):
                rangerows = range(pos[0].start or 0, pos[0].stop or self._rows, pos[0].step or 1)
                rangecols = range(pos[1].start or 0, pos[1].stop or self._cols, pos[1].step or 1)
                if rangerows.start < 0 or rangerows.start > self._rows or rangerows.stop < 0 or rangerows.stop > self._rows:
                    raise IndexError("Matrix size are incompatible")
                if rangecols.start < 0 or rangecols.start > self._cols or rangecols.stop < 0 or rangecols.stop > self._cols:
                    raise IndexError("Matrix size are incompatible")
                return Matrix([[self.c_get(i,j) for j in rangecols] for i in rangerows])
            elif isinstance(pos[0], int) and isinstance(pos[1], slice):
                if pos[0] >= self._rows or pos[0] < 0:
                    raise IndexError("Matrix size are incompatible")
                rangerows = [pos[0]]
                rangecols = range(pos[1].start or 0, pos[1].stop or self._cols, pos[1].step or 1)
                if rangecols.start < 0 or rangecols.start > self._cols or rangecols.stop < 0 or rangecols.stop > self._cols:
                    raise IndexError("Matrix size are incompatible")
                return Matrix([[self.c_get(i,j) for j in rangecols] for i in rangerows])
            elif isinstance(pos[0], slice) and isinstance(pos[1], int):
                rangerows = range(pos[0].start or 0, pos[0].stop or self._rows, pos[0].step or 1)
                if rangerows.start < 0 or rangerows.start > self._rows or rangerows.stop < 0 or rangerows.stop > self._rows:
                    raise IndexError("Matrix size are incompatible")
                if pos[1] >= self._cols or pos[1] < 0:
                    raise IndexError("Matrix size are incompatible")
                rangecols = [pos[1]]
                return Matrix([[self.c_get(i,j) for j in rangecols] for i in rangerows])
        elif isinstance(pos, int):
            if self._rows == 1:
                if pos >= self._cols:
                    raise IndexError("Matrix size are incompatible")
                return self._m[pos]
            else:
                if pos >= self._rows:
                    raise IndexError("Matrix size are incompatible")
                return Matrix([[self._m[pos*self._cols+j] for j in range(self._cols)]])

    def __setitem__(Matrix self, pos:tuple, value):
        """
        Set a matrix element using indexing notation.

        Args:
            pos: A tuple (row, col) specifying the position
            value: The value to set

        Examples:
            >>> m = Matrix([[1, 2], [3, 4]])
            >>> m[0, 1] = 99
        """
        self.c_set(pos[0],pos[1],value)


    def __json__(self):
        """
        Serialize matrix to JSON-compatible format.

        Returns:
            A dictionary containing matrix data and class information
        """
        return {'__jsoncls__': 'CLinAlg.AlgMatrix:Matrix.from_JSON', 'm': self.pythonized()}

    @classmethod
    def from_JSON(cls, jsondict):
        """
        Deserialize matrix from JSON format.

        Args:
            jsondict: Dictionary containing matrix data

        Returns:
            A new Matrix object
        """
        obj = cls(jsondict['m'])
        return obj

    cpdef Matrix copy(Matrix self):
        """
        Create a deep copy of the matrix.
        
        Returns:
            A new Matrix object with the same data
        
        Examples:
            >>> m1 = Matrix([[1, 2], [3, 4]])
            >>> m2 = m1.copy()
            >>> m2[0, 0] = 99
            >>> m1[0, 0]  # Original unchanged
            1.0
        """
        return Matrix(self)

    def size(Matrix self) -> tuple:
        """
        Get the dimensions of the matrix.

        Returns:
            A tuple (rows, cols) representing the matrix dimensions

        Examples:
            >>> m = Matrix([[1, 2, 3], [4, 5, 6]])
            >>> m.size()
            (2, 3)
        """
        return self._rows, self._cols

    def eigMPmath(self):
        """
        Calculate eigenvalues and eigenvectors using mpmath library.

        This method is used as a fallback for non-symmetric matrices or when
        the Jacobi method fails to converge.

        Returns:
            A tuple containing:
            - List of eigenvalues
            - Matrix object containing eigenvectors as columns

        Note:
            Requires the mpmath library to be installed.
        """
        import mpmath
        # print("AlgMatrix.Matrix cannot solve eigenvalue. Use mpmath")
        n = self.rows
        eigvalue, rot = mpmath.eig(mpmath.matrix(self.pythonized()))
        eigvalue = mpmath.diag(eigvalue)
        eigvalue = [float(eigvalue[i, i]) for i in range(n)]
        eigvect = [[float(x) for x in row] for row in rot.tolist()]
        return eigvalue, Matrix(eigvect)

    def eig(self, bint tosorted=True):
        """
        Calculate eigenvalues and eigenvectors of the matrix.

        Uses the Jacobi method for symmetric matrices, falls back to mpmath
        for non-symmetric matrices or if Jacobi fails to converge.

        Args:
            tosorted: If True, sort eigenvalues in descending order (default: True)

        Returns:
            A tuple containing:
            - List of eigenvalues
            - Matrix object containing eigenvectors as columns

        Examples:
            >>> m = Matrix([[4, 1], [1, 3]])
            >>> eigenvalues, eigenvectors = m.eig()

        Note:
            For non-symmetric matrices, requires mpmath library.
        """
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
        """
        Add two matrices or add a scalar to all matrix elements.

        Args:
            other: Either a Matrix object or a scalar (int/float)

        Returns:
            A new Matrix containing the result

        Raises:
            ValueError: If matrix dimensions are incompatible or operation not supported

        Examples:
            >>> m1 = Matrix([[1, 2], [3, 4]])
            >>> m2 = Matrix([[5, 6], [7, 8]])
            >>> m3 = m1 + m2  # [[6, 8], [10, 12]]
            >>> m4 = m1 + 10  # [[11, 12], [13, 14]]
        """
        cdef unsigned int i
        cdef Matrix newMat
        if isinstance(other, Matrix):
            if (<Matrix>other)._cols != self._cols or (<Matrix>other)._rows != self._rows:
                raise ValueError("Dimmensions shoud be compatible")
            newMat = Matrix()
            newMat._m = <double *> PyMem_Malloc(self._rows * self._cols * sizeof(double))
            if not newMat._m:
                raise MemoryError()
            add(newMat._m, self._m, (<Matrix>other)._m, self._rows, self._cols)
        elif isinstance(other, (int, float)):
            newMat = Matrix()
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
        """
        Subtract two matrices or subtract a scalar from all matrix elements.

        Args:
            other: Either a Matrix object or a scalar (int/float)

        Returns:
            A new Matrix containing the result

        Raises:
            ValueError: If matrix dimensions are incompatible or operation not supported

        Examples:
            >>> m1 = Matrix([[5, 6], [7, 8]])
            >>> m2 = Matrix([[1, 2], [3, 4]])
            >>> m3 = m1 - m2  # [[4, 4], [4, 4]]
            >>> m4 = m1 - 5   # [[0, 1], [2, 3]]
        """
        cdef unsigned int i
        cdef Matrix newMat
        if isinstance(other, Matrix):
            if (<Matrix>other)._rows != self._rows or (<Matrix>other)._cols != self._cols:
                raise ValueError("Dimmensions shoud be compatible")
            newMat = Matrix()
            newMat._m = <double *> PyMem_Malloc(self._rows * self._cols * sizeof(double))
            if not newMat._m:
                raise MemoryError()
            sub(newMat._m, self._m, (<Matrix>other)._m, self._rows, self._cols)
        elif isinstance(other, (int, float)):
            newMat = Matrix()
            newMat._m = <double *> PyMem_Malloc(self._rows * self._cols * sizeof(double))
            if not newMat._m:
                raise MemoryError()
            for i in range(self._cols * self._rows):
                newMat._m[i] = self._m[i] - other
        else:
            raise ValueError("Not implemented subtract")
        newMat._rows = self._rows
        newMat._cols = self._cols
        return newMat

    def __serialize__(self) -> list:
        """
        Serialize the matrix to a flat list.

        Returns:
            A 1D list containing all matrix elements in row-major order
        """
        cdef unsigned int i
        if self._rows == 0 or self._cols == 0 or self._m == NULL:
            return []
        return [self._m[i] for i in range(self._rows*self._cols)]

    def __eq__(Matrix self, Matrix other) -> bool:
        """
        Check if two matrices are equal (within numerical precision).

        Args:
            other: Another Matrix object

        Returns:
            True if matrices have same dimensions and all elements are equal
            within the precision threshold, False otherwise
        """
        cdef unsigned int i
        if self.size() != other.size():
            return False
        if self._rows == 0 or self._cols == 0:
            return True  # Both are empty with same dimensions
        if self._m == NULL or other._m == NULL:
            return self._m == other._m  # Both NULL or one is NULL
        return all([(fabs(self._m[i]-other._m[i]) < presition) for i in range(self._cols*self._rows)])

    def __ne__(Matrix self, Matrix other) -> bool:
        """
        Check if two matrices are not equal.

        Args:
            other: Another Matrix object

        Returns:
            True if matrices are different, False if they are equal
        """
        return not self.__eq__(other)

    def __neg__(Matrix self):
        """
        Negate all elements in the matrix (unary minus operator).

        Returns:
            A new Matrix with all elements negated

        Examples:
            >>> m = Matrix([[1, -2], [3, -4]])
            >>> m_neg = -m  # [[-1, 2], [-3, 4]]
        """
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

    def __mul__(Matrix self, other):
        """
        Multiply two matrices or multiply matrix by a scalar.

        Args:
            other: Either a Matrix object or a scalar (int/float)

        Returns:
            A new Matrix containing the result

        Raises:
            ArithmeticError: If matrix dimensions are incompatible for multiplication
            ValueError: If operation is not supported

        Examples:
            >>> m1 = Matrix([[1, 2], [3, 4]])
            >>> m2 = Matrix([[5, 6], [7, 8]])
            >>> m3 = m1 * m2  # Matrix multiplication
            >>> m4 = m1 * 2   # Scalar multiplication
        """
        cdef unsigned int row, col, j
        cdef double value
        cdef Matrix newMat

        if isinstance(other, Matrix):
            if self.cols != other.rows:
                raise ArithmeticError("The matrix dimmensions should be compatible")
            newMat = Matrix()
            newMat._m = <double *> PyMem_Malloc(other.cols * self.rows * sizeof(double))
            if not newMat._m:
                raise MemoryError()
            newMat._rows = self._rows
            newMat._cols = (<Matrix> other)._cols
            mult(newMat._m, self._m, self.rows, self.cols,
                 (<Matrix> other)._m, other.rows, other.cols)
        elif isinstance(other, (int, float)):
            newMat = Matrix()
            newMat._m = <double *> PyMem_Malloc(self.cols * self.rows * sizeof(double))
            if not newMat._m:
                raise MemoryError()
            newMat._rows = self._rows
            newMat._cols = self._cols
            multByScalar(newMat._m, self._m, self.rows, self.cols, other)
        else:
            raise ValueError("Not implemented multiply")
        return newMat

    def __rmul__(Matrix self, other):
        """
        Right multiplication (scalar * matrix or matrix * matrix).

        This method is called when the matrix is on the right side of the
        multiplication operator.

        Args:
            other: Either a Matrix object or a scalar (int/float)

        Returns:
            A new Matrix containing the result

        Examples:
            >>> m = Matrix([[1, 2], [3, 4]])
            >>> m2 = 2 * m  # Calls __rmul__
        """
        cdef unsigned int row, col, j
        cdef double value
        cdef Matrix newMat
        if isinstance(other, Matrix):
            if self.cols != other.rows:
                raise ArithmeticError("The matrix dimmensions should be compatible")
            newMat = Matrix()
            newMat._m = <double *> PyMem_Malloc(other.cols * self.rows * sizeof(double))
            if not newMat._m:
                raise MemoryError()
            newMat._rows = self._rows
            newMat._cols = (<Matrix> other)._cols
            mult(newMat._m, self._m, self.rows, self.cols,
                 (<Matrix> other)._m, other.rows, other.cols)
        elif isinstance(other, (int, float)):
            newMat = Matrix()
            newMat._m = <double *> PyMem_Malloc(self.cols * self.rows * sizeof(double))
            if not newMat._m:
                raise MemoryError()
            newMat._rows = self._rows
            newMat._cols = self._cols
            multByScalar(newMat._m, self._m, self.rows, self.cols, other)
        else:
            raise ValueError("Not implemented multiply")
        # newMat.pushdata(self._rows, self._cols, mm)
        return newMat


    def __bool__(Matrix self) -> bool:
        """
        Check if the matrix contains any non-zero elements.

        Returns:
            True if any element is non-zero (within precision), False if all zeros

        Examples:
            >>> m1 = Matrix([[0, 0], [0, 0]])
            >>> bool(m1)  # False
            >>> m2 = Matrix([[0, 1], [0, 0]])
            >>> bool(m2)  # True
        """
        cdef unsigned int i
        return any([fabs(self._m[i]) > presition for i in range(self._rows*self._cols)])

    def __div__(Matrix self, other):
        """
        Divide matrix by another matrix or by a scalar.

        For matrix division, this computes self * other^(-1).
        For scalar division, divides all elements by the scalar.

        Args:
            other: Either a Matrix object or a scalar (int/float)

        Returns:
            A new Matrix containing the result

        Raises:
            ValueError: If dividing by zero

        Examples:
            >>> m = Matrix([[2, 4], [6, 8]])
            >>> m2 = m / 2  # [[1, 2], [3, 4]]
        """
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
            newMat._rows = self._rows
            newMat._cols = self._cols
            for i in range(self._rows*self._cols):
                newMat._m[i] = self._m[i]/other
            return newMat

    def __truediv__(Matrix self, other):
        """
        True division operator (Python 3 style division).

        This is an alias for __div__.

        Args:
            other: Either a Matrix object or a scalar (int/float)

        Returns:
            A new Matrix containing the result
        """
        return self.__div__(other)

    def __copy__(Matrix self) -> Matrix:
        """
        Support for Python's copy module.

        Returns:
            A deep copy of the matrix
        """
        return self.copy()

    def __hash__(Matrix self):
        """
        Compute hash value for the matrix.

        Returns:
            Hash value based on matrix contents

        Warning:
            Matrix objects are mutable but hashable. If you use a Matrix as a
            dictionary key or in a set, DO NOT modify it afterwards, as this will
            break the dictionary/set. This is a design limitation.

            If the matrix is modified after being added to a dict/set, the hash
            will change and the object will become unreachable in the collection.

        Note:
            Allows matrices to be used as dictionary keys or in sets, but only
            if they are treated as immutable after hashing.
        """
        return hash(tuple(self.__serialize__()))


    @property
    def rows(Matrix self):
        """
        Get or set the number of rows in the matrix.

        When setting, the matrix is resized. New rows are filled with zeros,
        and excess rows are truncated.

        Returns:
            Number of rows in the matrix

        Examples:
            >>> m = Matrix([[1, 2], [3, 4]])
            >>> m.rows  # 2
            >>> m.rows = 3  # Adds a row of zeros
        """
        return self._rows

    @rows.setter
    def rows(Matrix self, unsigned int rows):
        cdef double * mat
        cdef unsigned int i,j,k
        if rows == 0 or self._cols == 0:
            raise ValueError("Cannot resize to zero dimensions")
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
        if self._m != NULL:
            PyMem_Free(self._m)
            self._m = NULL
        self._m = mat
        self._rows = rows

    @property
    def cols(Matrix self):
        """
        Get or set the number of columns in the matrix.

        When setting, the matrix is resized. New columns are filled with zeros,
        and excess columns are truncated.

        Returns:
            Number of columns in the matrix

        Examples:
            >>> m = Matrix([[1, 2], [3, 4]])
            >>> m.cols  # 2
            >>> m.cols = 3  # Adds a column of zeros
        """
        return self._cols

    @cols.setter
    def cols(Matrix self, unsigned int cols):
        cdef double * mat
        cdef unsigned int i, j, k
        if cols == 0 or self._rows == 0:
            raise ValueError("Cannot resize to zero dimensions")
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
        if self._m != NULL:
            PyMem_Free(self._m)
            self._m = NULL
        self._m = mat
        self._cols = cols


    def isSimetric(Matrix self):
        """
        Check if the matrix is symmetric.

        A matrix is symmetric if it equals its transpose, i.e., A[i,j] = A[j,i]
        for all i, j.

        Returns:
            True if the matrix is symmetric, False otherwise

        Examples:
            >>> m1 = Matrix([[1, 2], [2, 3]])
            >>> m1.isSimetric()  # True
            >>> m2 = Matrix([[1, 2], [3, 4]])
            >>> m2.isSimetric()  # False

        Note:
            Non-square matrices always return False.
        """
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
        """
        Compute the transpose of the matrix.
        
        The transpose swaps rows and columns: result[i,j] = original[j,i]
        
        Returns:
            A new Matrix object containing the transpose
        
        Raises:
            MemoryError: If memory allocation fails
        
        Examples:
            >>> m = Matrix([[1, 2, 3], [4, 5, 6]])
            >>> mt = m.transpose()
            >>> # mt is [[1, 4], [2, 5], [3, 6]]
        """
        cdef Matrix newMat = Matrix()
        newMat._rows = self._cols
        newMat._cols = self._rows
        newMat._m = <double *> PyMem_Malloc(self._rows*self._cols*sizeof(double))
        if not newMat._m:
            raise MemoryError()
        try:
            transpose(newMat._m, self._m, self._rows, self._cols)
            return newMat
        except Exception as err:
            del(newMat)
            raise err


    def adjugate(Matrix self):
        """
        Compute the adjugate (adjoint) matrix.

        The adjugate is the transpose of the cofactor matrix. It's used in
        calculating the matrix inverse.

        Returns:
            A new Matrix object containing the adjugate

        Raises:
            ArithmeticError: If the matrix is not square
            MemoryError: If memory allocation fails

        Examples:
            >>> m = Matrix([[1, 2], [3, 4]])
            >>> adj = m.adjugate()
        """
        cdef Matrix newMat
        if self._rows != self._cols:
            raise ArithmeticError("The matrix should be squared for adjugate calculation")
        newMat = Matrix()
        newMat._rows = self._rows
        newMat._cols = self._cols
        newMat._m = <double *> PyMem_Malloc(self._rows*self._cols*sizeof(double))
        if not newMat._m:
            raise MemoryError()
        try:
            adjugate(newMat._m, self._m, self._rows, self._cols)
            return newMat
        except Exception as err:
            del(newMat)
            raise err

    def determinant(Matrix self) -> float:
        """
        Calculate the determinant of the matrix.

        The determinant is a scalar value that can be computed from a square matrix.
        It has important properties in linear algebra.

        Returns:
            The determinant value as a float

        Raises:
            ArithmeticError: If the matrix is not square

        Examples:
            >>> m = Matrix([[1, 2], [3, 4]])
            >>> m.determinant()  # -2.0

        Note:
            Only defined for square matrices.
        """
        if self._rows != self._cols:
            raise ArithmeticError("Determinant is only defined for square matrices")
        return determinant(self._m, self._rows)

    def inverse(Matrix self):
        """
        Compute the inverse of the matrix.

        The inverse A^(-1) satisfies: A * A^(-1) = A^(-1) * A = I (identity matrix)

        Returns:
            A new Matrix object containing the inverse

        Raises:
            ArithmeticError: If matrix is not square or determinant is zero (singular)
            MemoryError: If memory allocation fails

        Examples:
            >>> m = Matrix([[1, 2], [3, 4]])
            >>> m_inv = m.inverse()
            >>> # m * m_inv should equal identity matrix
        """
        cdef Matrix newMat
        if self._rows != self._cols:
            raise ArithmeticError("The matrix should be squared for adjugate calculation")
        newMat = Matrix()
        newMat._rows = self._rows
        newMat._cols = self._cols
        newMat._m = <double *> PyMem_Malloc(self._rows*self._cols*sizeof(double))
        if not newMat._m:
            raise MemoryError()
        try:
            inverse(newMat._m, self._m, self._rows, self._cols)
            return newMat
        except Exception as err:
            del(newMat)
            raise err

    @staticmethod
    cdef c_zeros(int row, int col):
        """
        Create a zero matrix (C-level static method).
        
        Args:
            row: Number of rows
            col: Number of columns
        
        Returns:
            A new Matrix filled with zeros
        
        Raises:
            ValueError: If row or col is not positive
            MemoryError: If memory allocation fails
        """
        cdef int i
        cdef Matrix newMat
        if row <= 0 or col <= 0:
            raise ValueError("Row and col should be different than zero")
        newMat = Matrix()
        newMat._m = <double *> PyMem_Malloc(row*col*sizeof(double))
        if not newMat._m:
            raise MemoryError()
        for i in range(col*row):
            newMat._m[i] = 0.0
        newMat._rows = (<unsigned int>row)
        newMat._cols = (<unsigned int>col)
        return newMat

    @staticmethod
    def zeros(int row, int col):
        """
        Create a zero matrix.

        Args:
            row: Number of rows
            col: Number of columns

        Returns:
            A new Matrix filled with zeros

        Raises:
            ValueError: If row or col is not positive
            MemoryError: If memory allocation fails

        Examples:
            >>> m = Matrix.zeros(3, 4)
            >>> # Creates a 3x4 matrix filled with zeros
        """
        return Matrix.c_zeros(row, col)

    @staticmethod
    cdef c_identity(int n):
        """
        Create an identity matrix (C-level static method).
        
        Args:
            n: Size of the square identity matrix
        
        Returns:
            A new n×n identity Matrix (1s on diagonal, 0s elsewhere)
        
        Raises:
            ValueError: If n is not positive
            MemoryError: If memory allocation fails
        """
        cdef int i,j,k
        cdef Matrix newMat
        if n <= 0:
            raise ValueError("Dimmension should be a positive integer")
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

    @staticmethod
    def identity(int n):
        """
        Create an identity matrix.

        Args:
            n: Size of the square identity matrix

        Returns:
            A new n×n identity Matrix (1s on diagonal, 0s elsewhere)

        Raises:
            ValueError: If n is not positive
            MemoryError: If memory allocation fails

        Examples:
            >>> m = Matrix.identity(3)
            >>> # Creates [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        """
        return Matrix.c_identity(n)

