from cpython.mem cimport PyMem_Malloc, PyMem_Free
from libc.math cimport fabs, sqrt
from .AlgTool cimport presition

cdef class linspace:
    def __init__(self, start: float, stop: float, num: int = 50, endpoint =True):
        self.start = <double> start
        self._cur = <double> start
        self.stop = <double> stop
        self.delta = <double> ((stop - start) / num if not endpoint else (stop - start) / (num - 1))
        self.endpoint = <bint> endpoint
        self.num = <int> num
        # print("llego", self.start, self.stop, self.num, self.endpoint)

    def __iter__(self):
        return self

    def __len__(self) -> int:
        return self.num

    def __next__(self) -> float:
        if self.start <= self.stop:
            if self.endpoint:
                if self._cur <= (self.stop + presition):
                    self._tval = self._cur
                    self._cur += self.delta
                    return self._tval
            else:
                if self._cur < (self.stop - presition):
                    self._tval = self._cur
                    self._cur += self.delta
                    return self._tval
        else:
            if self.endpoint:
                if self._cur >= (self.stop - presition):
                    self._tval = self._cur
                    self._cur += self.delta
                    return self._tval
            else:
                if self._cur > (self.stop + presition):
                    self._tval = self._cur
                    self._cur += self.delta
                    return self._tval
        self._cur = self.start
        raise StopIteration

cdef class arange:
    def __init__(self, start: float, stop: float, step: float = 1.0):
        if (start > stop and step > 0) or (start < stop and step < 0):
            raise ValueError("start, stop and step should have compatible dimmensions")
        self.start = <double> start
        self.stop = <double> stop
        self.step = <double> fabs(step)
        self._cur = self.start

    def __iter__(self)-> arange:
        self._cur = self.start
        return self

    def __len__(self) -> int:
        return int(fabs(self.stop - self.start) / self.step) + 1

    def __next__(self) -> float:
        if self.start <= self.stop:  # ascendente
            if self._cur <= (self.stop + presition):
                self._tval = self._cur
                self._cur = self._cur + self.step
                return self._tval
        else:  #descendente
            if self._cur >= (self.stop - presition):
                self._tval = self._cur
                self._cur = self._cur - self.step
                return self._tval
        raise StopIteration

cdef class SplinePoint:
    def __init__(self, xData, yData, method='newton'):
        cdef unsigned int i
        if method not in ('newton', 'neville', 'cubic', 'linear', 'lagrange'):  # 'rational has been deactivated
            raise ValueError("Method %s NOT included as interpolation system" % method)
        if len(xData) != len(yData):
            raise ValueError("xValues and yValues arrays should have the same dimmension")
        self._xVal = <double *> PyMem_Malloc(len(xData) * sizeof(double))
        self._yVal = <double *> PyMem_Malloc(len(yData) * sizeof(double))
        if self._xVal == NULL or self._yVal == NULL:
            raise MemoryError()
        self._nsize = <unsigned int> len(xData)
        for i in range(self._nsize):
            self._xVal[i] = xData[i]
            self._yVal[i] = yData[i]
        if method == 'newton':
            self._met = 0
            self.updateNewtonCoef()
        elif method == 'neville':
            self._met = 1
        elif method == 'cubic':
            self._met = 2
            self.updateCubicCoef()
        elif method == 'linear':
            self._met = 3
        elif method == 'lagrange':
            self._met = 4
        else:
            raise NotImplementedError(
                "You should implement one of availables methods (newton, neville, cubic, linear, lagrange)")

    def __dealloc__(SplinePoint self):
        if self._xVal != NULL:
            PyMem_Free(self._xVal)
            self._xVal = NULL
        if self._yVal != NULL:
            PyMem_Free(self._yVal)
            self._yVal = NULL
        if self._tcoef != NULL:
            PyMem_Free(self._tcoef)
            self._tcoef = NULL

            #.....................C METHOD.....................................................
    cdef void updateNewtonCoef(SplinePoint self):
        cdef int i, k
        if self._tcoef != NULL:
            PyMem_Free(self._tcoef)
            self._tcoef = NULL
        self._tcoef = <double *> PyMem_Malloc(self._nsize * sizeof(double))
        for i in range(<int> self._nsize):
            self._tcoef[i] = self._yVal[i]
        for k in range(1, <int> self._nsize):
            for i in range(k, <int> self._nsize):
                self._tcoef[i] = (self._tcoef[i] - self._tcoef[k - 1]) / (self._xVal[i] - self._xVal[k - 1])

    cdef void updateCubicCoef(SplinePoint self):
        cdef double * c = <double *> PyMem_Malloc((self._nsize - 1) * sizeof(double))
        cdef double * d = <double *> PyMem_Malloc(self._nsize * sizeof(double))
        cdef double * e = <double *> PyMem_Malloc((self._nsize - 1) * sizeof(double))
        cdef double lam
        cdef int i
        if self._tcoef != NULL:
            PyMem_Free(self._tcoef)
            self._tcoef = NULL
        self._tcoef = <double *> PyMem_Malloc(self._nsize * sizeof(double))  #k
        for i in range(<int> (self._nsize - 1)):
            c[i] = 0.0
            e[i] = 0.0
        for i in range(<int> self._nsize):
            d[i] = 1.0
            self._tcoef[i] = 0.0
        try:
            if c == NULL or d == NULL or e == NULL:
                raise MemoryError()
            for i in range(0, self._nsize - 2):
                c[i] = self._xVal[i] - self._xVal[i + 1]
            for i in range(1, self._nsize - 1):
                d[i] = 2.0 * (self._xVal[i - 1] - self._xVal[i + 1])
                e[i] = self._xVal[i] - self._xVal[i + 1]
                self._tcoef[i] = 6.0 * (self._yVal[i - 1] - self._yVal[i]) / (self._xVal[i - 1] - self._xVal[i]) \
                                 - 6.0 * (self._yVal[i] - self._yVal[i + 1]) / (self._xVal[i] - self._xVal[i + 1])
            # descomposicion LU
            for i in range(1, <int> self._nsize):
                lam = c[i - 1] / d[i - 1]
                d[i] = d[i] - lam * e[i - 1]
                c[i - 1] = lam
            # resolucion LU(c,d,e,k=b)
            for i in range(1, <int> self._nsize):
                self._tcoef[i] = self._tcoef[i] - c[i - 1] * self._tcoef[i - 1]
            self._tcoef[self._nsize - 1] = self._tcoef[self._nsize - 1] / d[self._nsize - 1]
            for i in range(self._nsize - 2, -1, -1):
                self._tcoef[i] = (self._tcoef[i] - e[i] * self._tcoef[i + 1]) / d[i]
        finally:
            if c != NULL:
                PyMem_Free(c)
                c = NULL
            if d != NULL:
                PyMem_Free(d)
                d = NULL
            if e != NULL:
                PyMem_Free(e)
                e = NULL

    cdef double evalNewtonFunction(SplinePoint self, double x):
        cdef double p = self._tcoef[self._nsize - 1]
        cdef int i
        for i in range(1, <int> self._nsize):
            p = self._tcoef[self._nsize - 1 - i] + (x - self._xVal[self._nsize - 1 - i]) * p
        return p

    cdef double evalDerivative1NewtonFunction(SplinePoint self, double x):
        cdef double p = self._tcoef[self._nsize - 1]
        cdef double dp = 0.0
        cdef int i
        for i in range(1, <int> self._nsize):
            dp = p + dp * (x - self._xVal[self._nsize - 1 - i])
            p = self._tcoef[self._nsize - 1 - i] + (x - self._xVal[self._nsize - 1 - i]) * p
        return dp

    cdef double evalNevilleFunction(SplinePoint self, double x):
        cdef double * y = <double *> PyMem_Malloc(self._nsize * sizeof(double))
        cdef int i, k
        try:
            if y == NULL:
                raise MemoryError()
            for i in range(<int> self._nsize):
                y[i] = self._yVal[i]
            for k in range(1, <int> self._nsize):
                for i in range((<int> self._nsize) - k):
                    y[i] = ((x - self._xVal[i + k]) * y[i] + (self._xVal[i] - x) * y[i + 1]) / \
                           (self._xVal[i] - self._xVal[i + k])
            return y[0]
        finally:
            if y != NULL:
                PyMem_Free(y)
                y = NULL

    cdef double evalDerivative1NevilleFunction(SplinePoint self, double x):
        cdef int i, k
        cdef double * p = <double *> PyMem_Malloc(self._nsize * sizeof(double))
        cdef double * pp = <double *> PyMem_Malloc(self._nsize * sizeof(double))
        try:
            if p == NULL or pp == NULL:
                raise MemoryError()
            for i in range(<int> self._nsize):
                p[i] = self._yVal[i]
                pp[i] = 0.0
            for k in range(1, <int> self._nsize):
                for i in range((<int> self._nsize) - k):
                    #i==i;  j==k+i; i,j-1 == i ;; i+1,j == i+1
                    pp[i] = ((self._xVal[k + i] - x) * pp[i] - p[i] + (x - self._xVal[i]) * pp[i + 1] + p[i + 1]) / (
                            self._xVal[i + k] - self._xVal[i])
                    p[i] = ((x - self._xVal[i]) * p[i + 1] - (x - self._xVal[i + k]) * p[i]) / (
                            self._xVal[i + k] - self._xVal[i])
            return pp[0]
        finally:
            if p != NULL:
                PyMem_Free(p)
                p = NULL
            if pp != NULL:
                PyMem_Free(pp)
                pp = NULL

    cdef double evalLinearFunction(SplinePoint self, double x):
        cdef int i
        for i in range(self._nsize - 1):
            if self._xVal[i] <= x <= self._xVal[i + 1]:
                return self._yVal[i] + (self._yVal[i + 1] - self._yVal[i]) * (x - self._xVal[i]) / (
                        self._xVal[i + 1] - self._xVal[i])

    cdef double evalLagrangeFunction(SplinePoint self, double x):
        cdef int j, m
        cdef double l
        cdef double v = 0.0
        for j in range(<int> self._nsize):
            l = 1.0
            for m in range(<int> self._nsize):
                if j == m:
                    continue
                l *= (x - self._xVal[m]) / (self._xVal[j] - self._xVal[m])
            v += self._yVal[j] * l
        return v

    cdef double evalDerivative1LagrangeFunction(SplinePoint self, double x):
        cdef int j, i, m
        cdef double v = 0.0
        cdef double p, l
        for j in range(<int> self._nsize):
            l = 0.0
            for i in range(<int> self._nsize):
                if i == j:
                    continue
                p = 1 / (self._xVal[j] - self._xVal[i])
                for m in range(<int> self._nsize):
                    if m == i or m == j:
                        continue
                    p *= (x - self._xVal[m]) / (self._xVal[j] - self._xVal[m])
                l += p
            v += l * self._yVal[j]
        return v

    cdef double evalDerivative2LagrangeFunction(SplinePoint self, double x):
        cdef int j, i, m, n
        cdef double v = 0.0
        cdef double p, l, q
        for j in range(<int> self._nsize):
            p = 0.0
            for i in range(<int> self._nsize):
                if i == j:
                    continue
                l = 0.0
                for m in range(<int> self._nsize):
                    if m == i or m == j:
                        continue
                    q = 1 / (self._xVal[j] - self._xVal[m])
                    for n in range(<int> self._nsize):
                        if n in (j, i, m):
                            continue
                        q *= (x - self._xVal[n]) / (self._xVal[j] - self._xVal[n])
                    l += q
                p += l * (1 / (self._xVal[j] - self._xVal[i]))
            v += self._yVal[j] * p
        return v

    cdef double evalCubicFunction(SplinePoint self, double x):
        cdef int i, iLeft, iRight
        cdef double h
        iLeft = 0
        iRight = self._nsize - 1
        while True:
            if (iRight - iLeft) <= 1:
                i = iLeft
                break
            i = int((iRight + iLeft) / 2)
            if x < self._xVal[i]:
                iRight = i
            else:
                iLeft = i
        h = self._xVal[i] - self._xVal[i + 1]
        return ((x - self._xVal[i + 1]) ** 3 / h - (x - self._xVal[i + 1]) * h) * self._tcoef[i] / 6.0 \
            - ((x - self._xVal[i]) ** 3 / h - (x - self._xVal[i]) * h) * self._tcoef[i + 1] / 6.0 \
            + (self._yVal[i] * (x - self._xVal[i + 1]) - self._yVal[i + 1] * (x - self._xVal[i])) / h

    cpdef double derivative1(SplinePoint self, double x):
        """Function that return first derivate for the given point"""
        cdef int i
        if x < (self._xVal[0] - presition) or x > (self._xVal[self._nsize - 1] + presition):
            raise ValueError("X should be between x range")
        if self._met == 3:  #linear
            if fabs(x - self._xVal[self._nsize - 1]) < presition:
                return (self._yVal[self._nsize - 1] - self._yVal[self._nsize - 2]) / (
                        self._xVal[self._nsize - 1] - self._xVal[self._nsize - 2])

            for i in range(self._nsize - 1):
                if (self._xVal[i] - presition) <= x < self._xVal[i + 1]:
                    return (self._yVal[i + 1] - self._yVal[i]) / (
                            self._xVal[i + 1] - self._xVal[i])
        elif self._met == 2:
            for i in range(self._nsize - 1):
                if (self._xVal[i] - presition) <= x <= (self._xVal[i + 1] + presition):
                    return ((self._tcoef[i] / 6) * (
                            3 * (x - self._xVal[i + 1]) ** 2 / (self._xVal[i] - self._xVal[i + 1]) - (
                            self._xVal[i] - self._xVal[i + 1]))
                            - (self._tcoef[i + 1] / 6) * (
                                    3 * (x - self._xVal[i]) ** 2 / (self._xVal[i] - self._xVal[i + 1]) - (
                                    self._xVal[i] - self._xVal[i + 1]))
                            + (self._yVal[i] - self._yVal[i + 1]) / (self._xVal[i] - self._xVal[i + 1]))
        elif self._met == 1:
            # neville
            return self.evalDerivative1NevilleFunction(x)
        elif self._met == 0:
            # newton
            return self.evalDerivative1NewtonFunction(x)
        elif self._met == 4:
            return self.evalDerivative1LagrangeFunction(x)
        else:
            raise NotImplementedError("Not implemented derivative for non cubic function")

    cpdef double derivative2(SplinePoint self, double x):
        """Functions that returns second derivative for the given point"""
        cdef int i
        if self._met == 3:
            return 0.0
        elif self._met == 2:
            if x < self._xVal[0] or x > self._xVal[-1]:
                raise ValueError("X should be between x range")
            for i in range(<int> self._nsize - 1):
                if self._xVal[i] <= x <= self._xVal[i + 1]:
                    return (self._tcoef[i] * (x - self._xVal[i + 1]) / (self._xVal[i] - self._xVal[i + 1]) -
                            self._tcoef[i + 1] * (x - self._xVal[i]) / (self._xVal[i] - self._xVal[i + 1]))
        elif self._met == 4:
            return self.evalDerivative2LagrangeFunction(x)
        else:
            raise NotImplementedError("Not implemented derivative for non cubic function")

    #...................PYTHON METHODS..................................................

    @property
    def method(self):
        if self._met == 0:
            return 'newton'
        elif self._met == 1:
            return 'neville'
        elif self._met == 2:
            return 'cubic'
        elif self._met == 3:
            return 'linear'
        elif self._met == 4:
            return 'lagrange'

    def __call__(self, x: float) -> float:
        if self._met == 0:
            return self.evalNewtonFunction(x)
        elif self._met == 1:
            return self.evalNevilleFunction(x)
        elif self._met == 2:
            return self.evalCubicFunction(x)
        elif self._met == 3:
            return self.evalLinearFunction(x)
        elif self._met == 4:
            return self.evalLagrangeFunction(x)

    def root(self, method='incremental', **kwargs):
        "Function that return roots for selected curve"
        # methods: 'incremental': args: interval (a,b), step, presition
        # methods: 'mpmath' call to mpmath solver
        # f_eval is a 'y' value thath compute roots for y-f_eval = f(x)-f_eval
        # methods: 'newton
        # kwargs['roots'] is a lista that containts current roots
        if method == 'incremental':
            a = kwargs.get('a') if kwargs.get('a') is not None else self._xVal[0]
            b = kwargs.get('b') if kwargs.get('b') is not None else self._xVal[self._nsize - 1]
            dx = kwargs.get('step') or (self._xVal[self._nsize - 1] - self._xVal[0]) / 10
            presition = kwargs.get('presition') or 0.001
            roots = kwargs.get('roots') if kwargs.get('roots') is not None else []
            f_eval = kwargs.get('f_eval', 0.0)
            x1 = a
            f1 = self(x1) - f_eval
            x2 = a + dx
            f2 = self(x2) - f_eval
            while True:
                if x1 >= b:
                    return roots
                if sign(f1) != sign(f2):
                    if fabs(f1 - f2) < presition:
                        roots.append(x1 - (x2 - x1) * f1 / (f2 - f1))
                    else:
                        self.root(method='incremental', a=x1, b=x2, step=(x2 - x1) / 10, presition=presition,
                                  roots=roots, f_eval=f_eval)
                x1 = x2
                f1 = f2
                x2 = x1 + dx;
                f2 = self(x2) - f_eval
        elif method == 'mpmath':
            import mpmath
            a = kwargs.get('a') if kwargs.get('a') is not None else self._xVal[0]
            b = kwargs.get('b') if kwargs.get('b') is not None else self._xVal[self._nsize - 1]
            f_eval = kwargs.get('f_eval', 0.0)
            if f_eval:
                newfunct = SplinePoint(self.xData, [y - f_eval for y in self.yData], method=self.method)
            else:
                newfunct = self
            roots = kwargs.get('roots') if kwargs.get('roots') is not None else []
            step = (b - a) / 1000
            root = float(mpmath.findroot(self, a, solver='mnewton'))
            if a <= root <= b:
                roots.append(root)
                if root > a + step:
                    newfunct.root(method='mpmath', a=a, b=root - step, roots=roots)
                if root < b - step:
                    newfunct.root(method='mpmath', a=root + step, b=b, roots=roots)
            return roots

    def solve(self, y):
        "Return x value that y = self(x)"
        return self.root(f_eval=y)

    def copy(SplinePoint self):
        cdef unsigned int i
        return SplinePoint([self._xVal[i] for i in range(self._nsize)], [self._yVal[i] for i in range(self._nsize)],
                           self.method)

    @property
    def xData(SplinePoint self):
        cdef unsigned int i
        data = tuple(self._xVal[i] for i in range(self._nsize))
        return data

    @property
    def yData(SplinePoint self):
        cdef unsigned int i
        data = tuple(self._yVal[i] for i in range(self._nsize))
        return data

    def __repr__(SplinePoint self):
        return 'SplinePoint({})'.format(self.method)

    def __str__(SplinePoint self):
        return 'SplinePoint({})'.format(self.method)

    def __copy__(SplinePoint self):
        return self.copy()

    def __deepcopy__(self, memodict={}):
        return self.copy()

    def __json__(self):
        cdef unsigned int i
        return {'__jsoncls__': 'CLinAlg.ToolMath:SplinePoint.from_JSON', 'xVal': self.xData, 'yVal': self.yData,
                'method': self.method}

    @classmethod
    def from_JSON(cls, jsondict):
        cdef SplinePoint spline_func = SplinePoint(jsondict['xVal'], jsondict['yVal'], jsondict['method'])
        return spline_func

cpdef int sign(double x):
    return 1 if x >= 0 else -1

cdef class DoubleFuncWrapper:
    def __cinit__(self, pyfunction=None):
        self.pyfunct = pyfunction

    cdef double getValue(DoubleFuncWrapper self, double arg):
        if self.cyfunct:
            return self.cyfunct(arg)
        else:
            return self.pyfunct(arg)

    def __unsafe_set(self, ptr):
        "W.__unsafe_set(ctypes.addressof(cythonfunction)) can set pointer to c function within python"
        self.cyfunct = <func_t> <void *> <size_t> ptr

cdef DoubleFuncWrapper make_cyWrapper(func_t f):
    cdef DoubleFuncWrapper W = DoubleFuncWrapper()
    W.cyfunct = f
    return W

cdef DoubleFuncWrapper make_pyWrapper(object f):
    cdef DoubleFuncWrapper W = DoubleFuncWrapper(f)
    return W

cdef double cySolverNewton(DoubleFuncWrapper func, double Xmax, double Xmin, double tolerance) except? -1:
    #func: funcion que debe devolver 0 en el punto buscado
    #resuelve por el metodo newton
    cdef unsigned int count = 0
    cdef double V, Vmin, Vmax
    Vmin = func.getValue(Xmin)
    if fabs(Vmin) < tolerance:
        return Xmin
    Vmax = func.getValue(Xmax)
    if fabs(Vmax) < tolerance:
        return Xmax
    if not (min(Vmax, Vmin) < 0 < max(Vmax, Vmin)):  #no va a encontrar el 0 en el rango
        raise RuntimeError("Not possible to find 0 with this function")
    X = Xmin - Vmin * (Xmax - Xmin) / (Vmax - Vmin)
    V = func.getValue(X)
    while fabs(V) > tolerance and count < 200:  # and fabs(Xmax-Xmin) > tolerance2 :
        #Xo = X; Vo = V
        if sign(V) == sign(Vmin):
            Vmin = V
            Xmin = X
        elif sign(V) == sign(Vmax):
            Vmax = V
            Xmax = X
        X = Xmin - Vmin * (Xmax - Xmin) / (Vmax - Vmin)
        V = func.getValue(X)
        count += 1
    if count == 200:
        raise RuntimeError("Not solution finded")
    return X

def solverNewton(pyfunc, Xmax: float, Xmin: float, tolerance: float = 1e-6):
    cdef DoubleFuncWrapper wrapper = make_pyWrapper(pyfunc)
    return cySolverNewton(wrapper, Xmax, Xmin, tolerance)

cdef void _getIntersection2D(double * topoint, double * v1, double * dirV1, double * v2, double * dirV2) except*:
    "Unsafe because dont check if lines are paralel"
    cdef double mnu, lamb
    mnu = dirV2[1] - dirV1[1] * dirV2[0] / dirV1[0]
    mnu = (v1[1] + dirV1[1] * v2[0] / dirV1[0] - dirV1[1] * v1[0] / dirV1[0] - v2[1]) / mnu
    lamb = (v2[0] + dirV2[0] * mnu - v1[0]) / dirV1[0]
    topoint[0] = v1[0] + dirV1[0] * lamb
    topoint[1] = v1[1] + dirV1[1] * lamb

cdef void _normalize2D(double * v) except*:
    cdef double modulo
    modulo = sqrt(v[0] * v[0] + v[1] * v[1])
    v[0] = v[0] / modulo
    v[1] = v[1] / modulo

cdef double cySolverMinimumQuasiParabolic(DoubleFuncWrapper func, double Xmax, double Xmin, double tol) except? -1:
    cdef unsigned int count
    cdef double step
    cdef unsigned int i, j
    cdef double * v1 = <double *> PyMem_Malloc(2 * sizeof(double))
    cdef double * v2 = <double *> PyMem_Malloc(2 * sizeof(double))
    cdef double * dirV1 = <double *> PyMem_Malloc(2 * sizeof(double))
    cdef double * dirV2 = <double *> PyMem_Malloc(2 * sizeof(double))
    cdef double * vX = <double *> PyMem_Malloc(2 * sizeof(double))
    cdef double * dirVX = <double *> PyMem_Malloc(2 * sizeof(double))
    cdef double * p1 = <double *> PyMem_Malloc(2 * sizeof(double))
    cdef double * p2 = <double *> PyMem_Malloc(2 * sizeof(double))
    cdef bint spin  # 1 if +step, 0 if -step
    step = fabs(Xmax - Xmin) / 100
    count = 0
    try:
        if v1 == NULL or v2 == NULL or dirV1 == NULL or dirV2 == NULL or vX == NULL or dirVX == NULL or p1 == NULL or p2 == NULL:
            raise MemoryError()
        v1[0] = Xmax
        v1[1] = func.getValue(Xmax)
        dirV1[0] = -step
        dirV1[1] = func.getValue(Xmax - step) - v1[1]
        _normalize2D(dirV1)
        v2[0] = Xmin
        v2[1] = func.getValue(Xmin)
        dirV2[0] = step
        dirV2[1] = func.getValue(Xmin + step) - v2[1]
        _normalize2D(dirV2)
        if fabs(dirV1[1]) < tol or fabs(dirV2[1]) < tol:  # las tangentes se encuentran planas
            return Xmax if v1[1] < v2[1] else Xmin
        try:
            _getIntersection2D(vX, v1, dirV1, v2, dirV2)
        except:
            return Xmax if v1[1] < v2[1] else Xmin
        if not (Xmin <= vX[0] <= Xmax):
            return Xmax if v1[1] < v2[1] else Xmin
        vX[1] = func.getValue(vX[0])
        if vX[1] > min(v1[1], v2[1]):
            return Xmax if v1[1] < v2[1] else Xmin
        dirVX[0] = step
        dirVX[1] = func.getValue(vX[0] + step * 0.5) - func.getValue(vX[0] - step * 0.5)
        # print("vXinicial", vX[0], vX[1], dirVX[0], dirVX[1])
        _normalize2D(dirVX)
        while fabs(dirVX[1]) > tol and count < 200 and (fabs(v1[0] - vX[0]) > tol and fabs(v2[0] - vX[0]) > tol):
            try:
                _getIntersection2D(p1, v1, dirV1, vX, dirVX)
            except:
                return v1[0] if v1[1] < v2[1] else v2[0]
            try:
                _getIntersection2D(p2, v2, dirV2, vX, dirVX)
            except:
                return v1[0] if v1[1] < v2[1] else v2[0]
            if p1[1] < p2[1]:
                for i in range(2):
                    v2[i] = vX[i]
                    dirV2[i] = dirVX[i]
                vX[0] = p1[0]
                vX[1] = func.getValue(p1[0])
                dirVX[0] = -step
                dirVX[1] = func.getValue(p1[0] - step * 0.5) - func.getValue(p1[0] + step * 0.5)
            else:
                for i in range(2):
                    v1[i] = vX[i]
                    dirV1[i] = dirVX[i]
                vX[0] = p2[0]
                vX[1] = func.getValue(p2[0])
                dirVX[0] = step
                dirVX[1] = func.getValue(p2[0] + step * 0.5) - func.getValue(p2[0] - step * 0.5)
            _normalize2D(dirVX)
            print("vX", vX[0], vX[1], dirVX[0], dirVX[1])
            count += 1
        if count == 200:
            raise RuntimeError("Not founded minimun of quasiparabolic after 200 iterations")
        return vX[0]
    finally:
        if v1 != NULL:
            PyMem_Free(v1)
            v1 = NULL
        if v2 != NULL:
            PyMem_Free(v2)
            v2 = NULL
        if dirV1 != NULL:
            PyMem_Free(dirV1)
            dirV1 = NULL
        if dirV2 != NULL:
            PyMem_Free(dirV2)
            dirV2 = NULL
        if vX != NULL:
            PyMem_Free(vX)
            vX = NULL
        if dirVX != NULL:
            PyMem_Free(dirVX)
            dirVX = NULL
        if p1 != NULL:
            PyMem_Free(p1)
            p1 = NULL
        if p2 != NULL:
            PyMem_Free(p2)
            p2 = NULL

def solverMinimumQuasiParabolic(pyfunc: object, Xmax: float, Xmin: float, tol: float = 1e-6):
    cdef DoubleFuncWrapper wrapper = make_pyWrapper(pyfunc)
    return cySolverMinimumQuasiParabolic(wrapper, Xmax, Xmin, tol)

cdef tuple cyFindLimits(DoubleFuncWrapper func, double x_max, double x_min, int shift, double tolerance):
    "funcion que encuentra los limites de validez de una ecuacion"
    cdef list pos = []
    cdef double delta = (x_max - x_min)
    cdef object val
    while fabs(delta) > tolerance:
        delta = 0.5 * delta
        pos = [x + delta for x in pos]
        pos.insert(0, x_min + delta)
        if shift == 1:
            for x in pos:
                try:
                    val = func.getValue(x)
                except:
                    val = None
                if not val is None:
                    xfin, fin = cyFindLimits(func, x, x - delta, 1, tolerance)
                    return (xfin, fin) if not fin is None else (x, val)
        elif shift == -1:
            for x in reversed(pos):
                try:
                    val = func.getValue(x)
                except:
                    val = None
                if not val is None:
                    xfin, fin = cyFindLimits(func, x + delta, x, -1, tolerance)
                    return (xfin, fin) if not fin is None else (x, val)
    return 0.5 * (x_max + x_min), func(0.5 * (x_max + x_min))

def findLimits(pyfunc: object, Xmax, Xmin, shift, tolerance=1e-6):
    cdef DoubleFuncWrapper wrapper = make_pyWrapper(pyfunc)
    return cyFindLimits(wrapper, Xmax, Xmin, shift, tolerance)
