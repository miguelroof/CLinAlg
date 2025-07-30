cimport cython

cdef class linspace:
    cdef double start
    cdef double _cur
    cdef double _tval
    cdef double stop
    cdef double delta
    cdef bint endpoint
    cdef int num

cdef class arange:
    cdef double start
    cdef double stop
    cdef double _cur
    cdef double _tval
    cdef double step

cdef class SplinePoint:
    cdef double * _xVal
    cdef double * _yVal
    cdef double * _tcoef
    cdef unsigned int _met
    cdef unsigned int _nsize
    cdef void updateNewtonCoef(SplinePoint self)
    cdef void updateCubicCoef(SplinePoint self)
    cdef double evalNewtonFunction(SplinePoint self, double x)
    cdef double evalDerivative1NewtonFunction(SplinePoint self, double x)
    cdef double evalNevilleFunction(SplinePoint self, double x)
    cdef double evalDerivative1NevilleFunction(SplinePoint self, double x)
    cdef double evalLinearFunction(SplinePoint self, double x)
    cdef double evalCubicFunction(SplinePoint self, double x)
    cdef double evalLagrangeFunction(SplinePoint self, double x)
    cdef double evalDerivative1LagrangeFunction(SplinePoint self, double x)
    cdef double evalDerivative2LagrangeFunction(SplinePoint self, double x)
    cpdef double derivative1(SplinePoint self, double x)
    cpdef double derivative2(SplinePoint self, double x)

cpdef int sign(double x)

ctypedef double (* func_t)(double)

cdef class DoubleFuncWrapper:
    cdef func_t cyfunct
    cdef object pyfunct
    cdef double getValue(DoubleFuncWrapper self, double arg)

cdef DoubleFuncWrapper make_cyWrapper(func_t f)

cdef DoubleFuncWrapper make_pyWrapper(object f)

cdef double cySolverNewton(DoubleFuncWrapper func, double Xmax, double Xmin, double tolerance) except? -1

cdef double cySolverMinimumQuasiParabolic(DoubleFuncWrapper func, double Xmax, double Xmin, double tol) except? -1

cdef tuple cyFindLimits(DoubleFuncWrapper func, double x_max, double x_min, int shift, double tolerance)