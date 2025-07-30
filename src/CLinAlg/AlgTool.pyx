cdef double M_PI = 3.141592653589793238

cdef double presition = 1e-6

cpdef float getPresition():
    global presition
    return <float>presition

cpdef void setPresition(float presval):
    global presition
    presition = <double>presval

