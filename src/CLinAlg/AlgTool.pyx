from cpython.mem cimport PyMem_Malloc, PyMem_Free
cdef double presition = 1e-6

cpdef getPresition():
    return presition

cpdef setPresition(double presval):
    global presition
    presition = presval

