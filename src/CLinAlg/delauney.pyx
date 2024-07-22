from .cimport AlgVector, AlgVector2D
from cpython.mem cimport PyMem_Malloc, PyMem_Free

cdef void normalize_3d_to_2d(double * to_vec2d, double * vect3darray, int num_point, double * zero_3d,
                             double * axe_3d_u, double * axe_3d_v):
    "Function that normalize an array of 3d points into 2d points"
    cdef int i
    cdef double * v_sub = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double min_u, max_u, min_v, max_v, dX
    try:
        for i in range(num_point):
            AlgVector.sub(v_sub, &vect3darray[i * 3], zero_3d)
            to_vec2d[i * 2] = AlgVector.dot(v_sub, axe_3d_u)
            to_vec2d[i * 2 + 1] = AlgVector.dot(v_sub, axe_3d_v)
        # proceso de normalizado
        min_u = to_vec2d[0]
        max_u = to_vec2d[0]
        min_v = to_vec2d[1]
        max_v = to_vec2d[1]
        for i in range(1, num_point):
            min_u = min(min_u, to_vec2d[2 * i])
            max_u = max(max_u, to_vec2d[2 * i])
            min_v = min(min_v, to_vec2d[2 * i + 1])
            max_v = max(max_v, to_vec2d[2 * i + 1])
        dX = max(max_u - min_u, max_v - min_v)
        for i in range(num_point):
            to_vec2d[i * 2] = (to_vec2d[i * 2] - min_u) / dX
            to_vec2d[i * 2 + 1] = (to_vec2d[i * 2 + 1] - min_v) / dX
    finally:
        PyMem_Free(v_sub)

cdef bint should_swap(double * vCom, double * vA, double * vB, double * vP):
    "Return True if should swap the triangles vCom-vA-vB with oposite point vP"
    cdef double * v13 = <double *> PyMem_Malloc(2 * sizeof(double))
    cdef double * v23 = <double *> PyMem_Malloc(2 * sizeof(double))
    cdef double * v1p = <double *> PyMem_Malloc(2 * sizeof(double))
    cdef double * v2p = <double *> PyMem_Malloc(2 * sizeof(double))
    cdef double cosa, cosb, sinab
    try:
        AlgVector2D.sub(v13, vA, vCom)
        AlgVector2D.sub(v23, vB, vCom)
        AlgVector2D.sub(v1p, vA, vP)
        AlgVector2D.sub(v2p, vB, vP)
        cosa = v13[0] * v23[0] + v13[1] * v23[1]
        cosb = v2p[0] * v1p[0] + v1p[1] * v2p[1]
        if cosa >= 0 and cosb >= 0:
            return False
        if cosa < 0 and cosb < 0:
            return True
        sinab = (v13[0] * v23[1] - v23[0] * v13[1]) * cosb + (v2p[0] * v1p[1] - v1p[0] * v2p[1]) * cosa
        if sinab < 0:
            return True
        return False
    finally:
        PyMem_Free(v13)
        PyMem_Free(v23)
        PyMem_Free(v1p)
        PyMem_Free(v2p)

cdef void split_triangle_into_3(double * plist2d, int num_point, int * tri_vertex, int * tri_adjoin, list vertex_in_tri,
                                int triIndex, int pointIndex, int lastTriIndex):
    """
    Function that divide the triangle given in 3 with the inner point given by index
    :param plist2d: array with num_point*2 * double 
    :param num_point: number of elements
    :param tri_vertex: array with 3 * num_tri * int indexes of the vertex of triangle
    :param tri_adjoin: array with 3 * num_tri * int indexes of adjoin to each triangle
    :param vertex_in_tri: list of int with the vertexes within the initial triangle
    :param triIndex: Index of triangle that should be divied
    :param pointIndex: Index of point that will divide the triangle
    :param lastTriIndex: Last triangle added to tri_vertex array (after this index, the array is not initializated
    """
    cdef int i
    cdef int v0, v1,v2 # gets the ver
    cdef double * vP
    cdef double * v1 = <double *> PyMem_Malloc(2*sizeof(double))
    cdef double * v2 = <double *> PyMem_Malloc(2*sizeof(double))

    try:
        v0 = tri_vertex[triIndex*3]
        v1 = tri_vertex[triIndex*3+1]
        v2 = tri_vertex[triIndex*3+2]
        new_tri = [[tri_vertex[triIndex*3], tri_vertex[triIndex*3+1], pointIndex],
                   [tri_vertex[triIndex*3+1], tri_vertex[triIndex*3+2], pointIndex],
                   [tri_vertex[triIndex*3+2], tri_vertex[triIndex*3], pointIndex]]
        new_tri_index = [triIndex, lastTriIndex+1, lastTriIndex+2] # he agregado nuevos puntos
        vertex_in_new_tri = [[],[],[]]
        ptocheck = [i for i,x in enumerate(vertex_in_tri) if x == triIndex and i != pointIndex]
        vP = &plist2d[pointIndex*2]
        for k in range(2):
            if not ptocheck:
                break

    finally:
        PyMem_Free(v1)
        PyMem_Free(v2)




cdef (int * , int *) delauney_2d_constricted(double * plist2d, int num_point, int * contour, int num_contour):
    """
    Function that construct a 2d delauney triangulation with constrictions
    :param plist2d: points 2d 
    :param num_point: num points in array
    :param contour: index of points that create a contours. The index -1 divide the contour. The contour should be closed (last index == first)
    :param num_contour: num of points in the contour array
    """
    cdef int i, j, k
    cdef double * plist2d_ext = <double *> PyMem_Malloc((num_point+3) * 2 * sizeof(double))
    cdef int * avoid_swap = <int *> PyMem_Malloc((num_point + 3) * (num_point + 3) * sizeof(int))
    cdef int * used_vertex = <int *> PyMem_Malloc(num_point * sizeof(int))
    cdef int num_tri = 2 * (num_point + 3) - 2 - 3  # esquema de triangulo envolvente
    cdef int * tri_vertex = <int *> PyMem_Malloc (num_tri*3*sizeof(int))
    cdef int * tri_adjoin = <int *> PyMem_Malloc (num_tri*3*sizeof(int))

    if not num_contour:
        raise RuntimeError("This functions needs contour points")
    try:
        for i in range(num_point):
            for j in range(2):
                plist2d_ext[2*i+j] = plist2d[2*i+j]
        plist2d_ext[2*num_point] = -100.0
        plist2d_ext[2*num_point+1] = -100.0
        plist2d_ext[2*num_point+2] = 100.0
        plist2d_ext[2*num_point+3] = -100.0
        plist2d_ext[num_point * 2 + 4] = 0.0
        plist2d_ext[num_point * 2 + 5] = 100.0

        for j in range(num_point):
            used_vertex[j] = 0
            for k in range(num_point):
                avoid_swap[j * (num_point + 3) + k] = 0
        for j in range(num_point, num_point + 3):
            for k in range(num_point, num_point + 3):
                avoid_swap[j * (num_point + 3) + k] = 1

        for i in range(num_contour - 1):
            if contour[i] == -1 or contour[i + 1] == -1:
                continue
            avoid_swap[contour[i] * (num_point + 3) + contour[i + 1]] = 1
            avoid_swap[contour[i + 1] * (num_point + 3) + contour[i]] = 1

        # genero el array del tri_vertex
        tri_vertex[0] = num_point
        tri_vertex[1] = num_point + 1
        tri_vertex[2] = num_point + 2
        tri_adjoin[0] = -1
        tri_adjoin[1] = -1
        tri_adjoin[2] = -1
        vertex_in_tri = [0 for _ in range(num_point)] # lista con todos los puntos en el triangulo que los contine
        last_tri_index = 0
        i = 0
        while i < num_contour:
            if contour[i] == -1:
                i += 1
                continue
            if not used_vertex[contour[i]]:
                new_tri_index = split_triangle_into_3(plist2d, num_point, tri_vertex, tri_adjoin, vertex_in_tri, vertex_in_tri[contour[i]],
                                                      contour[i], last_tri_index)



    finally:
        PyMem_Free(avoid_swap)
        PyMem_Free(used_vertex)
        PyMem_Free(tri_adjoin)
        PyMem_Free(plist2d_ext)
