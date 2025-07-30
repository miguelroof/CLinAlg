from Cython.Compiler.ExprNodes import NullNode
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from . cimport AlgVector, AlgVector2D
from .AlgMatrix cimport Matrix
from .ToolMath cimport sign
from .AlgTool cimport presition

cdef double *  normalize_3d_to_2d(double * to_vec2d, double * vect3darray, int num_point,
                             double * zero_3d, double * axe_3d_u, double * axe_3d_v):
    "Function that normalize an array of 3d points into 2d points, using plane defined by zero, axe U and axe V"
    cdef int i
    cdef double * v_sub
    cdef double min_u, max_u, min_v, max_v, dX
    cdef double to_vec2d
    try:
        to_vec2d = <double *> PyMem_Malloc(2 * num_point * sizeof(double))
        if to_vec2d == NULL:
            raise MemoryError()
        v_sub = <double *> PyMem_Malloc(3 * sizeof(double))
        if v_sub == NULL:
            raise MemoryError()
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
        return to_vec2d
    except Exception as err:
        if to_vec2d != NULL:
            PyMem_Free(to_vec2d)
            to_vec2d = NULL
        raise err
    finally:
        if v_sub != NULL:
            PyMem_Free(v_sub)
            v_sub = NULL


cdef bint is_inside(double * vA, double * vB, double * vC, double * vP):
    "Check if point p is inside triangle vA, vB vC in 2d"
    cdef double v1 = <double *> PyMem_Malloc(2* sizeof(double))
    cdef double v2 = <double *> PyMem_Malloc(2*sizeof(double))
    cdef double v3 = <double *> PyMem_Malloc(2*sizeof(double))
    cdef double v12, v23, v31
    try:
        if v1 == NULL or v2 == NULL or v3 == NULL:
            raise MemoryError()
        AlgVector2D.sub(v1, vA, vP)
        AlgVector2D.sub(v2, vB, vP)
        AlgVector2D.sub(v3, vC, vP)
        v12 = AlgVector2D.cross(v1,v2)
        v23 = AlgVector2D.cross(v2,v3)
        v31 = AlgVector2D.cross(v3,v1)
        return sign(v12) == sign(v23) and sign(v12) == sign(v31)
    finally:
        if v1 != NULL:
            PyMem_Free(v1)
            v1 = NULL
        if v2 != NULL:
            PyMem_Free(v2)
            v2 = NULL
        if v3 != NULL:
            PyMem_Free(v3)
            v3 = NULL



cdef bint should_swap(double * vCom, double * vA, double * vB, double * vP):
    "Return True if should swap the triangles vCom-vA-vB with oposite point vP"
    cdef double * v13 = <double *> PyMem_Malloc(2 * sizeof(double))
    cdef double * v23 = <double *> PyMem_Malloc(2 * sizeof(double))
    cdef double * v1p = <double *> PyMem_Malloc(2 * sizeof(double))
    cdef double * v2p = <double *> PyMem_Malloc(2 * sizeof(double))
    cdef double cosa, cosb, sinab
    try:
        if v13 == NULL or v23 == NULL or v1p == NULL or v2p == NULL:
            raise MemoryError()
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
        if v13 != NULL:
            PyMem_Free(v13)
            v13 = NULL
        if v23 != NULL:
            PyMem_Free(v23)
            v23 = NULL
        if v1p != NULL:
            PyMem_Free(v1p)
            v1p = NULL
        if v2p != NULL:
            PyMem_Free(v2p)
            v2p = NULL

cdef void split_triangle_into_3(double * plist2d, int num_point,
                                int * tri_vertex, int num_tri, int last_tri,
                                int * tri_adjoin, int num_adj, int last_adj,
                                int * vertex_in_tri,
                                int triIndex,
                                int pointIndex):
    """
    Function that divide the triangle given in 3 with the inner point given by index
    :param plist2d: array with num_point*2 * double 
    :param num_point: number of elements
    :param tri_vertex: array with 3 * num_tri * int indexes of the vertex of triangle
    :param tri_adjoin: array with 3 * num_tri * int indexes of adjoin to each triangle
    :param vertex_in_tri: matrix with num_tri rows * num_point columns
    :param triIndex: Index of triangle that should be divied
    :param pointIndex: Index of point that will divide the triangle
    :param lastTriIndex: Last triangle added to tri_vertex array (after this index, the array is not initializated
    """
    cdef int i
    cdef int v0, v1,v2 # gets the vertices of triangle
    cdef double * vP
    cdef double * v1 = <double *> PyMem_Malloc(2*sizeof(double))
    cdef double * v2 = <double *> PyMem_Malloc(2*sizeof(double))
    try:
        if v1 == NULL or v2 == NULL:
            raise MemoryError()
        if not is_inside(&plist2d[2*tri_vertex[triIndex*3]],
                         &plist2d[2*tri_vertex[triIndex*3+1]],
                         &plist2d[2*tri_vertex[triIndex*3+2]],
                         &plist2d[2*pointIndex]):
            raise RuntimeError("Point to split should be within the triangle")

        newTri = [[tri_vertex[triIndex*3], tri_vertex[triIndex*3 + 1], pointIndex],
                  [tri_vertex[triIndex*3+1], tri_vertex[triIndex*3 + 2], pointIndex],
                  [tri_vertex[triIndex*3+2], tri_vertex[triIndex*3], pointIndex]]
        # reemplazo las posiciones en la matriz de triangulos

        newTriIndex = [triIndex, lastTriIndex + 1, lastTriIndex + 2]
        vertex_in_new_Tri = [[], [], []]
        ptocheck = [i for i, x in enumerate(vertex_in_tri) if x == triIndex and i != pointIndex]
        vP = plist2d[pointIndex]
        for k in range(2):
            if not ptocheck:
                break
            v1 = sub2D(plist2d[newTri[k][0]], vP)
            v2 = sub2D(plist2d[newTri[k][1]], vP)
            j = 0
            while j < len(ptocheck):
                v3 = sub2D(plist2d[ptocheck[j]], plist2d[newTri[k][2]])
                if cross2D(v1, v3) >= 0 and cross2D(v3, v2) > 0:
                    vertex_in_new_Tri[k].append(ptocheck.pop(j))
                else:
                    j += 1
        vertex_in_new_Tri[2] = ptocheck  # por descarte, todos los demas
        # ahora chequeo los trianculos adjacentes
        adjoin_to_newTri = [[tri_adjoin[triIndex][0], lastTriIndex + 1, lastTriIndex + 2],
                            [tri_adjoin[triIndex][1], lastTriIndex + 2, triIndex],
                            [tri_adjoin[triIndex][2], triIndex, lastTriIndex + 1]]
        for k in range(3):
            tri_vertex[newTriIndex[k]] = newTri[k]
            tri_adjoin[newTriIndex[k]] = adjoin_to_newTri[k]
            if adjoin_to_newTri[k][0] != -1:  # el antiguo colega de barrio
                old_adjoin = tri_adjoin[adjoin_to_newTri[k][0]]
                old_adjoin[old_adjoin.index(triIndex)] = newTriIndex[k]  # reemplazo el nombre
            for i in vertex_in_new_Tri[k]:
                vertex_in_tri[i] = newTriIndex[k]
        vertex_in_tri[pointIndex] = -1  # lo saco del listado sin usar
        return newTriIndex

    finally:
        if v1 != NULL:
            PyMem_Free(v1)
            v1 = NULL
        if v2 != NULL:
            PyMem_Free(v2)
            v2 = NULL

cdef void check_and_swap(plist2d, tri_vertex, tri_adjoin, vertex_in_tri, triP, avoid_swap=None):
    pass

cdef void swap_triangles(plist2d, tri_vertex, tri_adjoin, triK, triL, vertex_in_tri=None):
    pass

cdef void inspect_grouped_triangles(tri_vertex, tri_adjoin, avoid_swap, triangle_group, start_tri):
    pass



cdef void force_edge(double * plist2d, int num_point,
                     int * tri_vertex, int num_tri, int last_tri,
                     int * tri_adjoin, int last_adj,
                     int iA, int iB, vertex_in_tri=None):
    "Funcion que fuerza la conexion entre dos indices dados"
    vLine = normalize2D(sub2D(plist2d[iB], plist2d[iA]))
    crossEdges = []
    for i in range(len(tri_vertex)):
        if not iA in tri_vertex[i]:
            continue
        curI = i
        clockWise = True
        while True:
            if iB in tri_vertex[curI]:
                return True
            pA = tri_vertex[curI].index(iA)
            v1 = sub2D(plist2d[tri_vertex[curI][(pA + 1) % 3]], plist2d[iA])
            v2 = sub2D(plist2d[tri_vertex[curI][(pA + 2) % 3]], plist2d[iA])
            crossV1 = cross2D(v1, vLine)
            crossV2 = cross2D(vLine, v2)
            if abs(crossV1) < tolerance and dot2D(v1, vLine) > tolerance:
                return forceEdge(plist2d, tri_vertex, tri_adjoin, tri_vertex[curI][(pA + 1) % 3], iB)
            if abs(crossV2) < tolerance and dot2D(vLine, v2) > tolerance:
                return forceEdge(plist2d, tri_vertex, tri_adjoin, tri_vertex[curI][(pA + 2) % 3], iB)
            if crossV1 > tolerance and crossV2 > tolerance:
                triA = curI
                break
            if clockWise:
                curI = tri_adjoin[curI][pA]
                if curI == -1:  # El lado adjacente es de borde. tengo que invertir el sentido
                    pA = tri_vertex[i].index(iA)
                    curI = tri_adjoin[i][(pA + 2) % 3]
                    clockWise = False
            else:
                curI = tri_adjoin[curI][(pA + 2) % 3]
                if curI == -1:
                    raise RuntimeError("not founded solution for force edge algorithm")
        break
    crossEdges.append((triA, (pA + 1) % 3))
    while True:
        iK = tri_vertex[crossEdges[-1][0]][crossEdges[-1][1]]
        triK = tri_adjoin[crossEdges[-1][0]][crossEdges[-1][1]]
        pK = (tri_vertex[triK].index(iK) + 1) % 3
        if tri_vertex[triK][pK] == iB:
            break
        v1 = sub2D(point2d[tri_vertex[triK][pK]], plist2d[iA])
        crossV1 = cross2D(vLine, v1)
        if crossV1 > tolerance:
            crossEdges.append((triK, (pK + 2) % 3))
        elif crossV1 < -tolerance:
            crossEdges.append((triK, pK))
        else:  # el punto pK se encuentra en la linea
            forceEdge(plist2d, tri_vertex, tri_adjoin, iA, iK)
            forceEdge(plist2d, tri_vertex, tri_adjoin, iK, iB)
            return True
    i = 0
    newEdgesToAnalize = []
    while crossEdges:
        triK, pA = crossEdges.pop(i)
        triL = tri_adjoin[triK][pA]
        pP = (pA - 1) % 3
        pOP = (tri_vertex[triL].index(tri_vertex[triK][pA]) + 1) % 3
        v1 = sub2D(plist2d[tri_vertex[triK][pA]], plist2d[tri_vertex[triK][pP]])
        v2 = sub2D(plist2d[tri_vertex[triK][(pA + 1) % 3]], plist2d[tri_vertex[triK][pP]])
        vL = sub2D(plist2d[tri_vertex[triL][pOP]], plist2d[tri_vertex[triK][pP]])
        if cross2D(v1, vL) < tolerance:  # no convexo, vuelvo a poner
            crossEdges.insert(i, (triK, pA))
            i += 1
        elif cross2D(vL, v2) < tolerance:  # no convexo, vuelvo a poner el contador a 0
            crossEdges.insert(i, (triK, pA))
            i += 1
        else:
            # PB=1, PA=2, BOP=3, AOP = 6
            fval = 0
            if tri_vertex[triK][pP] == iA:
                pass
            else:
                if tri_vertex[crossEdges[i - 1][0]][crossEdges[i - 1][1]] == tri_vertex[triK][pP]:
                    fval += 1  # viene por PB
                else:
                    fval += 2  # viene por PA
            if tri_vertex[triL][pOP] == iB:
                pass
            else:
                if crossEdges[i][1] == pOP:  # sale por B-OP
                    fval += 3
                else:
                    fval += 6
            if fval == 0:
                pass
            elif fval == 1:
                newEdgesToAnalize.append((triL, 1))
            elif fval == 2:
                newEdgesToAnalize.append((triK, 2))
            elif fval == 3:
                crossEdges.pop(i)
                crossEdges.insert(0, (triK, 0))
                newEdgesToAnalize.append((triL, 1))
            elif fval == 6:
                crossEdges.pop(i)
                crossEdges.insert(0, (triL, 0))
                newEdgesToAnalize.append((triK, 2))
            elif fval == 4:
                crossEdges.pop(i)
                crossEdges.insert(i, (triK, 0))
                newEdgesToAnalize.append((triL, 1))
            elif fval == 5:
                crossEdges.pop(i)
                crossEdges.insert(i, (triL, 1))
                crossEdges.insert(i + 1, (triK, 0))
            elif fval == 7:
                crossEdges.pop(i)
                crossEdges.insert(i, (triK, 2))
                crossEdges.insert(i + 1, (triL, 0))
            elif fval == 8:
                crossEdges.pop(i)
                crossEdges.insert(i, (triL, 0))
                newEdgesToAnalize.append((triK, 2))
            _swapTriangles(plist2d, tri_vertex, tri_adjoin, triK, triL, vertex_in_tri)
            i = 0  # reinicio el contador a cero
    # ahora para cada uno de los vertices creados, analizo la condicion de delauney
    while newEdgesToAnalize:
        triK, pA = newEdgesToAnalize.pop(0)
        triL = tri_adjoin[triK][pA]
        pP = (pA - 1) % 3
        pOP = (tri_vertex[triL].index(tri_vertex[triK][pA]) + 1) % 3
        pB = (pA + 1) % 3
        if not shouldSwap(plist2d[tri_vertex[triK][pP]],
                          plist2d[tri_vertex[triK][pA]],
                          plist2d[tri_vertex[triK][pB]],
                          plist2d[tri_vertex[triL][pOP]]):
            continue
        _swapTriangles(plist2d, tri_vertex, tri_adjoin, triK, triL, vertex_in_tri)
    return



cdef (int * , int *) delauney_2d_constricted(double * plist2d, int num_point,
                                             int * contour, int num_contour):
    """
    Function that construct a 2d delauney triangulation with constrictions
    :param plist2d: points 2d 
    :param num_point: num points in array
    :param contour: index of points that create a contours. The index -1 divide the contour. The contour should be closed (last index == first)
    :param num_contour: num of points in the contour array
    """
    cdef int i, j, k
    cdef int n = num_point + 3
    cdef int num_tri = 2 * n - 2 - 3  # esquema de triangulo envolvente
    cdef double * plist2d_ext = <double *> PyMem_Malloc(n * 2 * sizeof(double))
    cdef int * avoid_swap = <int *> PyMem_Malloc(n * n * sizeof(int))
    cdef int * used_vertex = <int *> PyMem_Malloc(num_point * sizeof(int))
    cdef int * tri_vertex = <int *> PyMem_Malloc (num_tri*3*sizeof(int))
    cdef int * tri_adjoin = <int *> PyMem_Malloc (num_tri*3*sizeof(int))
    cdef int * vertex_in_tri = <int *> PyMem_Malloc(num_point * sizeof(int))
    cdef int last_tri_index = 0
    cdef int last_adj_index = 0


    try:
        if plist2d_ext == NULL or avoid_swap == NULL or used_vertex == NULL or tri_vertex == NULL or tri_adjoin == NULL or vertex_in_tri == NULL :
            raise MemoryError()

        for i in range(2*num_point):
            plist2d_ext[i] = plist2d[i]
        plist2d_ext[2*num_point] = -100.0
        plist2d_ext[2*num_point+1] = -100.0
        plist2d_ext[2*num_point+2] = 100.0
        plist2d_ext[2*num_point+3] = -100.0
        plist2d_ext[2*num_point+4] = 0.0
        plist2d_ext[2*num_point+5] = 100.0
        # inicializo avoid_swap

        for i in range(n*n):
            avoid_swap[i] = 0
        for i in range(num_point, n):
            for j in range(num_point, n):
                avoid_swap[n*i+j] = 1
        for i in range(num_contour - 1):
            if contour[i] == -1 or contour[i + 1] == -1:
                continue
            avoid_swap[contour[i]*n + contour[i+1]] = 1
            avoid_swap[contour[i+1]*n + contour[i]] = 1
        for i in range(num_point):
            used_vertex[i] = 0
            vertex_in_tri[i] = 0
        for i in range(num_tri*3):
            tri_vertex[i] = -2
            tri_adjoin[i] = -2
        tri_vertex[0] = num_point
        tri_vertex[1] = num_point+1
        tri_vertex[2] = num_point+2 # el triangulo principal
        tri_adjoin[0] = -1
        tri_adjoin[1] = -1
        tri_adjoin[2] = -1 # periferia



        # agrego primero todos los puntos del contorno, y mantengo la linea recta.
        i = 0
        while i < num_contour:
            if contour[i] == -1:
                i += 1
                continue
            if used_vertex[contour[i]] == 0:  # el punto no se habia usado antes. Lo agrego
                newTriIndex = split_triangle_into_3(plist2d_ext, tri_vertex, tri_adjoin, vertex_in_tri,
                                                  vertex_in_tri[contour[i]],
                                                  contour[i], lastTriIndex)
                used_vertex[contour[i]] = 1
                lastTriIndex += 2
                for triP in newTriIndex:
                    _checkAndSwap(plist2d, tri_vertex, tri_adjoin, vertex_in_tri, triP, avoid_swap)
            if i > 0 and contour[i - 1] != -1:
                forceEdge(plist2d, tri_vertex, tri_adjoin, contour[i - 1], contour[i], vertex_in_tri)
            i += 1

        # agrego el resto de los puntos
        for i in range(n):
            if used_vertex[i]:
                continue
            newTriIndex = _splitTriangleInto3(plist2d, tri_vertex, tri_adjoin, vertex_in_tri, vertex_in_tri[i], i,
                                              lastTriIndex)
            used_vertex[i] = 1
            lastTriIndex += 2
            for triP in newTriIndex:
                _checkAndSwap(plist2d, tri_vertex, tri_adjoin, vertex_in_tri, triP,
                              avoid_swap)  # siempre el valor P es el tercer dato

        # busco los triangulos que son a eliminar. Esto solo sera recursivo para los
        # elimino primero todos lost triangulos que esten en contacto con los extremos
        for i in reversed(range(ntri)):
            if tri_vertex[i]:
                break
        if ntri != i + 1:
            raise RuntimeError("This solution has not the same number of triangles as spected")
        triangles_to_pop = []
        triangles_to_keep = []
        if contour:
            for i in range(ntri):
                if any([x >= n for x in
                        tri_vertex[i]]):
                    triangles_to_pop.append(i)
                    _inspectGroupedTriangles(tri_vertex, tri_adjoin, avoid_swap, triangles_to_pop, i)
                    break
            triangles_to_inspect = set(range(ntri)) - set(triangles_to_pop)
            grouped_triangles = []
            while triangles_to_inspect:
                tri_group = [triangles_to_inspect.pop()]
                _inspectGroupedTriangles(tri_vertex, tri_adjoin, avoid_swap, tri_group, tri_group[0])
                triangles_to_inspect.difference_update(tri_group)
                grouped_triangles.append(tri_group)
            groups_index = [0] * len(grouped_triangles)  # 0 =unknown, 1=keep, -1=remove
            is_to_keep = True
            while any([x == 0 for x in groups_index]):
                for i in range(len(groups_index)):
                    if groups_index[i] != 0:
                        continue
                    founded = False
                    if is_to_keep:
                        for tri in grouped_triangles[i]:
                            for j in range(3):
                                if tri_adjoin[tri][j] in triangles_to_pop:
                                    founded = True
                                    break
                            if founded:
                                break
                        if founded:
                            triangles_to_keep.extend(grouped_triangles[i])
                            groups_index[i] = 1
                    else:
                        for tri in grouped_triangles[i]:
                            for j in range(3):
                                if tri_adjoin[tri][j] in triangles_to_keep:
                                    founded = True
                                    break
                            if founded:
                                break
                        if founded:
                            triangles_to_pop.extend(grouped_triangles[i])
                            groups_index[i] = -1
                is_to_keep = not is_to_keep
        else:
            triangles_to_keep = list(set(range(ntri)) - set(triangles_to_pop))

        # ya tengo los triangulos a guardar y a borrar
        # los indices de tri_index no cambian, pero al eliminar triangulos los indices de tri_adjoin si cambian
        index_to_discount = [0] * ntri
        triangles_to_pop = list(sorted(triangles_to_pop))
        triangles_to_pop.append(0)
        j = 0
        for i in range(ntri):
            if i == triangles_to_pop[j]:
                j += 1
                index_to_discount[i] = (i + 1)
            else:
                index_to_discount[i] = j
        returned_tri_vertex = []
        returned_tri_adjoin = []
        for i in range(ntri):
            if i - index_to_discount[i] < 0:  # este triangulo no existe
                continue
            returned_tri_vertex.append(tri_vertex[i])
            returned_tri_adjoin.append([x - index_to_discount[x] for x in tri_adjoin[i]])

        return returned_tri_vertex, returned_tri_adjoin

    finally:
        if avoid_swap != NULL:
            PyMem_Free(avoid_swap)
            avoid_swap = NULL
        if used_vertex != NULL:
            PyMem_Free(used_vertex)
            used_vertex = NULL
        if tri_vertex != NULL:
            PyMem_Free(tri_vertex)
            tri_vertex = NULL
        if tri_adjoin != NULL:
            PyMem_Free(tri_adjoin)
            tri_adjoin = NULL
        if plist2d_ext != NULL:
            PyMem_Free(plist2d_ext)
            plist2d_ext = NULL
        if vertex_in_tri != NULL:
            PyMem_Free(vertex_in_tri)
            vertex_in_tri = NULL
