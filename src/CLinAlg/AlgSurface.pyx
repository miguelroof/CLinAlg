import sys
from cpython.mem cimport PyMem_Malloc, PyMem_Free, PyMem_Realloc
from libc.math cimport M_PI, fabs
from .AlgTool cimport presition
from .cimport AlgVector, AlgWire, AlgSegment, AlgMatrix, AlgLine, AlgQuaternion

cdef list tessellate(double * surfPoint, int * indPer, unsigned int numPer, double * axis):
    """
    Tessellate a surface into triangles. This functions dont verify the dimmensions for the triangle
    :param surfPoint: pointer to array that keeps points in rows
    :param indPer: int pointer to the index for the points of perimeters. Each edge is separated by -1. The edges cannot be ordered
    :param numPer: number of points of perimeters (size of indPer)
    :return: list of index of points for triangles.
    """
    # toda la algebra de la funcion raiz se definia a partir de un contorno CERRADO
    cdef int i, j, k, ci, pi, ni, npoints
    cdef double modulo
    cdef double * ndir = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * pdir = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * vtemp = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef int * itemp = <int *> PyMem_Malloc(numPer * sizeof(int))
    cdef list triangles = []
    try:
        npoints = numPer
        for j in range(npoints):
            itemp[j] = indPer[j]
        i = 0
        indices = [indPer[j] for j in range(numPer - 1)]
        while npoints > 4 and i < npoints - 1:
            pi = itemp[npoints - 2] if i == 0 else itemp[i - 1]
            ci = itemp[i]
            ni = itemp[i + 1]
            AlgVector.vdir(pdir, &surfPoint[pi * 3], &surfPoint[ci * 3])
            AlgVector.vdir(ndir, &surfPoint[ci * 3], &surfPoint[ni * 3])
            AlgVector.cross(vtemp, pdir, ndir)
            modulo = AlgVector.module(vtemp)
            if modulo > presition:
                for j in range(3):
                    vtemp[j] /= modulo
            else:
                i += 1
                continue
            AlgVector.add(vtemp, vtemp, axis)
            modulo = AlgVector.module(vtemp)
            if modulo < 0.5:
                i += 1
                continue
            AlgVector.vdir(vtemp, &surfPoint[pi * 3], &surfPoint[ni * 3])
            linecortes = AlgWire.getIntersectionPointWithLine(surfPoint, itemp, npoints, &surfPoint[ci * 3], vtemp,
                                                              <bint> False)
            cortes = []
            valid = True
            for j in range(len(linecortes)):
                if AlgSegment.isInside(&surfPoint[pi * 3], &surfPoint[ni * 3],
                                       (<AlgVector.Vector> (linecortes[i]))._v, <bint> False):
                    valid = False
                    break
            if not valid:
                i += 1
                continue
            triangles.extend([pi, ci, ni])
            k = 0
            for j in range(npoints):
                if j == i:
                    continue
                itemp[k] = itemp[j]
                k += 1
            npoints -= 1
            itemp[npoints - 1] = itemp[0]
            i = 0
        if npoints > 4:
            raise ValueError("Not possible to tessellate surface")
        triangles.extend([itemp[0], itemp[1], itemp[2]])
        return triangles
    finally:
        PyMem_Free(ndir)
        PyMem_Free(pdir)
        PyMem_Free(vtemp)
        PyMem_Free(itemp)

cdef void sortTrianglesByNormal(double * surfPoint, int * indTri, int numTri, double * direction):
    "Function that sort triangles for a given direction"
    cdef int i, ptemp
    cdef double mult
    cdef double * v1 = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * v2 = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * v3 = <double *> PyMem_Malloc(3 * sizeof(double))
    try:
        for i in range(numTri):
            AlgVector.sub(v1, &surfPoint[indTri[i * 3 + 1] * 3], &surfPoint[indTri[i * 3] * 3])
            AlgVector.sub(v2, &surfPoint[indTri[i * 3 + 2] * 3], &surfPoint[indTri[i * 3] * 3])
            AlgVector.cross(v3, v1, v2)
            mult = AlgVector.dot(direction, v3)
            if mult < 0:
                ptemp = indTri[i * 3 + 1]
                indTri[i * 3 + 1] = indTri[i * 3 + 2]
                indTri[i * 3 + 2] = ptemp
    finally:
        PyMem_Free(v1)
        PyMem_Free(v2)
        PyMem_Free(v3)

cdef double getArea(Surface surf):
    """
    Returns area for the given surface
    :param surf: Surface object
    :return: double (area)
    """
    cdef double area = 0.0
    cdef double * vAB = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * vAC = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * vCross = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef int i
    try:
        for i in range(0, surf.indTri._ntriangle * 3, 3):
            AlgVector.sub(vAB, &surf.vectMat._m[surf.indTri._m[i + 1] * 3], &surf.vectMat._m[surf.indTri._m[i] * 3])
            AlgVector.sub(vAC, &surf.vectMat._m[surf.indTri._m[i + 2] * 3], &surf.vectMat._m[surf.indTri._m[i] * 3])
            AlgVector.cross(vCross, vAB, vAC)
            area += AlgVector.module(vCross) / 2
        return area
    finally:
        PyMem_Free(vAB)
        PyMem_Free(vAC)
        PyMem_Free(vCross)

cdef void getCDG(double * cdg, Surface surf):
    """
    Set de cdg of given Surface
    :param cdg: pointer to the returned cdg
    :param surf: Surface object
    """
    cdef int i, j
    cdef double ox, oy, oz, area
    cdef double * vAB = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * vAC = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * vCross = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double totalArea
    totalArea = 0
    try:
        for i in range(3):
            cdg[i] = 0
        for i in range(0, surf.indTri._ntriangle * 3, 3):
            ox = sum([surf.vectMat._m[surf.indTri._m[i + k] * 3 + 0] for k in range(3)]) / 3.0
            oy = sum([surf.vectMat._m[surf.indTri._m[i + k] * 3 + 1] for k in range(3)]) / 3.0
            oz = sum([surf.vectMat._m[surf.indTri._m[i + k] * 3 + 2] for k in range(3)]) / 3.0
            AlgVector.sub(vAB, &surf.vectMat._m[surf.indTri._m[i + 1] * 3], &surf.vectMat._m[surf.indTri._m[i] * 3])
            AlgVector.sub(vAC, &surf.vectMat._m[surf.indTri._m[i + 2] * 3], &surf.vectMat._m[surf.indTri._m[i] * 3])
            AlgVector.cross(vCross, vAB, vAC)
            area = AlgVector.module(vCross) * 0.5
            cdg[0] += ox * area
            cdg[1] += oy * area
            cdg[2] += oz * area
            totalArea += area
        for i in range(3):
            cdg[i] = cdg[i] / totalArea
    finally:
        PyMem_Free(vAB)
        PyMem_Free(vAC)
        PyMem_Free(vCross)

cdef AlgMatrix.Matrix getInertia(Surface surf):
    """
    Returns inertia tensor of the given surface referenced to CDG in main axis
    :param surf: Surface object
    :return: Matrix 3x3 with the tensor (Ixx, Jxy, Jxz...)
    """
    cdef int i, j, k
    cdef double alfa
    cdef AlgMatrix.Matrix inertia = AlgMatrix.zeros(3, 3)
    cdef AlgMatrix.Matrix V = AlgMatrix.zeros(3, 3)
    cdef AlgMatrix.Matrix S = AlgMatrix.Matrix(
        [[2.0 / 24, 1.0 / 24, 1.0 / 24], [1.0 / 24, 2.0 / 24, 1.0 / 24], [1.0 / 24, 1.0 / 24, 2.0 / 24]])
    # cdef AlgMatrix.Matrix J = AlgMatrix.Matrix()
    # cdef AlgMatrix.Matrix IDD = AlgMatrix.identity(3)
    cdef AlgMatrix.Matrix C
    cdef double * v0
    cdef double * v1
    cdef double * v2
    cdef double * vtemp1 = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * vtemp2 = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * vtemp3 = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * pzero = <double *> PyMem_Malloc(3 * sizeof(double))
    try:
        pzero[0] = 0.0
        pzero[1] = 0.0
        pzero[2] = 0.0
        for i in range(0, surf.indTri._ntriangle * 3, 3):
            v0 = &surf.vectMat._m[surf.indTri._m[i] * 3]
            v1 = &surf.vectMat._m[surf.indTri._m[i + 1] * 3]
            v2 = &surf.vectMat._m[surf.indTri._m[i + 2] * 3]
            for j in range(3):
                for k in range(3):
                    V._m[j * 3 + k] = surf.vectMat._m[surf.indTri._m[i + j] * 3 + k]
            AlgVector.sub(vtemp1, v1, v0)
            AlgVector.sub(vtemp2, v2, v0)
            AlgVector.cross(vtemp3, vtemp1, vtemp2)
            alfa = AlgVector.module(vtemp3)
            C = alfa * (V.transpose() * (S * V))
            for j in range(3):
                vtemp1[j] = v0[j] + v1[j] + v2[j]
            alfa = (alfa / 24) * (
                    AlgVector.dot(v0, v0) + AlgVector.dot(v1, v1) + AlgVector.dot(v2, v2) + AlgVector.dot(vtemp1,
                                                                                                          vtemp1))
            for j in range(3):
                for k in range(3):
                    if j == k:
                        inertia._m[j * 3 + k] = inertia._m[j * 3 + k] + alfa - C._m[j * 3 + k]
                    else:
                        inertia._m[j * 3 + k] = inertia._m[j * 3 + k] - C._m[j * 3 + k]
        inertia = stenierInGlobalAxis_p(inertia, surf.area(), pzero, (<AlgVector.Vector> surf.cdg())._v)
        return inertia
    finally:
        PyMem_Free(vtemp1)
        PyMem_Free(vtemp2)
        PyMem_Free(vtemp3)
        PyMem_Free(pzero)

cdef double staticMoment(Surface surf, double * point, double * dir):
    """
    Returns static moment of the given surface to a line (point and direction)
    :param surf: Surface object
    :param point: pointer with the point within the line 
    :param dir: pointer to line direction
    :return: double with the static moment
    """
    cdef double mom = 0
    cdef double dist
    cdef double * vAB = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * vAC = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * vCross = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * cdg = <double *> PyMem_Malloc(3 * sizeof(double))
    try:
        for i in range(0, surf.indTri._ntriangle * 3, 3):
            for j in range(3):
                cdg[j] = sum([surf.vectMat._m[surf.indTri._m[i + k] * 3 + j] for k in range(3)]) / 3.0
            AlgVector.sub(vAB, &surf.vectMat._m[surf.indTri._m[i + 1] * 3], &surf.vectMat._m[surf.indTri._m[i] * 3])
            AlgVector.sub(vAC, &surf.vectMat._m[surf.indTri._m[i + 2] * 3], &surf.vectMat._m[surf.indTri._m[i] * 3])
            AlgVector.cross(vCross, vAB, vAC)
            area = AlgVector.module(vCross) * 0.5
            dist = AlgLine.distanceToPoint(point, dir, cdg)
            mom += area * dist
        return mom
    finally:
        PyMem_Free(vAB)
        PyMem_Free(vAC)
        PyMem_Free(vCross)
        PyMem_Free(cdg)

cdef bint isInside(Surface surf, double * point, bint incEdge):
    """
    Check if one point is inside the surface
    :param surf: Surface object
    :param point: pointer to the point to be checked
    :param incEdge: Include edge in check condition
    :return: True if the point is within the surface, otherwise false
    """
    cdef double * mm = <double *> PyMem_Malloc(9 * sizeof(double))
    cdef double * ang = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double lmod
    cdef int i, j
    try:
        indices = [[0, 1], [0, 2], [1, 2]]
        for i in range(0, len(surf.indTri), 3):
            for j in range(3):
                AlgVector.sub(&mm[j * 3], &surf.vectMat._m[surf.indTri._m[i + j] * 3], point)
                lmod = AlgVector.module(&mm[j * 3])
                if lmod < presition:  # puntos coincidentes.
                    return incEdge
            for j in range(3):
                ang[j] = AlgVector.angle(&mm[indices[j][0] * 3], &mm[indices[j][1] * 3])
            if fabs(ang[0] + ang[1] + ang[2] - 2 * M_PI) < presition:
                # puede que este en una linelinea o no
                if incEdge:
                    return True
                else:
                    # puede que la linea no sea de indices consecutivos
                    for j in range(3):
                        if fabs(ang[j] - M_PI) < presition:
                            # Se encuentra sobre la linea
                            if abs(surf.indTri._m[i + indices[j][0]] - surf.indTri._m[
                                i + indices[j][1]]) == 1:  # son dos puntos consecutivos del contorno
                                return False
                    return True
        return False  # llegado a este punto, si ningun trianglo lo contempla, que salga
    finally:
        PyMem_Free(mm)
        PyMem_Free(ang)

cdef tuple dictIntersectByLine(Surface surf, double * point, double * dir, bint incEdge):
    """
    Returns a tuple of values (point, segment) with all intersection points of a line with the surface
    :param surf: Surface object 
    :param point: pointer with the point within the line
    :param dir: pointer with the direction of the line
    :param incEdge: Include edge in check condition
    :return: Tuple(list(index of edge), list(points of intersection))
    """
    cdef int i, n, pp, pc, pn
    cdef double * ipoint = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * tpoint = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef bint sclosed
    cdef double * dirMIn = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * dirMOut = <double *> PyMem_Malloc(3 * sizeof(double))
    wireInd = surf.indPer.__serialize__()
    n = <int> (len(wireInd))
    sclosed = True
    try:
        cortes = []
        indices = []
        for i in range(n - 1):
            pc = i
            pn = i + 1
            if not AlgSegment.getIntersectionPointWithLine(ipoint, &surf.vectMat._m[pc * 3],
                                                           &surf.vectMat._m[pn * 3],
                                                           point,
                                                           dir,
                                                           incEdge) or any(
                [AlgVector.isEqual(ipoint, (<AlgVector.Vector> v)._v) for v in cortes]):
                continue
            if (incEdge or (not AlgVector.isEqual(ipoint, &surf.vectMat._m[pc * 3])
                            and not AlgVector.isEqual(ipoint, &surf.vectMat._m[pn * 3]))):
                newVect = AlgVector.Vector(None)
                (<AlgVector.Vector> newVect)._v = ipoint
                cortes.append(newVect)
                indices.append(i)
                ipoint = <double *> PyMem_Malloc(3 * sizeof(double))
                continue
            if AlgVector.isEqual(ipoint, &surf.vectMat._m[pc * 3]):  # punto inicial
                if i == 0:
                    if not sclosed:
                        continue
                    pp = n - 2
                    AlgVector.vdir(dirMIn, &surf.vectMat._m[pc * 3], &surf.vectMat._m[pp * 3])
                else:
                    pp = <int> (i - 1)
                    AlgVector.vdir(dirMIn, &surf.vectMat._m[pc * 3], &surf.vectMat._m[pp * 3])
                AlgVector.vdir(dirMOut, &surf.vectMat._m[pc * 3], &surf.vectMat._m[pn * 3])
                AlgVector.add(dirMIn, &surf.vectMat._m[pc * 3], dirMIn)
                AlgVector.add(dirMOut, &surf.vectMat._m[pc * 3], dirMOut)
            elif AlgVector.isEqual(ipoint, &surf.vectMat._m[pn * 3]):  #punto final
                if i == n - 2:
                    if not sclosed:
                        continue
                    pp = 1
                    AlgVector.vdir(dirMOut, &surf.vectMat._m[pn * 3], &surf.vectMat._m[pp * 3])
                else:
                    pp = i + 2
                    AlgVector.vdir(dirMOut, &surf.vectMat._m[pn * 3], &surf.vectMat._m[pp * 3])
                AlgVector.vdir(dirMIn, &surf.vectMat._m[pn * 3], &surf.vectMat._m[pc * 3])
                AlgVector.add(dirMIn, &surf.vectMat._m[pn * 3], dirMIn)
                AlgVector.add(dirMOut, &surf.vectMat._m[pn * 3], dirMOut)
            else:  #ninguna de las anteriores
                raise ValueError('No es posible esta combinacion')
            if AlgSegment.getIntersectionPointWithLine(tpoint, dirMIn, dirMOut, point,
                                                       dir, <bint> False):
                newVect = AlgVector.Vector(None)
                (<AlgVector.Vector> newVect)._v = ipoint
                cortes.append(newVect)
                indices.append(i)
                ipoint = <double *> PyMem_Malloc(3 * sizeof(double))
                continue
        return indices, cortes
    finally:
        PyMem_Free(dirMIn)
        PyMem_Free(dirMOut)
        PyMem_Free(tpoint)
        PyMem_Free(ipoint)

cdef int _splitTriangleByLine(double * v1, double * v2, double * wirePoint, int * indexTriangle, double * linepnt,
                              double * linedir):
    cdef int i, j
    #1 = 0-1, 2 = 0-2, 4 = 1-2. Devuelvo un valor sumando de los cortes, que sera 3, 5 o 6. En caso contrario, no se han producido dos cortes.
    cdef int cut1 = 0
    cdef int cut2 = 0
    cdef AlgSegment.PhantomSegment seg

    if AlgSegment.getIntersectionPointWithLine(v1, &wirePoint[indexTriangle[0] * 3], &wirePoint[indexTriangle[1] * 3],
                                               linepnt, linedir, <bint> True):  # corte a-b
        cut1 = 1
    if AlgSegment.getIntersectionPointWithLine(v2, &wirePoint[indexTriangle[0] * 3], &wirePoint[indexTriangle[2] * 3],
                                               linepnt, linedir, <bint> True):  # corte a-c
        cut2 = 2
    if not cut1 and AlgSegment.getIntersectionPointWithLine(v1, &wirePoint[indexTriangle[1] * 3],
                                                            &wirePoint[indexTriangle[2] * 3], linepnt, linedir,
                                                            <bint> True):  # corte b-c
        cut1 = 4
    elif not cut2 and AlgSegment.getIntersectionPointWithLine(v2, &wirePoint[indexTriangle[1] * 3],
                                                              &wirePoint[indexTriangle[2] * 3], linepnt, linedir,
                                                              <bint> True):  # corte b-c
        cut2 = 4
    if cut1 + cut2 in (3, 5, 6):  # tengo que verificar si los dos puntos de corte se encuentran sobre el mismo segmento
        for i in range(3):
            if AlgSegment.isInside(&wirePoint[indexTriangle[i] * 3], &wirePoint[indexTriangle[(i + 1) % 3] * 3], v1,
                                   <bint> True) and AlgSegment.isInside(&wirePoint[indexTriangle[i] * 3],
                                                                        &wirePoint[indexTriangle[(i + 1) % 3] * 3], v2,
                                                                        <bint> True):
                cut1 = cut2 = 0
                break
    return cut1 + cut2

cdef list _getEdgeByTriangleIndex(list indTri):
    "Devuelve los edges de una superficie. No tienen porque estar ordenados"
    cdef int i
    cdef int pa, pb, pc
    cdef dict laterales = {}
    cdef list tupLat = []
    for i in range(0, <int> len(indTri), 3):
        pa, pb, pc = list(sorted([indTri[i + j] for j in range(3)]))
        laterales[(pa, pb)] = laterales.get((pa, pb), []) + [i]
        laterales[(pa, pc)] = laterales.get((pa, pc), []) + [i]
        laterales[(pb, pc)] = laterales.get((pb, pc), []) + [i]
    for lat, vind in laterales.items():
        if len(vind) > 1:
            continue
        tupLat.extend(list(lat))
        tupLat.append(-1)  # [1,2,-1,2,3,-1,2,4,-1]
    return tupLat

# cdef list _getPerimeterByTriangleIndex(list indTri):
#     """
#     Function that computes outer for the given points and triangles
#     :param surfpoints: double pointer with the points
#     :param indTri: int pointer with the index of triangles, given by 3 in 3
#     :returns list: list with the points that defines perimeteres
#     """
#     cdef int i, j, k, currTri
#     cdef int pa, pb, pc
#     cdef list outerIndex
#     cdef dict laterales = {}
#     cdef list tupLat
#     cdef int numTri = <int> len(indTri) // 3
#     cdef bint founded
#
#     for i in range(numTri):
#         pa, pb, pc = list(sorted([indTri[i * 3 + j] for j in range(3)]))
#         laterales[(pa, pb)] = laterales.get((pa, pb), []) + [i]
#         laterales[(pa, pc)] = laterales.get((pa, pc), []) + [i]
#         laterales[(pb, pc)] = laterales.get((pb, pc), []) + [i]
#
#     tupLat = []
#     for lat in list(laterales.keys()):
#         if len(laterales[lat]) > 1:
#             continue
#         tupLat.append(lat)
#     outerList = []
#     outerIndex = [tupLat[0][0], tupLat[0][1]]
#     tupLat.pop(0)
#     while tupLat:
#         founded = False
#         i = 0
#         while i < len(tupLat):
#             if tupLat[i][0] == outerIndex[-1]:
#                 outerIndex.append(tupLat[i][1])
#                 tupLat.pop(i)
#                 founded = True
#                 break
#             elif tupLat[i][1] == outerIndex[-1]:
#                 outerIndex.append(tupLat[i][0])
#                 tupLat.pop(i)
#                 founded = True
#                 break
#             i += 1
#
#         if not founded and tupLat:
#             if outerIndex[0] != outerIndex[-1]:
#                 raise ValueError("Not closed contour with this triangles indexes and tupLat value")
#             outerList.append(outerIndex)
#             outerIndex = [tupLat[0][0], tupLat[0][1]]
#             tupLat.pop(0)
#     if outerIndex[0] != outerIndex[-1]:
#         raise ValueError("Not closed contour with this triangles indexes")
#     outerList.append(outerIndex)
#     return outerList

cdef list _getSurfaceByTriangleIndex(list triIndex):
    cdef dict laterales = {}
    cdef dict triangulos = {}
    cdef int numTri = <int> len(triIndex) // 3
    cdef int pa, pb, pc, itri, i, j
    cdef list surfaces, surf, psurf, listTri
    cdef bint founded
    for i in range(numTri):
        pa, pb, pc = list(sorted([triIndex[i * 3 + j] for j in range(3)]))
        triangulos[i] = [(pa, pb), (pa, pc), (pb, pc)]
        laterales[(pa, pb)] = laterales.get((pa, pb), []) + [i]
        laterales[(pa, pc)] = laterales.get((pa, pc), []) + [i]
        laterales[(pb, pc)] = laterales.get((pb, pc), []) + [i]
    surfaces = []
    listTri = list(range(numTri))
    surf = [listTri.pop(0)]
    psurf = [triIndex[surf[0] * 3 + j] for j in range(3)]
    while True:
        i = 0
        while i < len(listTri):
            founded = False
            itri = listTri[i]
            for lat in triangulos[itri]:
                for stri in laterales[lat]:
                    if itri != stri and stri in surf:
                        surf.append(itri)
                        psurf.extend([triIndex[itri * 3 + j] for j in range(3)])
                        listTri.pop(i)
                        founded = True
                        break
                if founded:
                    break
            if not founded:
                i += 1
            else:
                i = 0
        surfaces.append(psurf)
        if listTri:
            surf = [listTri.pop(0)]
            psurf = [triIndex[surf[0] * 3 + j] for j in range(3)]
        else:
            break
    return surfaces

cdef list splitByLine(Surface surf, double * linepnt, double * linedir):
    "returns a list of surfaces returned from split"
    cdef double * v1 = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * v2 = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * normal = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * subnormal = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * vpos = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * vtemp1 = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef int vcorte
    cdef list pointsLeft = []
    cdef list triIndexLeft = []
    cdef list pointsRight = []
    cdef list triIndexRight = []
    cdef list surfaces = []
    cdef dict replDict
    cdef int i, j, pa, pb, pc
    cdef AlgVector.PhantomVector phv0, phv1, phv2, phc0, phc1
    cdef double vposit
    cdef dict laterales
    cdef Surface newSurf
    cdef list tsurf, surfacesIndex, perimIndex, tPerimIndex, tindex, tPerim
    cdef int minIndex
    cdef set pdist
    cdef AlgWire.IndexPerimeter indPer = AlgWire.IndexPerimeter()
    try:
        # la normal sera la noral de uno de los triangulos
        AlgVector.sub(v1, &surf.vectMat._m[surf.indTri._m[1] * 3], &surf.vectMat._m[surf.indTri.m[0] * 3])
        AlgVector.sub(v2, &surf.vectMat._m[surf.indTri._m[2] * 3], &surf.vectMat._m[surf.indTri.m[0] * 3])
        AlgVector.cross(normal, v1, v2)
        vposit = AlgVector.module(normal)
        for i in range(3):
            normal[i] = normal[i] / vposit
        subnormal[0] = -normal[0]
        subnormal[1] = -normal[1]
        subnormal[2] = -normal[2]
        AlgVector.cross(vpos, normal, linedir)
        for itri in range(0, surf.indTri._ntriangle * 3, 3):
            vcorte = _splitTriangleByLine(v1, v2, surf.vectMat._m, &surf.indTri._m[itri], linepnt, linedir)
            if vcorte == 3:  # corte entre el los vectores 0-1 (v1) y 0-2 (v2). El vector aislado es el cero
                pindep = 0
                phv0 = AlgVector.PhantomVector()
                phv0._v = &surf.vectMat._m[surf.indTri._m[itri] * 3]
                phv1 = AlgVector.PhantomVector()
                phv1._v = &surf.vectMat._m[surf.indTri._m[itri + 1] * 3]
                phv2 = AlgVector.PhantomVector()
                phv2._v = &surf.vectMat._m[surf.indTri._m[itri + 2] * 3]
                phc0 = AlgVector.PhantomVector()
                phc0._v = v1
                phc1 = AlgVector.PhantomVector()
                phc1._v = v2
            elif vcorte == 5:  # corte entre los vectores 0-2 (v2) y 1-2 (v1). El vector aislado es el 2
                pindep = 1
                phv0 = AlgVector.PhantomVector()
                phv0._v = &surf.vectMat._m[surf.indTri._m[itri + 1] * 3]
                phv1 = AlgVector.PhantomVector()
                phv1._v = &surf.vectMat._m[surf.indTri._m[itri] * 3]
                phv2 = AlgVector.PhantomVector()
                phv2._v = &surf.vectMat._m[surf.indTri._m[itri + 2] * 3]
                phc0 = AlgVector.PhantomVector()
                phc0._v = v1
                phc1 = AlgVector.PhantomVector()
                phc1._v = v2
            elif vcorte == 6:  # corte entre los vectores 0-1 (v1) y 1-2 (v2). El vector aislado es el 1
                pindep = 2
                phv0 = AlgVector.PhantomVector()
                phv0._v = &surf.vectMat._m[surf.indTri._m[itri + 2] * 3]
                phv1 = AlgVector.PhantomVector()
                phv1._v = &surf.vectMat._m[surf.indTri._m[itri + 1] * 3]
                phv2 = AlgVector.PhantomVector()
                phv2._v = &surf.vectMat._m[surf.indTri._m[itri] * 3]
                phc0 = AlgVector.PhantomVector()
                phc0._v = v1
                phc1 = AlgVector.PhantomVector()
                phc1._v = v2
            else:
                phv0 = AlgVector.PhantomVector()
                phv0._v = &surf.vectMat._m[surf.indTri._m[itri] * 3]
                phv1 = AlgVector.PhantomVector()
                phv1._v = &surf.vectMat._m[surf.indTri._m[itri + 1] * 3]
                phv2 = AlgVector.PhantomVector()
                phv2._v = &surf.vectMat._m[surf.indTri._m[itri + 2] * 3]

            if vcorte in (3, 5, 6):
                if AlgLine.isInside(linepnt, linedir, &surf.vectMat._m[surf.indTri._m[itri + pindep] * 3]):
                    AlgVector.sub(vtemp1, &surf.vectMat._m[surf.indTri._m[itri + ((pindep + 1) % 3)] * 3], linepnt)
                    vposit = -AlgVector.dot(vtemp1, vpos)
                else:
                    AlgVector.sub(vtemp1, &surf.vectMat._m[surf.indTri._m[itri + pindep] * 3], linepnt)
                    vposit = AlgVector.dot(vtemp1, vpos)

                if vposit > 0:
                    # izquierdo
                    for ph in (phv0, phc0, phc1):
                        if ph in pointsLeft:
                            triIndexLeft.append(pointsLeft.index(ph))
                        else:
                            pointsLeft.append(ph.copy())
                            triIndexLeft.append(len(pointsLeft) - 1)
                    # derecho
                    if phc0 in (phv1, phv2):
                        for ph in (phc1, phv1, phv2):
                            if ph in pointsRight:
                                triIndexRight.append(pointsRight.index(ph))
                            else:
                                pointsRight.append(ph.copy())
                                triIndexRight.append(len(pointsRight) - 1)
                    elif phc1 in (phv1, phv2):
                        for ph in (phc0, phv1, phv2):
                            if ph in pointsRight:
                                triIndexRight.append(pointsRight.index(ph))
                            else:
                                pointsRight.append(ph.copy())
                                triIndexRight.append(len(pointsRight) - 1)
                    else:
                        for ph in (phc0, phv1, phv2, phc0, phv2, phc1):
                            if ph in pointsRight:
                                triIndexRight.append(pointsRight.index(ph))
                            else:
                                pointsRight.append(ph.copy())
                                triIndexRight.append(len(pointsRight) - 1)
                elif vposit < 0:
                    # izquierdo
                    if phc0 in (phv1, phv2):
                        for ph in (phc1, phv1, phv2):
                            if ph in pointsLeft:
                                triIndexLeft.append(pointsLeft.index(ph))
                            else:
                                pointsLeft.append(ph.copy())
                                triIndexLeft.append(len(pointsLeft) - 1)
                    elif phc1 in (phv1, phv2):
                        for ph in (phc0, phv1, phv2):
                            if ph in pointsLeft:
                                triIndexLeft.append(pointsLeft.index(ph))
                            else:
                                pointsLeft.append(ph.copy())
                                triIndexLeft.append(len(pointsLeft) - 1)
                    else:
                        for ph in (phc0, phv1, phv2, phc0, phv2, phc1):
                            if ph in pointsLeft:
                                triIndexLeft.append(pointsLeft.index(ph))
                            else:
                                pointsLeft.append(ph.copy())
                                triIndexLeft.append(len(pointsLeft) - 1)
                    # derecho
                    for ph in (phv0, phc0, phc1):
                        if ph in pointsRight:
                            triIndexRight.append(pointsRight.index(ph))
                        else:
                            pointsRight.append(ph.copy())
                            triIndexRight.append(len(pointsRight) - 1)
            else:
                for i in range(3):
                    vposit = 0
                    if not AlgLine.isInside(linepnt, linedir, &surf.vectMat._m[surf.indTri._m[itri + i] * 3]):
                        AlgVector.sub(vtemp1, &surf.vectMat._m[surf.indTri._m[itri + i] * 3], linepnt)
                        vposit = AlgVector.dot(vtemp1, vpos)
                        break
                if vposit > 0:  # todos los putnos estan a la izquierda
                    for ph in (phv0, phv1, phv2):
                        if ph in pointsLeft:
                            triIndexLeft.append(pointsLeft.index(ph))
                        else:
                            pointsLeft.append(ph.copy())
                            triIndexLeft.append(len(pointsLeft) - 1)
                elif vposit < 0:
                    for ph in (phv0, phv1, phv2):
                        if ph in pointsRight:
                            triIndexRight.append(pointsRight.index(ph))
                        else:
                            pointsRight.append(ph.copy())
                            triIndexRight.append(len(pointsRight) - 1)
        # leftside
        if pointsLeft:
            surfacesIndex = _getSurfaceByTriangleIndex(triIndexLeft)
            for tsurf in surfacesIndex:
                pdist = set(tsurf)
                minIndex = min(pdist)
                newSurf = Surface()
                newSurf.vectMat._m = <double *> PyMem_Malloc(3 * len(pdist) * sizeof(double))
                newSurf.vectMat._rows = <unsigned int> len(pdist)
                newSurf.vectMat._cols = 3
                replDict = {}
                i = 0
                for tindex in tsurf:
                    if replDict.get(tindex) is None:
                        for j in range(3):
                            newSurf.vectMat._m[i * 3 + j] = pointsLeft[tindex][j]
                        replDict[tindex] = i
                        i += 1
                newSurf.indTri.setList([replDict[i] for i in tsurf])
                tEdgeIndex = _getEdgeByTriangleIndex([replDict[i] for i in tsurf])  #[0,1,-1,1,2,-1...
                # perimIndex = []
                # tPerimIndex = _getPerimeterByTriangleIndex(tsurf)
                # tPerimIndex = list(
                #     zip(tPerimIndex, [min([pointsLeft[y][1] for y in tPerim]) for tPerim in tPerimIndex]))
                # tPerimIndex = [x[0] for x in sorted(tPerimIndex, key=_sortkey)]
                # i = 0
                # while tPerimIndex:
                #     tPerim = [replDict[j] for j in tPerimIndex.pop(0)]
                #     indPer.setList(tPerim)
                #     AlgWire.getNormal(vtemp1, newSurf.vectMat._m, indPer._m, indPer._npoint)
                #     if i == 0:
                #         if not AlgVector.isEqual(vtemp1, normal):
                #             tPerim = list(reversed(tPerim))
                #         i += 1
                #     else:
                #         if not AlgVector.isEqual(vtemp1, subnormal):
                #             tPerim = list(reversed(tPerim))
                #     perimIndex.extend(tPerim)
                #     if tPerimIndex:
                #         perimIndex.append(-1)
                # newSurf.indPer.setList(perimIndex)
                newSurf.indPer.setList(tEdgeIndex)
                sortTrianglesByNormal(newSurf.vectMat._m, newSurf.indTri._m, newSurf.indTri._ntriangle, normal)
                surfaces.append(newSurf)
        # rightside
        if pointsRight:
            surfacesIndex = _getSurfaceByTriangleIndex(triIndexRight)
            for tsurf in surfacesIndex:
                pdist = set(tsurf)
                minIndex = min(pdist)
                newSurf = Surface()
                newSurf.vectMat._m = <double *> PyMem_Malloc(3 * len(pdist) * sizeof(double))
                newSurf.vectMat._rows = <unsigned int> len(pdist)
                newSurf.vectMat._cols = 3
                replDict = {}
                i = 0
                for tindex in tsurf:
                    if replDict.get(tindex) is None:
                        for j in range(3):
                            newSurf.vectMat._m[i * 3 + j] = pointsRight[tindex][j]
                        replDict[tindex] = i
                        i += 1
                newSurf.indTri.setList([replDict[i] for i in tsurf])
                sortTrianglesByNormal(newSurf.vectMat._m, newSurf.indTri._m, newSurf.indTri._ntriangle, normal)
                tEdgeIndex = _getEdgeByTriangleIndex([replDict[i] for i in tsurf])  #[0,1,-1,1,2,-1...
                newSurf.indPer.setList(tEdgeIndex)
                # perimIndex = []
                # tPerimIndex = _getPerimeterByTriangleIndex(tsurf)
                # tPerimIndex = list(
                #     zip(tPerimIndex, [min([pointsRight[y][1] for y in tPerim]) for tPerim in tPerimIndex]))
                # tPerimIndex = [x[0] for x in sorted(tPerimIndex, key=_sortkey)]
                # i = 0
                # while tPerimIndex:
                #     tPerim = [replDict[j] for j in tPerimIndex.pop(0)]
                #     indPer.setList(tPerim)
                #     AlgWire.getNormal(vtemp1, newSurf.vectMat._m, indPer._m, indPer._npoint)
                #     if i == 0:
                #         if not AlgVector.isEqual(vtemp1, normal):
                #             tPerim = list(reversed(tPerim))
                #         i += 1
                #     else:
                #         if not AlgVector.isEqual(vtemp1, subnormal):
                #             tPerim = list(reversed(tPerim))
                #     perimIndex.extend(tPerim)
                #     if tPerimIndex:
                #         perimIndex.append(-1)
                # newSurf.indPer.setList(perimIndex)
                surfaces.append(newSurf)
        return surfaces
    finally:
        PyMem_Free(v1)
        PyMem_Free(v2)
        PyMem_Free(normal)
        PyMem_Free(subnormal)
        PyMem_Free(vpos)
        PyMem_Free(vtemp1)

cdef (double, double, double, double) boundRectangle(Surface surf, int v1, int v2):
    """
    Calculate boundRectangle for a surface
    :param surf: Surface object
    :param v1: Integer that defines first direction (x=0, y=1, z=2)
    :param v2: Integer that defines second direction (x=0, y=1, z=2)
    :return: tuple with (minV1, maxV1, minV2, maxV2)
    """
    cdef int i
    cdef double minV1, maxV1, minV2, maxV2
    minV1 = sys.float_info.max
    maxV1 = -sys.float_info.max
    minV2 = sys.float_info.max
    maxV2 = -sys.float_info.max
    for i in range(surf.vectMat._rows - 1):
        minV1 = min(minV1, surf.vectMat._m[i * 3 + v1])
        maxV1 = max(maxV1, surf.vectMat._m[i * 3 + v1])
        minV2 = min(minV2, surf.vectMat._m[i * 3 + v2])
        maxV2 = max(maxV2, surf.vectMat._m[i * 3 + v2])
    return (minV1, maxV1, minV2, maxV2)

cdef Matrix stenierInGlobalAxis_p(Matrix I, double area, double * cdg, double * to_point):
    """Apply Steniner in global Axis
       PARAMETERS:
       -----------
       arg1: I == tensor of inertia
           numpy array 3x3 with tensor of inertia given at 'cdg' Vector in Global coordinates
       arg2: Area
           int or float with the area of surface
       arg2: Vector with the CDG
           represent the CDG of surface that has the tensor given
       arg3: Vector
           Point to translate the inertie tensor
       RETURN:
       ----------
       numpy array 3x3 with the transformed tensor of inertia
       """
    cdef double modulo
    cdef AlgMatrix.Matrix r = AlgMatrix.Matrix()
    cdef AlgMatrix.Matrix mident = AlgMatrix.identity(3)
    cdef Matrix J
    cdef double * vdir = <double *> PyMem_Malloc(3 * sizeof(double))
    try:
        r._rows = 3
        r._cols = 3
        r._m = <double *> PyMem_Malloc(9 * sizeof(double))
        AlgVector.sub(vdir, cdg, to_point)
        modulo = AlgVector.dot(vdir, vdir)
        r._m[0] = vdir[0] * vdir[0]
        r._m[1] = vdir[0] * vdir[1]
        r._m[2] = vdir[0] * vdir[2]
        r._m[3] = vdir[1] * vdir[0]
        r._m[4] = vdir[1] * vdir[1]
        r._m[5] = vdir[1] * vdir[2]
        r._m[6] = vdir[2] * vdir[0]
        r._m[7] = vdir[2] * vdir[1]
        r._m[8] = vdir[2] * vdir[2]
        mident._m[0] = modulo
        mident._m[4] = modulo
        mident._m[8] = modulo
        mident = area * (mident - r)
        J = I - mident
        return J
    finally:
        PyMem_Free(vdir)

cdef int _sign(double v):
    return 1 if v >= -presition else -1

cdef list _getP(list trap, bint sup):
    cdef double _b, _e
    if sup:
        _b = trap[1]
        _e = trap[3]
    else:
        _b = trap[2]
        _e = trap[4]
    if _b < presition:
        return [_e]
    if trap[5]:
        if abs(abs(_e) - 0.5 * _b) > presition:
            return [-abs(_e) - 0.5 * _b, -abs(_e) + 0.5 * _b, abs(_e) - 0.5 * _b, abs(_e) + 0.5 * _b]
        else:
            return [-abs(_e) - 0.5 * _b, 0, abs(_e) + 0.5 * _b]
    else:
        return [_e - 0.5 * _b, _e + 0.5 * _b]

cdef double _distanceYZ(double y1, double z1, double y2, double z2):
    return ((y1 - y2) ** 2 + (z1 - z2) ** 2) ** 0.5

cdef double _areaYZ(double y1, double z1, double y2, double z2, double y3, double z3):
    return abs(0.5 * (y1 * (z2 - z3) + y2 * (z3 - z1) + y3 * (z1 - z2)))

cdef double _radiusYZ(double y1, double z1, double y2, double z2, double y3, double z3):
    "Function that computes the raduis for circumference that cross the 3 points"
    ar = _areaYZ(y1, z1, y2, z2, y3, z3)
    if ar:
        return _distanceYZ(y1, z1, y2, z2) * _distanceYZ(y2, z2, y3, z3) * _distanceYZ(y3, z3, y1, z1) / (4 * ar)
    else:
        return sys.float_info.max

cdef list _triangularize(double tn, double tp, list tlist, double bn, double bp, list blist, double h, int nstart):
    cdef int itn, itp, ibn, ibp, i, lentList
    cdef list triIndex = []
    cdef double r1, r2

    itn = itp = -1
    ibn = ibp = -1
    for i in range(len(tlist)):
        if itn == -1 and abs(tn - tlist[i]) < presition:
            itn = i
        if abs(tp - tlist[i]) < presition:
            itp = i
            break
    for i in range(len(blist)):
        if ibn == -1 and abs(bn - blist[i]) < presition:
            ibn = i
        if abs(bp - blist[i]) < presition:
            ibp = i
            break
    if ibp == -1 or itp == -1:
        raise ArithmeticError("Not founded limits in algorithm")
    triIndex = []
    lentlist = len(tlist)
    while itn < itp or ibn < ibp:
        if itn < itp and ibn < ibp:
            r1 = _radiusYZ(tlist[itn], h, tlist[itn + 1], h, blist[ibn], 0)
            r2 = _radiusYZ(tlist[itn], h, blist[ibn], 0, blist[ibn + 1], 0)
            if r1 <= r2:
                triIndex.extend([itn + nstart, itn + 1 + nstart, nstart + lentlist + ibn])
                itn += 1
            else:
                triIndex.extend([nstart + itn, nstart + lentlist + ibn, nstart + lentlist + ibn + 1])
                ibn += 1
        elif itn < itp:
            triIndex.extend([itn + nstart, itn + 1 + nstart, nstart + lentlist + ibn])
            itn += 1
        elif ibn < ibp:
            triIndex.extend([nstart + itn, nstart + lentlist + ibn, nstart + lentlist + ibn + 1])
            ibn += 1
    return triIndex

cdef double _sortkey(tuple v):
    return v[1]

cpdef list trapezeToSurfaceYZ(list trapeces):
    "Function that convert trapeze list [[h 0, bs 1, bi 2, es 3, ei 4, sym 5]... ] into surface object. The trapezes will be defined from upper to lower"
    # primero chequeo todos los trapecios y si hubiera algun simetrico con incompatibilidad de desarrollo, lo divido por la mitad
    cdef int i, npoints, nstart, tindex
    cdef double h
    cdef list trapMod, hpoints, tpoints, llist, triIndex, trap, pointList
    cdef list surfaces, perimIndex, tPerimIndex, tsurf, surfacesIndex, tPerim
    cdef dict replDict
    cdef double * vtemp1 = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * normal = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * subnormal = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef AlgWire.IndexPerimeter indPer = AlgWire.IndexPerimeter()
    cdef Surface newSurf
    if not vtemp1 or not normal or not subnormal:
        raise MemoryError()
    try:
        i = 0
        h = 0.0
        trapMod = []
        normal[0] = 1
        normal[1] = 0
        normal[2] = 0
        subnormal[0] = -1
        subnormal[1] = 0
        subnormal[2] = 0
        for trap in trapeces:
            if trap[5]:
                ts1 = trap[3] + 0.5 * trap[1]  # linea exterior f1 = ts1-(x/h)*(ts1-ti1)
                ti1 = trap[4] + 0.5 * trap[2]
                ts2 = trap[3] - 0.5 * trap[1]  # linea interior f2 = ts2-(x/h)*(ts2-ti2)
                ti2 = trap[4] - 0.5 * trap[1]
                # evaluo el corte mas superior
                if _sign(ts2) == _sign(ti2):
                    trapMod.append(trap[:])
                    continue
                xa = ts2 * trap[0] / (ts2 - ti2)
                tsa = ts1 - (xa / trap[0]) * (ts1 - ti1)
                trapMod.append([xa, abs(ts1 - ts2), abs(tsa), abs(0.5 * (ts1 + ts2)), abs(0.5 * tsa), 1])
                # evaluo ahora el corte intermedio
                if abs((ts1 - ti1) + (ts2 - ti2)) > presition:
                    xb = trap[0] * (ts1 + ts2) / ((ts1 - ti1) + (ts2 - ti2))
                    if xa < xb <= (trap[0] + presition):
                        tsb = ts1 - (xb / trap[0]) * (ts1 - ti1)
                        trapMod.append([xb - xa, abs(2 * tsa), abs(2 * tsb), 0, 0, 0])
                    else:
                        trapMod.append([trap[0] - xa, abs(2 * tsa), abs(2 * ti1), 0, 0, 0])
                        continue
                else:
                    trapMod.append([trap[0] - xa, abs(2 * tsa), abs(2 * ti1), 0, 0, 0])
                    continue
                # evaluo corte inferior
                if _sign(ts1) == _sign(ti1):
                    trapMod.append([trap[0] - xb, abs(tsb * 2), abs(2 * ti2), 0, 0, 0])
                    continue
                xc = ts1 * trap[0] / (ts1 - ti1)
                tsc = ts2 - (xc / trap[0]) * (ts2 - ti2)
                trapMod.append([xc - xb, abs(2 * tsb), abs(2 * tsc), 0, 0, 0])
                trapMod.append([trap[0] - xc, abs(tsc), abs(ti1 - ti2), 0.5 * abs(tsc), abs(0.5 * (ti1 + ti2)), 1])
            else:
                trapMod.append(trap[:])

        hpoints = [sum([x[0] for x in trapMod])]
        tpoints = [_getP(trapMod[0], True)]
        npoints = <int> len(tpoints[0])

        for i in range(<int> (len(trapMod)) - 1):
            llist = list(sorted(_getP(trapMod[i], False) + _getP(trapMod[i + 1], True)))
            j = 0
            while j < len(llist) - 1:
                if abs(llist[j] - llist[j + 1]) < presition:
                    llist.pop(j)
                    continue
                j += 1
            tpoints.append(llist)
            hpoints.append(hpoints[-1] - trapMod[i][0])
            npoints += <int> len(llist)

        hpoints.append(hpoints[-1] - trapMod[-1][0])
        tpoints.append(_getP(trapMod[-1], False))
        triIndex = []
        nstart = 0
        for i, trap in enumerate(trapMod):
            if trap[5]:
                triIndex.extend(_triangularize(-abs(trap[3]) - 0.5 * trap[1], -abs(trap[3]) + 0.5 * trap[1], tpoints[i],
                                               -abs(trap[4]) - 0.5 * trap[2], -abs(trap[4]) + 0.5 * trap[2],
                                               tpoints[i + 1],
                                               trap[0],
                                               nstart))
                triIndex.extend(_triangularize(abs(trap[3]) - 0.5 * trap[1], abs(trap[3]) + 0.5 * trap[1], tpoints[i],
                                               abs(trap[4]) - 0.5 * trap[2], abs(trap[4]) + 0.5 * trap[2],
                                               tpoints[i + 1],
                                               trap[0],
                                               nstart))
            else:
                triIndex.extend(_triangularize(trap[3] - 0.5 * trap[1], trap[3] + 0.5 * trap[1], tpoints[i],
                                               trap[4] - 0.5 * trap[2], trap[4] + 0.5 * trap[2], tpoints[i + 1],
                                               trap[0],
                                               nstart))

            nstart += <int> len(tpoints[i])
        # genero los puntos
        pointList = []
        for ttrap, h in zip(tpoints, hpoints):
            for py in ttrap:
                pointList.append([0, py, h])
        surfacesIndex = _getSurfaceByTriangleIndex(triIndex)
        surfaces = []
        for tsurf in surfacesIndex:
            pdist = set(tsurf)
            newSurf = Surface()
            newSurf.vectMat._m = <double *> PyMem_Malloc(3 * len(pdist) * sizeof(double))
            newSurf.vectMat._rows = <unsigned int> len(pdist)
            newSurf.vectMat._cols = 3
            replDict = {}
            i = 0
            for tindex in tsurf:
                if replDict.get(tindex) is None:
                    for j in range(3):
                        newSurf.vectMat._m[i * 3 + j] = pointList[tindex][j]
                    replDict[tindex] = i
                    i += 1
            newSurf.indTri.setList([replDict[i] for i in tsurf])
            sortTrianglesByNormal(newSurf.vectMat._m, newSurf.indTri._m, newSurf.indTri._ntriangle, normal)
            tPerimIndex = _getEdgeByTriangleIndex([replDict[i] for i in tsurf])
            newSurf.indPer.setList(tPerimIndex)
            # perimIndex = []
            # tPerimIndex = _getPerimeterByTriangleIndex(tsurf)
            # # ordeno los perimetros de izquierda a derecha. El que tenga el valor mas a la izquierda sera necesariamente
            # tPerimIndex = list(zip(tPerimIndex, [min([pointList[y][1] for y in tPerim]) for tPerim in tPerimIndex]))
            # tPerimIndex = [x[0] for x in sorted(tPerimIndex, key=_sortkey)]
            # i = 0
            # while tPerimIndex:
            #     tPerim = [replDict[j] for j in tPerimIndex.pop(0)]
            #     indPer.setList(tPerim)
            #     AlgWire.getNormal(vtemp1, newSurf.vectMat._m, indPer._m, indPer._npoint)
            #     if i == 0:
            #         if not AlgVector.isEqual(vtemp1, normal):
            #             tPerim = list(reversed(tPerim))
            #         i += 1
            #     else:
            #         if not AlgVector.isEqual(vtemp1, subnormal):
            #             tPerim = list(reversed(tPerim))
            #     perimIndex.extend(tPerim)
            #     if tPerimIndex:
            #         perimIndex.append(-1)
            # newSurf.indPer.setList(perimIndex)
            surfaces.append(newSurf)
        return surfaces
    except Exception as err:
        print(err)
        raise
    finally:
        PyMem_Free(vtemp1)
        PyMem_Free(normal)
        PyMem_Free(subnormal)

cpdef tuple mainAxisInertia(Matrix tensor):
    """Returns Inertia in main axis.
    Parameters:
    -----------
    args: tensor
        numpy.array(3x3) tensor of inertia

    Output:
    -----------
    out1: diagonal matrix with inertia values
    out2: main axis in rows
    """
    eig, rot = tensor.eig(True)
    return eig, rot.transpose()

cdef class EdgeIterator():
    def __init__(self, vectMat: AlgMatrix.Matrix, indPer: AlgWire.IndexPerimeter):
        self.vectMat = vectMat
        self.indPer = indPer

    def __iter__(self):
        self.niterator = 0
        return self

    def __next__(self):
        cdef AlgWire.Wire w
        cdef list li = []
        while self.niterator < self.indPer._npoint:
            if self.indPer._m[self.niterator] == -1:
                w = AlgWire.Wire()
                w.vectMat = self.vectMat
                w.indPer.setList(li)
                self.niterator += 1
                return w
            else:
                li.append(self.indPer._m[self.niterator])
                self.niterator += 1
        if li:
            w = AlgWire.Wire()
            w.vectMat = self.vectMat
            w.indPer.setList(li)
            return w
        raise StopIteration

    def __getitem__(self, index: int):
        cdef unsigned int i
        cdef AlgWire.Wire w
        cdef list li = []
        cdef unsigned int curInd = 0
        cdef int count = 0
        while count != index and curInd < self.indPer._npoint:
            if self.indPer._m[curInd] == -1:
                count += 1
            curInd += 1
        while self.indPer._m[curInd] != -1 and curInd < self.indPer._npoint:
            li.append(self.indPer._m[curInd])
            curInd += 1
        if li:
            w = AlgWire.Wire()
            w.vectMat = self.vectMat
            w.indPer.setList(li)
            return w

    def __len__(self):
        cdef unsigned int i, curInd
        if self.indPer._npoint > 1:
            i = 1
        else:
            return 0
        for curInd in range(self.indPer._npoint):
            if self.indPer._m[curInd] == -1 and curInd < <unsigned int> (self.indPer._npoint - 1):
                i += 1
        return i

cdef class IndexTriangle():
    "Objeto simple para los indices del triangulo de una superficie"
    def __cinit__(self, data=None):
        if data:
            self.setList(data)

    def __dealloc__(self):
        PyMem_Free(self._m)

    cpdef void setList(IndexTriangle self, list data):
        cdef int i
        if self._m:
            mem = <int *> PyMem_Realloc(self._m, len(data) * sizeof(int))
        else:
            mem = <int *> PyMem_Malloc(len(data) * sizeof(int))
        if not mem:
            raise MemoryError()
        self._m = mem
        self._ntriangle = <unsigned int> (len(data) // 3)
        for i in range(3 * self._ntriangle):
            self._m[i] = <int> (data[i])

    def copy(self):
        return IndexTriangle(self.__serialize__())

    def __serialize__(self):
        cdef int i
        return [self._m[i] for i in range(self._ntriangle * 3)]

    def __len__(self):
        return self._ntriangle

    def __json__(self):
        return {'__jsoncls__': 'CLinAlg.AlgSurface:MatrixTriangle.from_JSON', 'data': self.__serialize__()}

    @classmethod
    def from_JSON(cls, jsondict):
        obj = cls(jsondict['data'])
        return obj

cdef class Surface():
    def __cinit__(self, *args, **kwargs):
        cdef unsigned int i, j
        self.vectMat = Matrix()
        self.indTri = IndexTriangle()
        self.indPer = IndexPerimeter()  # contiene -1 para separar los distintos edges.

    def __init__(self, pointList=[], **kwargs):
        cdef int numpoints = 0
        cdef int i, j
        cdef AlgWire.Wire twire
        cdef AlgVector.Vector normal
        cdef list indPer

        if pointList:
            if isinstance(pointList, (list, tuple)):
                numpoints = <int> len(pointList)
                if pointList[0] == pointList[-1]:
                    numpoints -= 1  #
                self.vectMat._m = <double *> PyMem_Malloc(numpoints * 3 * sizeof(double))
                self.vectMat._rows = <unsigned int> numpoints
                self.vectMat._cols = 3
                for i in range(numpoints):
                    for j in range(3):
                        self.vectMat._m[i * 3 + j] = pointList[i][j]
            elif isinstance(pointList, Matrix):

                numpoints = (<Matrix> pointList)._rows
                if AlgVector.isEqual(&(<Matrix> pointList)._m[0], &(<Matrix> pointList)._m[(numpoints - 1) * 3]):
                    numpoints -= 1
                if (<Matrix> pointList)._cols != 3:
                    raise TypeError("The matrix input for surface shold have 3 columns")
                self.vectMat._m = <double *> PyMem_Malloc(numpoints * 3 * sizeof(double))
                self.vectMat._rows = <unsigned int> numpoints
                self.vectMat._cols = 3
                for i in range(numpoints * 3):
                    self.vectMat._m[i] = (<Matrix> pointList)._m[i]

            if kwargs.get('indPer') is None:
                indPer = []
                for i in range(1, <int> self.vectMat._rows):
                    indPer.extend([i - 1, i, -1])
                indPer.extend([self.vectMat._rows - 1, 0, -1])
            else:
                indPer = kwargs['indPer']
            self.indPer.setList(indPer)
            if kwargs.get('indTri') is None:
                twire = AlgWire.Wire()
                twire.vectMat = self.vectMat
                tperim = list(range(<int>self.vectMat._rows))
                tperim.append(0)
                twire.indPer.setList(tperim)
                normal = twire.normal()
                tlist = tessellate(self.vectMat._m, twire.indPer._m, twire.indPer._npoint, normal._v)
                self.indTri.setList(tlist)
                sortTrianglesByNormal(self.vectMat._m, self.indTri._m, self.indTri._ntriangle, normal._v)
            else:
                self.indTri.setList(kwargs['indTri'])

    def __dealloc__(self):
        pass  #no hago nada, ya se encaragarn los objetos constituyentes de borrarse

    def __json__(self):
        return {'__jsoncls__': 'CLinAlg:Surface.from_JSON', 'mat': self.vectMat, 'indTri': self.indTri,
                'indPer': self.indPer}

    @classmethod
    def from_JSON(cls, jsondict):
        obj = cls()
        obj.vectMat = jsondict['mat']
        obj.indTri = jsondict['indTri']
        obj.indPer = jsondict['indPer']
        return obj

    @property
    def Wire(self):
        """
        Return wire of contour as unified Wire
        """
        cdef int i
        edges = []
        indices = self.indPer.__serialize__()
        if not indices[-1] == -1:
            indices.append(-1)
        cortes = [idx for idx, val in enumerate(indices) if val == -1]
        cortes.insert(0, -1)
        edges = [indices[cortes[i - 1] + 1: cortes[i]] for i in range(1, len(cortes))]
        wires = []
        founded = False
        curWire = edges.pop(0)
        i = 0
        while edges:
            if edges[i][0] == curWire[-1]:
                curWire.extend(edges[i][1:])
                edges.pop(i)
                i = 0
                continue
            elif edges[i][-1] == curWire[-1]:
                curWire.extend(list(reversed([i][:-1])))
                edges.pop(i)
                i = 0
                continue
            i += 1
            if i == len(edges):
                wires.extend(curWire)
                wires.append(-1)
                curWire = edges.pop(0)
                i = 0
        wires.extend(curWire)
        wires.append(-1)
        indPer = AlgWire.IndexPerimeter()
        indPer.setList(wires)
        return EdgeIterator(self.vectMat, indPer)

    @property
    def Edge(self):
        """
        Returns edge iterator
        """
        return EdgeIterator(self.vectMat, self.indPer)

    @property
    def Point(self):
        "Return points as iterator"
        cdef AlgWire.IndexPerimeter indper = AlgWire.IndexPerimeter()
        indper.setList(list(range(self.vectMat._rows)))
        return AlgWire.PointIterator(self.vectMat, indper)

    def area(self):
        return getArea(self)

    def cdg(self):
        cdef AlgVector.Vector cdg = AlgVector.Vector()
        getCDG(cdg._v, self)
        return cdg

    def perimeter(self):
        return sum([x.perimeter() for x in self.Wire])

    def inertia(self):
        return getInertia(self)

    def staticMoment(self, line: AlgLine.Line):
        "Returns static moment to a line (point and direction)"
        return staticMoment(self, line._pnt, line._dir)

    def modulus(self):
        "Return section modulus at main axis of inertia"
        cdef list inertia
        cdef Matrix axis
        inertia, axis = mainAxisInertia(self.inertia())  #lista de autovalores y matriz de autovectores
        cdg = <AlgVector.Vector> self.cdg()
        modulus = []
        for i in range(3):
            maxDist = 0
            for j in range(self.vectMat._rows - 1):
                maxDist = max(maxDist, AlgLine.distanceToPoint(cdg._v, &axis._m[i * 3], &self.vectMat._m[j * 3]))
            modulus.append(inertia[i] / maxDist)
        return modulus, axis

    def radius(self):
        "Return radius of section in main axis of inertia"
        cdef list inertia
        cdef Matrix axis
        inertia, axis = mainAxisInertia(self.inertia())
        area = self.area()
        _radiusYZ = []
        for i in range(3):
            _radiusYZ.append((inertia[i] / area) ** 0.5)
        return _radiusYZ, axis

    def isInside(self, point: AlgVector.Vector, incEdge=True):
        "Funcion que calcula si un punto esta dentro de la superficie"
        return isInside(self, point._v, <bint> incEdge)

    def getRotate(self, center: AlgVector.Vector, axis: AlgVector.Vector, angle: float):
        cdef double * quat = <double *> PyMem_Malloc(4 * sizeof(double))
        cdef Surface newSurface = Surface()
        newSurface.vectMat._m = <double *> PyMem_Malloc(3 * self.vectMat._rows * sizeof(double))
        newSurface.vectMat._rows = self.vectMat._rows
        newSurface.vectMat._cols = 3
        AlgQuaternion.quatFromAxisAngle(quat, axis._v, angle)
        AlgQuaternion.rotateVectorMatrix(newSurface.vectMat._m, quat, center._v, self.vectMat._m, self.vectMat._rows)
        newSurface.indTri.setList(self.indTri.__serialize__())
        newSurface.indPer.setList(self.indPer.__serialize__())
        PyMem_Free(quat)
        return newSurface

    def getTranslate(self, translateVector: AlgVector.Vector):
        """
        Returns a copy of surface translate by the given vector
        :param translateVector: Vector with the translation
        """
        cdef Surface newSurface = Surface()
        cdef int i, j
        newSurface.vectMat._m = <double *> PyMem_Malloc(self.vectMat._rows * 3 * sizeof(double))
        newSurface.vectMat._rows = self.vectMat._rows
        newSurface.vectMat._cols = self.vectMat._cols
        for i in range(<int> self.vectMat._rows):
            for j in range(3):
                newSurface.vectMat._m[j + i * 3] = self.vectMat._m[j + i * 3] + translateVector._v[j]
        newSurface.indTri.setList(self.indTri.__serialize__())
        newSurface.indPer.setList(self.indPer.__serialize__())
        return newSurface

    def splitByLine(self, line: AlgLine.Line):
        return splitByLine(self, line._pnt, line._dir)

    def intersectByLine(self, line: AlgLine.Line, incEdge=True):
        "Returns a list of intersection points with a line"
        # puntos = list(range(self.vectMat._rows))
        cortes = AlgWire.getIntersectionPointWithLine(self.vectMat._m, self.indPer._m, self.indPer._npoint, line._pnt,
                                                      line._dir,
                                                      <bint> incEdge)
        if len(cortes) < 2:
            return []
        return cortes

    def intersects(self, other: Surface) -> bool:
        "Returns True if surfaces intersect with each other"
        p: AlgVector.Vector
        for p in other.Points:
            if isInside(self, p._v, <bint> True):
                return True
        for p in self.Points:
            if isInside(other, p._v, <bint> True):
                return True
        intPoints = AlgWire.getIntersectionPointWithWire(self.vectMat._m, self.indPer._m, self.indPer._npoint,
                                                         other.vectMat._m, other.indPer._m, other.indPer._npoint,
                                                         <bint> True)
        return intPoints != []

    def boundRectangle(self, b0, b1) -> tuple:
        "Returns the bounding rectangle of the surface"
        if isinstance(b0, str):
            b0 = ['x', 'y', 'z'].index(b0)
            b1 = ['x', 'y', 'z'].index(b1)
        return boundRectangle(self, b0, b1)

    def copy(self):
        "Returns a copy of the surface"
        newSurf = Surface()
        newSurf.vectMat = self.vectMat.copy()
        newSurf.indTri.setList(self.indTri.__serialize__())
        newSurf.indPer.setList(self.indPer.__serialize__())
        return newSurf

    cpdef Surface transform(Surface self, Transformation transf):
        "return surface with transformation"
        cdef unsigned int v, i, j
        cdef double val
        cdef Surface newSurface = Surface()
        newSurface.vectMat._m = <double *> PyMem_Malloc(self.vectMat._rows * self.vectMat._cols * sizeof(double))
        newSurface.vectMat._rows = self.vectMat._rows
        newSurface.vectMat._cols = self.vectMat._cols
        for v in range(self.vectMat._rows):
            for i in range(3):
                val = 0
                for j in range(3):
                    val += transf._mat._m[j + 4 * i] * self.vectMat._m[j + v * 3]
                val += transf._mat._m[3 + 4 * i]
                newSurface.vectMat._m[i + v * 3] = val
        newSurface.indPer.setList(self.indPer.__serialize__())
        newSurface.indTri.setList(self.indTri.__serialize__())
        return newSurface

    def normal(Surface self):
        cdef double * v1 = <double *> PyMem_Malloc(3 * sizeof(double))
        cdef double * v2 = <double *> PyMem_Malloc(3 * sizeof(double))
        cdef AlgVector.Vector normal = AlgVector.Vector()
        AlgVector.sub(v1, &self.vectMat._m[self.indTri._m[1] * 3], &self.vectMat._m[self.indTri._m[0] * 3])
        AlgVector.sub(v2, &self.vectMat._m[self.indTri._m[2] * 3], &self.vectMat._m[self.indTri._m[0] * 3])
        AlgVector.cross(normal._v, v1, v2)
        normal.normalize()
        return normal

    def invertOrientation(Surface self):
        "Function that invert current orientation for surface"
        cdef AlgVector.Vector normal = self.normal()
        cdef int i
        for i in range(3):
            normal._v[i] = -normal._v[i]
        sortTrianglesByNormal(self.vectMat._m, self.indTri._m, self.indTri._ntriangle, normal._v)

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memodict=None):
        return self.copy()

    def __repr__(self):
        cdg: AlgVector.Vector = self.cdg()
        return "Surface(Area=%.2f; Perimeter=%.2f; CDG=(%.2f,%.2f,%.2f))" % (
            self.area(), self.perimeter(), cdg.x, cdg.y, cdg.z)

    def __repr__(self):
        cdg: AlgVector.Vector = self.cdg()
        return "Surface(Area=%.2f; Perimeter=%.2f; CDG=(%.2f,%.2f,%.2f))" % (
            self.area(), self.perimeter(), cdg.x, cdg.y, cdg.z)

    def indexTriangles(self):
        return self.indTri.__serialize__()

    def indexPerimeter(self):
        return self.indPer.__serialize__()
