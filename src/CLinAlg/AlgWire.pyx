import sys
from cpython.mem cimport PyMem_Malloc, PyMem_Free, PyMem_Realloc
from libc.math cimport M_PI, sin, fabs
from .AlgVector cimport Vector
from .AlgTransformation cimport Transformation
from . cimport AlgVector, AlgSegment, AlgQuaternion, AlgMatrix, AlgLine
from .AlgTool cimport presition


cdef double getLength(double * wirepoint, int * indexes, unsigned int numvect):
    cdef double length = 0
    cdef unsigned int i
    for i in range(<unsigned int>(numvect - 1)):
        length += AlgVector.distance(&wirepoint[indexes[i] * 3], &wirepoint[indexes[i + 1] * 3])
    return length


cdef bint isInside(double * wirepoint, int * indexes, unsigned int numvect, double * punto, bint incEdges):
    cdef unsigned int i
    for i in range(<unsigned int> (numvect - 1)):
        if AlgSegment.isInside(&wirepoint[indexes[i] * 3], &wirepoint[indexes[i + 1] * 3], punto, incEdges):
            return True
    return False

cdef bint isClosed(double * wirePoint, int * indexes, int numpoints):
    if indexes[0] == indexes[numpoints-1]:
        return True
    else:
        return all([fabs(wirePoint[indexes[0]*3+i] - wirePoint[indexes[numpoints - 1] * 3 + i]) < presition for i in range(3)])


cdef void getRotate(double * toPointer, double * wirePoint, unsigned int numpoints, double * center, double * axis,
                            double angle):
    cdef double * quat = <double *> PyMem_Malloc (4 * sizeof(double))
    AlgQuaternion.quatFromAxisAngle(quat, axis, angle)
    AlgQuaternion.rotateVectorMatrix(toPointer, quat, center, wirePoint, numpoints)
    PyMem_Free(quat)

cdef void getCDG(double * cdg, double * wirePoint, int * indexes, unsigned int numpoints):
    cdef unsigned int i, j
    cdef unsigned int numrows = numpoints
    if isClosed(wirePoint, indexes, numpoints):
        numrows -= 1
    for i in range(3):
        cdg[i] = 0
        for j in range(numrows):
            cdg[i] += wirePoint[i + indexes[j] * 3]
        cdg[i] = cdg[i] / numrows

cdef bint getPointFromOrigin(double * toPoint, double * wirePoint, int * indexes, unsigned int numpoints, double length):
    cdef unsigned int i, j
    cdef double tlen = 0
    cdef double prop, slen
    if length <= 0:
        for i in range(3):
            toPoint[i] = wirePoint[indexes[0]+i]
        return True
    elif length > getLength(wirePoint, indexes, numpoints):
        for i in range(3):
            toPoint[i] = wirePoint[indexes[numpoints - 1] * 3 + i]
        return True
    for i in range(<unsigned int>(numpoints-1)):
        slen = AlgVector.distance(&wirePoint[indexes[i]*3], &wirePoint[indexes[i+1]*3])
        if slen and tlen + slen >= length:
            prop = (length - tlen) / slen
            for j in range(3):
                toPoint[j] = wirePoint[indexes[i]*3+j] + prop * (wirePoint[indexes[i+1]*3+j] - wirePoint[indexes[i]*3+j])
            return True
        tlen += slen
    return False

cdef double getLengthFromOrigin(double * wirePoint, int * indexes, unsigned int numpoints, double * point):
    cdef double dist = 0
    cdef double ldist
    cdef unsigned int i
    for i in range(<unsigned int> (numpoints - 1)):
        if AlgSegment.isInside(&wirePoint[indexes[i] * 3],
                                    &wirePoint[indexes[i + 1] * 3],
                                    point, <bint>True):
            dist += AlgVector.distance(&wirePoint[indexes[i] * 3], point)
            return dist
        else:
            dist += AlgVector.distance(&wirePoint[indexes[i] * 3],
                                 &wirePoint[indexes[i + 1] * 3])
    return sys.float_info.max

cdef bint getDirectionAtPoint(double * toDir, double * wirePoint, int * indexes, unsigned int numpoints, double * point):
    cdef int i
    for i in range(numpoints-1):
        if AlgVector.isEqual(&wirePoint[indexes[i]*3], point) or AlgSegment.isInside(&wirePoint[indexes[i]*3], &wirePoint[indexes[i+1]*3], point, <bint>False):
            AlgVector.vdir(toDir, &wirePoint[indexes[i]*3], &wirePoint[indexes[i+1]*3])
            return True
    return False

cdef void getNearestPoint(double * toPoint, double * wirePoint, int * indexes, unsigned int numpoints, double * point):
    cdef unsigned int i, j
    cdef double val, curval
    cdef double * distpoint = <double *> PyMem_Malloc(3 * sizeof(double))
    val = sys.float_info.max
    for i in range(<unsigned int> (numpoints - 1)):
        AlgSegment.getNearestPoint(distpoint, &wirePoint[indexes[i] * 3],
                                                    &wirePoint[indexes[i + 1] * 3],
                                                    point)
        curval = AlgVector.distance(distpoint, point)
        if curval < val:
            val = curval
            AlgVector.copy(toPoint, distpoint)
    PyMem_Free(distpoint)

cdef bint isClockWise(double * wirePoint, int * indexes, unsigned int numpoints, double * obsPoint):
    "Determine if points are clockwise in reference to an observer point"
    cdef double * cdg = <double *> PyMem_Malloc (3 * sizeof(double))
    cdef double * v = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * vt1 = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * vt2 = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * vt3 = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef unsigned int i, j
    cdef double modulo, angle
    getCDG(cdg, wirePoint, indexes,  numpoints)
    try:
        for j in range(3):
            v[j] = 0
        for i in range(<unsigned int> (numpoints - 1)):
            AlgVector.sub(vt1, &wirePoint[indexes[i + 1] * 3], &wirePoint[indexes[i] * 3])
            AlgVector.sub(vt2, &wirePoint[indexes[i] * 3], cdg)
            AlgVector.cross(vt3, vt1, vt2)
            for j in range(3):
                v[j] = v[j] + vt3[j]
        modulo = AlgVector.module(v)
        for j in range(3):
            v[j] = v[j] / modulo
        AlgVector.sub(vt1, obsPoint, cdg)
        modulo = AlgVector.module(vt1)
        for j in range(3):
            vt1[j] = vt1[j] / modulo
        angle = AlgVector.angle(vt1, v)
        if angle < M_PI / 2:
            return True
        elif angle > M_PI / 2:
            return False
    finally:
        PyMem_Free(cdg)
        PyMem_Free(v)
        PyMem_Free(vt1)
        PyMem_Free(vt2)
        PyMem_Free(vt3)


cdef list getShifted(double * wirePoint, int * indexes, unsigned int numpoints, double value):
    """
    Retrun list of vectors that has the shifted points for a wire
    :param wirePoint: double pointer to vectors
    :param indexes: int pointer for sequential index for wire
    :param numpoints: number of points of wire
    :param value: shift value, positiv expand, negativ is contract
    :return: 
    """
    cdef unsigned int i, nump, iprev, inext, j
    cdef double * vt1 = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * vt2 = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * vt3 = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * vt4 = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * obs = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * cdg = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double ang, shp
    getCDG(cdg, wirePoint, indexes, numpoints)
    try:
        vlist = []
        if numpoints == 2:
            for i in range(2):
                vlist.append(Vector(wirePoint[indexes[i]* 3], wirePoint[indexes[i ]* 3 + 1], wirePoint[indexes[i] * 3 + 2]))
            return vlist
        for j in range(3):
            obs[j] = 0
        for i in range(<unsigned int> (numpoints - 1)):
            AlgVector.sub(vt1, &wirePoint[indexes[i+1]*3], &wirePoint[indexes[i]*3])
            AlgVector.sub(vt2, &wirePoint[indexes[i]*3], cdg)
            AlgVector.cross(vt3, vt1, vt2)
            AlgVector.add(obs, vt3, obs)
        shp = AlgVector.module(obs)
        for j in range(3):
            obs[j] = obs[j] / shp
        if isClosed(wirePoint, indexes, numpoints):
            nump = numpoints - 1
            iprev = <unsigned int> ((-1) % nump)
            for i in range(nump):
                inext = <unsigned int> ((i + 1) % nump)
                AlgVector.vdir(vt1, &wirePoint[indexes[iprev]*3], &wirePoint[indexes[i]*3])
                AlgVector.vdir(vt2, &wirePoint[indexes[i]*3], &wirePoint[indexes[inext]*3])
                AlgVector.cross(vt4, vt2, vt1)
                if AlgVector.module(vt4) == 0:
                    continue
                AlgVector.sub(vt3, vt2, vt1)
                ang = AlgVector.angle(vt3, vt2)
                shp = sin(ang) * AlgVector.module(vt3)
                for j in range(3):
                    vt3[j] = vt3[j] * value / shp

                ang = AlgVector.angle(vt4, obs)
                if ang <= M_PI / 2:
                    vlist.append(Vector(wirePoint[indexes[i] * 3] + vt3[0], wirePoint[indexes[i] * 3 + 1] + vt3[1],
                                        wirePoint[indexes[i]* 3 + 2] + vt3[2]))
                else:
                    vlist.append(Vector(wirePoint[indexes[i] * 3] - vt3[0], wirePoint[indexes[i] * 3 + 1] - vt3[1],
                                        wirePoint[indexes[i] * 3 + 2] - vt3[2]))
                iprev = i
            vlist.append(vlist[0].copy())
        else:
            nump = numpoints
            iprev = 0
            for i in range(1, <unsigned int> (nump - 1)):
                inext = <unsigned int> ((i + 1) % nump)
                AlgVector.vdir(vt1, &wirePoint[3*indexes[iprev]], &wirePoint[3*indexes[i]])
                AlgVector.vdir(vt2, &wirePoint[3*indexes[i]], &wirePoint[3*indexes[inext]])
                if AlgVector.isEqual(vt1, vt2):
                    continue
                AlgVector.sub(vt3, vt2, vt1)
                ang = AlgVector.angle(vt3, vt2)
                shp = sin(ang) * AlgVector.module(vt3)
                for j in range(3):
                    vt3[j] = vt3[j] * value / shp
                AlgVector.cross(vt4, vt2, vt1)
                ang = AlgVector.angle(vt4, obs)
                if len(vlist) == 0:
                    AlgVector.cross(vt2, obs, vt1)
                    shp = AlgVector.module(vt2)
                    for j in range(3):
                        vt2[j] = vt2[j] * value / shp
                    if ang > M_PI / 2:
                        vlist.append(Vector(wirePoint[indexes[0]*3] + vt2[0], wirePoint[indexes[0]*3+1] + vt2[1],
                                            wirePoint[indexes[0]*3+2] + vt2[2]))
                    else:
                        vlist.append(Vector(wirePoint[indexes[0]*3] - vt2[0], wirePoint[indexes[0]*3+1] - vt2[1],
                                            wirePoint[indexes[0]*3+2] - vt2[2]))
                if ang <= M_PI / 2:
                    vlist.append(Vector(wirePoint[indexes[i]*3] + vt3[0], wirePoint[indexes[i]*3 + 1] + vt3[1],
                                        wirePoint[indexes[i]*3 + 2] + vt3[2]))
                else:
                    vlist.append(Vector(wirePoint[indexes[i]*3] - vt3[0], wirePoint[indexes[i]*3 + 1] - vt3[1],
                                        wirePoint[indexes[i]*3 + 2] - vt3[2]))
                iprev = i
            #agrego el ultimo
            AlgVector.cross(vt1, obs, vt2)
            shp = AlgVector.module(vt1)
            for j in range(3):
                vt1[j] = vt1[j] * value / shp
            iprev = <unsigned int> (numpoints- 1)
            if ang > M_PI / 2:
                vlist.append(Vector(wirePoint[indexes[iprev] * 3] + vt1[0], wirePoint[indexes[iprev] * 3 + 1] + vt1[1],
                                    wirePoint[indexes[iprev] * 3 + 2] + vt1[2]))
            else:
                vlist.append(Vector(wirePoint[indexes[iprev] * 3] - vt1[0], wirePoint[indexes[iprev] * 3 + 1] - vt1[1],
                                    wirePoint[indexes[iprev] * 3 + 2] - vt1[2]))
        return vlist
    finally:
        PyMem_Free(vt1)
        PyMem_Free(vt2)
        PyMem_Free(vt3)
        PyMem_Free(vt4)
        PyMem_Free(obs)
        PyMem_Free(cdg)

cdef list getIntersectionPointWithWire(double * wirePoint, int * indexes, unsigned int numpoints,
                                       double * otherWirePoint,  int * otherIndexes, unsigned int otherNumpoints,
                                       bint incEdge):
    "Return the intersection point between this wire and another wire"
    cdef int i, j
    cdef double * ipoint = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * tpoint = <double *> PyMem_Malloc (3 * sizeof(double))
    cdef bint sclosed = isClosed(wirePoint, indexes, numpoints)
    cdef bint oclosed = isClosed(otherWirePoint, otherIndexes, otherNumpoints)
    cdef double * dirMIn = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * dirMOut = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * dirOIn = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * dirOOut = <double *> PyMem_Malloc(3 * sizeof(double))
    try:
        cortes = []
        for i in range(<int>numpoints - 1):
            for j in range(<int>otherNumpoints - 1):
                if not AlgSegment.getIntersectionPointWithSegment(ipoint, &wirePoint[indexes[i] * 3],
                                                                         &wirePoint[indexes[i + 1] * 3],
                                                                         &otherWirePoint[otherIndexes[j] * 3],
                                                                         &otherWirePoint[otherIndexes[j + 1] * 3],
                                                                         <bint>True) or any([AlgVector.isEqual(ipoint, (<Vector> v)._v) for v in cortes]):

                    continue
                if (incEdge or (not AlgVector.isEqual(ipoint, &wirePoint[indexes[i] * 3])
                                and not AlgVector.isEqual(ipoint, &wirePoint[indexes[i + 1] * 3])
                                and not AlgVector.isEqual(ipoint, &otherWirePoint[otherIndexes[j] * 3])
                                and not AlgVector.isEqual(ipoint, &otherWirePoint[otherIndexes[j + 1] * 3]))):
                    newVect = Vector(None)
                    newVect._v = ipoint
                    cortes.append(newVect)
                    ipoint = <double *> PyMem_Malloc(3 * sizeof(double))
                    continue

                #vectors iniciales
                if AlgVector.isEqual(ipoint, &wirePoint[indexes[i] * 3]):  # punto inicial
                    if i == 0:
                        if not sclosed:
                            continue
                        AlgVector.vdir(dirMIn, &wirePoint[indexes[numpoints - 1] * 3],
                              &wirePoint[indexes[numpoints - 2] * 3])
                    else:
                        AlgVector.vdir(dirMIn, &wirePoint[i * 3], &wirePoint[(i - 1) * 3])
                    AlgVector.vdir(dirMOut, &wirePoint[i * 3], &wirePoint[(i + 1) * 3])
                elif AlgVector.isEqual(ipoint, &wirePoint[(i + 1) * 3]):  #punto final
                    if i == (numpoints - 2):
                        if not sclosed:
                            continue
                        AlgVector.vdir(dirMOut, &wirePoint[0], &wirePoint[3])
                    else:
                        AlgVector.vdir(dirMOut, &wirePoint[i * 3], &wirePoint[(i + 1) * 3])
                    AlgVector.vdir(dirMIn, &wirePoint[i * 3], &wirePoint[(i - 1) * 3])
                else:  #ninguna de las anteriores
                    AlgVector.vdir(dirMOut, &wirePoint[i * 3], &wirePoint[(i + 1) * 3])
                    dirMIn[0] = -dirMOut[0]
                    dirMIn[1] = -dirMOut[1]
                    dirMIn[2] = -dirMOut[2]
                #vectores finales
                if AlgVector.isEqual(ipoint, &otherWirePoint[j * 3]):  # punto inicial
                    if j == 0:
                        if not oclosed:
                            continue
                        AlgVector.vdir(dirOIn, &otherWirePoint[(otherNumpoints - 1) * 3],
                              &otherWirePoint[(otherNumpoints - 2) * 3])
                    else:
                        AlgVector.vdir(dirOIn, &otherWirePoint[j * 3], &otherWirePoint[(j - 1) * 3])
                    AlgVector.vdir(dirOOut, &otherWirePoint[j * 3], &otherWirePoint[(j + 1) * 3])
                elif AlgVector.isEqual(ipoint, &otherWirePoint[(j + 1) * 3]):  #punto final
                    if j == (otherNumpoints - 2):
                        if not oclosed:
                            continue
                        AlgVector.vdir(dirOOut, &otherWirePoint[0], &otherWirePoint[3])
                    else:
                        AlgVector.vdir(dirOOut, &otherWirePoint[j * 3], &otherWirePoint[(j + 1) * 3])
                    AlgVector.vdir(dirOIn, &otherWirePoint[j * 3], &otherWirePoint[(j - 1) * 3])
                else:  #ninguna de las anteriores
                    AlgVector.vdir(dirOOut, &otherWirePoint[j * 3], &otherWirePoint[(j + 1) * 3])
                    dirOIn[0] = -dirOOut[0]
                    dirOIn[1] = -dirOOut[1]
                    dirOIn[2] = -dirOOut[2]

                if AlgSegment.getIntersectionPointWithSegment(tpoint, dirMIn, dirMOut, dirOIn, dirOOut, <bint>False):
                    newVect = Vector(None)
                    newVect._v = ipoint
                    cortes.append(newVect)
                    ipoint = <double *> PyMem_Malloc(3 * sizeof(double))
        return cortes
    finally:
        PyMem_Free(dirMIn)
        PyMem_Free(dirMOut)
        PyMem_Free(dirOIn)
        PyMem_Free(dirOOut)
        PyMem_Free(tpoint)
        PyMem_Free(ipoint)


cdef list getIntersectionPointWithLine(double * wirePoint, int * indexes, unsigned int numpoints, double * linePoint, double * lineDir,
                                           bint incEdge):
    # Funcion que genera el mismo algoritmo que BaseWire.c_getINtersecitonPointWidthLine pero usando un secuenciador de puntos
    cdef int i, pp, pc, pn
    cdef double * ipoint = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * tpoint = <double *> PyMem_Malloc (3 * sizeof(double))
    cdef bint sclosed = isClosed(wirePoint, indexes, numpoints)
    cdef double * dirMIn = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * dirMOut = <double *> PyMem_Malloc(3 * sizeof(double))

    try:
        cortes = []
        for i in range(numpoints - 1):
            pc = indexes[i]
            pn = indexes[i+1]
            if not AlgSegment.getIntersectionPointWithLine(ipoint, &wirePoint[pc * 3],
                                                                  &wirePoint[pn * 3],
                                                                  linePoint,
                                                                  lineDir,
                                                                  incEdge) or any([AlgVector.isEqual(ipoint, (<Vector> v)._v) for v in cortes]):

                continue
            if (incEdge or (not AlgVector.isEqual(ipoint, &wirePoint[pc * 3])
                            and not AlgVector.isEqual(ipoint, &wirePoint[pn * 3]))):
                newVect = Vector(None)
                newVect._v = ipoint
                cortes.append(newVect)
                ipoint = <double *> PyMem_Malloc(3 * sizeof(double))
                continue

            #vectors iniciales. Condicion a la que salto porque me ha coincidido sobre un vector
            if AlgVector.isEqual(ipoint, &wirePoint[pc * 3]):  # punto inicial
                if i == 0:
                    if not sclosed:
                        continue
                    pp = indexes[numpoints-2]
                    AlgVector.vdir(dirMIn, &wirePoint[pc * 3], &wirePoint[pp * 3])
                else:
                    pp = indexes[numpoints-1]
                    AlgVector.vdir(dirMIn, &wirePoint[pc * 3], &wirePoint[pp * 3])
                AlgVector.vdir(dirMOut, &wirePoint[pc * 3], &wirePoint[pn * 3])
                AlgVector.add(dirMIn, &wirePoint[pc * 3], dirMIn)
                AlgVector.add(dirMOut, &wirePoint[pc * 3], dirMOut)
            elif AlgVector.isEqual(ipoint, &wirePoint[pn * 3]):  #punto final
                if i == numpoints - 2:
                    if not sclosed:
                        continue
                    pp = indexes[1]
                    AlgVector.vdir(dirMOut, &wirePoint[pn * 3], &wirePoint[pp * 3])
                else:
                    pp = indexes[i+2]
                    AlgVector.vdir(dirMOut, &wirePoint[pn * 3], &wirePoint[pp * 3])
                AlgVector.vdir(dirMIn, &wirePoint[pn * 3], &wirePoint[pc * 3])
                AlgVector.add(dirMIn, &wirePoint[pn * 3], dirMIn)
                AlgVector.add(dirMOut, &wirePoint[pn * 3], dirMOut)
            else:  #ninguna de las anteriores
                raise ValueError('No es posible esta combinacion')
            if AlgSegment.getIntersectionPointWithLine(tpoint, dirMIn, dirMOut, linePoint,
                                                                  lineDir, <bint>False):
                newVect = Vector(None)
                newVect._v = ipoint
                cortes.append(newVect)
                ipoint = <double *> PyMem_Malloc (3 * sizeof(double))

        return cortes
    finally:
        PyMem_Free(dirMIn)
        PyMem_Free(dirMOut)
        PyMem_Free(tpoint)
        PyMem_Free(ipoint)


cdef list getIntersectionPointWithSegment(double * wirePoint, int * indexes, unsigned int numpoints,
                                          double * vini, double * vfin,
                                          bint incEdge):
    "Return the intersection point between this wire and another wire"
    cdef double * segDir = <double *> PyMem_Malloc (3 * sizeof(double))
    try:
        AlgVector.vdir(segDir, vini, vfin)
        linecortes = getIntersectionPointWithLine(wirePoint, indexes, numpoints, vini, segDir, incEdge)
        cortes = []
        for i in range(len(linecortes)):
            if AlgSegment.isInside(vini, vfin, (<Vector>linecortes[i])._v, incEdge):
                cortes.append(linecortes[i])
        return cortes
    finally:
        PyMem_Free(segDir)


cdef Wire getSubWire(double * wirePoint, int * indexes, unsigned int numpoints, double length1, double length2):
    "Funcion que devuelve un subwire desde los puntos inicial y final. Destruye ciclos cerrados"
    cdef bint revert = False
    cdef double _l = 0
    cdef double elength, prop
    cdef int i, j
    if length2 < length1:
        revert = True
        length1, length2 = length2, length1
    length1 = min(<double>0, length1)
    length2 = max(getLength(wirePoint, indexes, numpoints), length2)
    puntos = []
    for i in range(<int>numpoints - 1):
        elength = AlgVector.distance(&wirePoint[indexes[i] * 3], &wirePoint[indexes[i + 1] * 3])
        if elength < 2 * presition:
            continue
        if _l > length2:
            break
        if _l - presition <= length1 < _l + elength - presition:  # estamos en el punto inicial
            prop = max((length1 - _l) / elength, 0)
            puntos.append(
                [wirePoint[indexes[i] * 3 + j] + prop * (wirePoint[indexes[i + 1] * 3 + j] - wirePoint[indexes[i] * 3 + j])
                 for j in range(3)])
        if length1 + presition < _l < length2 - presition:
            puntos.append([wirePoint[indexes[i] * 3 + j] for j in range(3)])
        if _l + presition < length2 <= _l + elength + presition:  # estamos en el punto final
            propfin = min(1, <int>((length2 - _l) / elength))
            puntos.append(
                [wirePoint[indexes[i] * 3 + j] + prop * (wirePoint[indexes[i + 1] * 3 + j] - wirePoint[indexes[i] * 3 + j])
                 for j in range(3)])
            break
        _l += elength
    if not puntos:
        return
    if revert:
        return Wire(list(reversed(puntos)))
    else:
        return Wire(puntos)
    

cdef void getNormal(double * normal, double * wirePoint, int * indexes, unsigned int numpoints):
    "Funcion que devuelve la normal de un wire. Siempre define la direccion contrarreloj"
    cdef double * cdg = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * vtemp1 = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * vtemp2 = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * vtemp3 = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double modulo
    cdef int i, j, numrows
    try:
        getCDG(cdg, wirePoint, indexes, numpoints)
        numrows = <int>numpoints
        if isClosed(wirePoint, indexes, numpoints):
            numrows -= 1
        for i in range(3):
            normal[i] = 0
        # eje principal
        for i in range(numrows - 1):
            AlgVector.sub(vtemp1, &wirePoint[indexes[i + 1] * 3], &wirePoint[indexes[i] * 3])
            AlgVector.sub(vtemp2, &wirePoint[indexes[i] * 3], cdg)
            AlgVector.cross(vtemp3, vtemp2, vtemp1)
            AlgVector.add(normal, normal, vtemp3)
        modulo = AlgVector.module(normal)
        if modulo:
            for i in range(3):
                normal[i] = normal[i] / modulo
    finally:
        PyMem_Free(vtemp1)
        PyMem_Free(vtemp2)
        PyMem_Free(vtemp3)
        PyMem_Free(cdg)

#..........................................CLASES .....................................................................

cdef class IndexPerimeter():
    "Class matrix that stores the index of perimeter for a group of doubles that contains vectors"

    def __cinit__(self, data=None):
        if data:
            self.setList(data)

    def __dealloc__(self):
        PyMem_Free(self._m)

    cpdef void setList(IndexPerimeter self, list cdata):
        cdef int i
        if self._m:
            mem = <int *> PyMem_Realloc(self._m, len(cdata) * sizeof(int))
            if not mem:
                raise MemoryError()
            self._m = mem
        else:
            self._m = <int *> PyMem_Malloc(len(cdata) * sizeof(int))
            if not self._m:
                raise MemoryError()
        self._npoint = <unsigned int>(len(cdata))
        for i in range(len(cdata)):
            self._m[i] = <int>(cdata[i])

    def __getitem__(self, item:int):
        return self._m[item % self._npoint]

    def __iter__(self):
        self.niterator = 0
        return self

    def __next__(self):
        if self.niterator < self._npoint:
            try:
                return self._m[self.niterator]
            finally:
                self.niterator += 1
        else:
            raise StopIteration

    def __contains__(self, item):
        cdef unsigned int i
        for i in range(self._npoint):
            if self._m[i] == item:
                return True
        return False

    def copy(self):
        return IndexPerimeter(self.__serialize__())

    def __serialize__(self):
        cdef unsigned int i
        return [self._m[i] for i in range(self._npoint)]

    def __len__(self):
        return self._npoint

    def __json__(self):
        return {'__jsoncls__': 'CLinAlg.AlgWire:IndexPerimeter.from_JSON', 'data': self.__serialize__()}

    @classmethod
    def from_JSON(cls, jsondict):
        obj = cls(jsondict['data'])
        return obj


cdef class PointIterator():

    def __init__(self, vectMat:AlgMatrix.Matrix, indPer:IndexPerimeter):
        self.vectMat = vectMat # keep the vector point
        self.indPer = indPer # keep the secuential index

    # def __dealloc__(self):
    #     pass

    def __iter__(PointIterator self):
        self.niterator = 0
        return self

    def __next__(PointIterator self):
        cdef AlgVector.PhantomVector v
        if self.niterator < self.indPer._npoint:
            v = AlgVector.PhantomVector()
            v._v = &self.vectMat._m[self.indPer[self.niterator] * 3]
            self.niterator += 1
            return v
        else:
            raise StopIteration

    def __getitem__(self, index):
        cdef unsigned int i
        cdef int start, stop, step
        cdef AlgVector.PhantomVector phv
        if isinstance(index, slice):
            start = 0 if index.start is None else index.start if index.start >= 0 else self.indPer._npoint + index.start
            stop = self.indPer._npoint if index.stop is None else index.stop if index.stop >= 0 else self.indPer._npoint + index.stop
            step = 1 if index.step is None else index.step
            phlist = []
            for i in range(start, stop, step):
                phv = AlgVector.PhantomVector()
                phv._v = &self.vectMat._m[self.indPer[i % self.indPer._npoint]*3]
                phlist.append(phv)
            return phlist
        elif isinstance(index, int):
            i = self.indPer[index % self.indPer._npoint] * 3
            phv = AlgVector.PhantomVector()
            phv._v = &self.vectMat._m[i]
            return phv

    def __len__(self):
        return self.indPer._npoint


cdef class SegmentIterator():

    def __init__(self, vectMat:AlgMatrix.Matrix, indPer:IndexPerimeter):
        self.vectMat = vectMat # keep the vector point
        self.indPer = indPer # keep the secuential index

    # def __dealloc__(self):
    #     pass

    def __iter__(self):
        self.niterator = 0
        return self

    def __next__(self):
        cdef AlgSegment.PhantomSegment v
        if self.niterator < <unsigned int>(<int>(self.indPer._npoint) - 1):
            v = AlgSegment.PhantomSegment()
            v._vini = &self.vectMat._m[self.indPer._m[self.niterator] * 3]
            v._vfin = &self.vectMat._m[self.indPer._m[self.niterator + 1] * 3]
            self.niterator += 1
            return v
        else:
            raise StopIteration

    def __getitem__(self, index):
        cdef unsigned int i
        cdef int start, stop, step
        cdef AlgSegment.PhantomSegment phs
        if isinstance(index, slice):
            start = 0 if index.start is None else index.start if index.start >= 0 else <int>self.indPer._npoint + index.start - 1
            stop = self.indPer._npoint if index.stop is None else index.stop if index.stop >= 0 else self.indPer._npoint + index.stop - 1
            step = 1 if index.step is None else index.step
            phslist = []
            for i in range(start, stop, step):
                phs = AlgSegment.PhantomSegment()
                phs._vini = &self.vectMat._m[self.indPer._m[i]*3]
                phs._vfin = &self.vectMat._m[self.indPer._m[i+1]*3]
                phslist.append(phs)
            return phslist
        elif isinstance(index, int):
            phs = AlgSegment.PhantomSegment()
            phs._vini = &self.vectMat._m[self.indPer._m[index % self.indPer._npoint]*3]
            phs._vfin = &self.vectMat._m[self.indPer._m[(index+1) % self.indPer._npoint]*3]
            return phs

    def __len__(self):
        return self.indPer._npoint - 1



cdef class Wire:

    def __cinit__(self, pointList=[], **kwargs):
        # self.phantom = <bint>kwargs.get('phantom', False)
        self.vectMat = AlgMatrix.Matrix()
        self.indPer = IndexPerimeter()

    def __init__(self, pointList=[], **kwargs):
        cdef unsigned int i, j
        if pointList:
            if isinstance(pointList, (list, tuple)):
                self.vectMat._m = <double *> PyMem_Malloc(len(pointList) * 3 * sizeof(double))
                self.vectMat._rows = <unsigned int> len(pointList)
                self.vectMat._cols = 3
                for i in range(len(pointList)):
                    for j in range(3):
                        self.vectMat._m[i*3 + j] = pointList[i][j]
                self.indPer.setList(list(range(len(pointList))))
            elif isinstance(pointList, AlgMatrix.Matrix):
                self.vectMat._m = <double *> PyMem_Malloc((<AlgMatrix.Matrix> pointList)._rows * (<AlgMatrix.Matrix> pointList)._cols * sizeof(double))
                self.vectMat._rows = (<AlgMatrix.Matrix> pointList)._rows
                self.vectMat._cols = (<AlgMatrix.Matrix> pointList)._cols
                for i in range((<AlgMatrix.Matrix> pointList)._rows * (<AlgMatrix.Matrix> pointList)._cols):
                    self.vectMat._m[i] = (<AlgMatrix.Matrix> pointList)._m[i]
                self.indPer.setList(list(range(self.vectMat._rows)))

    def __dealloc__(self):
        pass

    #.........................C-METHODS.................................

    def __json__(self):
        return {'__jsoncls__': 'CLinAlg.AlgWire:Wire.from_JSON', 'mat': self.vectMat, 'indexes': self.indPer}

    @classmethod
    def from_JSON(cls, jsondict):
        obj = Wire()
        obj.vectMat = jsondict['mat']
        obj.indPer = jsondict['indexes']
        return obj

    #........................PYTHON METHODS.............................


    @property
    def Point(self):
        return PointIterator(self.vectMat, self.indPer)

    @property
    def Segment(self):
        return SegmentIterator(self.vectMat, self.indPer)

    def indexPerimeter(self):
        return self.indPer.__serialize__()

    def perimeter(self):
        return getLength(self.vectMat._m, self.indPer._m, self.indPer._npoint)

    cpdef double length(self):
        return getLength(self.vectMat._m, self.indPer._m, self.indPer._npoint)

    cpdef list getBounds(self):
        cdef list bound = [sys.float_info.max, -sys.float_info.max, sys.float_info.max, -sys.float_info.max, sys.float_info.max, -sys.float_info.max]
        cdef int j
        for p in self.Point:
            for j in range(3):
                bound[j*2] = min(bound[j*2], p[j])
                bound[j*2+1] = max(bound[j*2+1], p[j])
        return bound

    def isInside(self, point:Vector, bint incEdges=True):
        return isInside(self.vectMat._m, self.indPer._m, self.indPer._npoint, point._v, <bint>incEdges)

    def cdg(self):
        cdef Vector v = Vector()
        getCDG(v._v, self.vectMat._m, self.indPer._m, self.indPer._npoint)
        return v

    def normal(self):
        cdef Vector v = Vector()
        getNormal(v._v, self.vectMat._m, self.indPer._m, self.indPer._npoint)
        return v

    def isClosed(self):
        return isClosed(self.vectMat._m, self.indPer._m, self.indPer._npoint)


    def getRotate(self, center:Vector, axis:Vector, double angle):
        cdef double * toPointer = <double *> PyMem_Malloc (self.vectMat._rows*self.vectMat._cols*sizeof(double))
        cdef Wire newWire = Wire()
        getRotate(toPointer, self.vectMat._m, self.vectMat._rows, center._v, axis._v, angle)
        newWire.vectMat._m = toPointer
        newWire.vectMat._rows = self.vectMat._rows
        newWire.vectMat._cols = 3
        newWire.indPer.setList(self.indPer.__serialize__())
        return newWire

    def getPointFromOrigin(self, length):
        cdef Vector newvect = Vector()
        if getPointFromOrigin(newvect._v, self.vectMat._m, self.indPer._m, self.indPer._npoint, length):
            return newvect


    def getLengthFromOrigin(self, point:Vector):
        return getLengthFromOrigin(self.vectMat._m, self.indPer._m, self.indPer._npoint, point._v)

    def getDirectionAtPoint(self, point:Vector):
        cdef Vector newVect = Vector()
        if getDirectionAtPoint(newVect._v, self.vectMat._m, self.indPer._m, self.indPer._npoint, point._v):
            return newVect

    def getNearestPoint(Wire self, Vector point):
        cdef Vector newVect = Vector()
        getNearestPoint(newVect._v, self.vectMat._m, self.indPer._m, self.indPer._npoint, point._v)
        return newVect

    def isClockWise(self, obsPoint:Vector):
        return isClockWise(self.vectMat._m, self.indPer._m, self.indPer._npoint, obsPoint._v)

    def getShifted(self, value):
        "Shift a wire by the given value (+) to inner (-) to outer"
        cdef Wire newWire = Wire(getShifted(self.vectMat._m, self.indPer._m, self.indPer._npoint, value))
        return newWire

    def getIntersectionPoint(self, other, incEdge=False):
        if isinstance(other, Wire):
            return getIntersectionPointWithWire(self.vectMat._m, self.indPer._m, self.indPer._npoint,
                                                (<Wire>other).vectMat._m, (<Wire>other).indPer._m,(<Wire>other).indPer._npoint,
                                                <bint>incEdge)
        elif isinstance(other, AlgLine.Line):
            return getIntersectionPointWithLine(self.vectMat._m, self.indPer._m, self.indPer._npoint,
                                                (<AlgLine.Line>other)._pnt, (<AlgLine.Line>other)._dir,
                                                <bint>incEdge)
        elif isinstance(other, AlgSegment.Segment):
            return getIntersectionPointWithSegment(self.vectMat._m, self.indPer._m, self.indPer._npoint,
                                                   (<AlgSegment.Segment>other)._vini, (<AlgSegment.Segment>other)._vfin,
                                                   <bint>incEdge)
        else:
            raise ValueError("Intersection point arguments are wrong")

    def simplify(self):
        "Simplifica el modelo"
        cdef int i, nexti, previ, nump
        if self.isClosed():
            nump = (<int>self.indPer._npoint - 1)
        else:
            nump = <int>self.indPer._npoint
        i = 1
        while i < nump:
            previ = <int> ((i - 1) % nump)
            nexti = <int> ((i + 1) % nump)
            if AlgSegment.isInside(&self.vectMat._m[previ * 3], &self.vectMat._m[nexti * 3],
                                        &self.vectMat._m[i * 3], <bint>True):
                self.vectMat.deleteRow(i)
                nump -= 1
                continue
            i += 1

    def getSubWire(self, length1, length2):
        return getSubWire(self.vectMat._m,self.indPer._m, self.indPer._npoint, length1, length2)

    def copy(Wire self) -> Wire:
        return Wire(self.vectMat)

    cpdef Wire transform(Wire self, Transformation transf):
        "function that apply transformation matrix to current points. Returns a new wire transformed."
        cdef unsigned int v, i, j
        cdef double val
        cdef Wire newWire = Wire()
        newWire.vectMat._m = <double *> PyMem_Malloc (self.vectMat._rows * self.vectMat._cols * sizeof(double))
        newWire.vectMat._rows = self.vectMat._rows
        newWire.vectMat._cols = self.vectMat._cols
        for v in range(self.vectMat._rows):
            for i in range(3):
                val = 0
                for j in range(3):
                    val += transf._mat._m[j+4*i] * self.vectMat._m[j + v*3]
                val += transf._mat._m[3+4*i]
                newWire.vectMat._m[i+v*3] = val
        newWire.indPer.setList(self.indPer.__serialize__())
        return newWire




    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memodict={}):
        return self.copy()

    def __repr__(self):
        cdef double * cdg = <double *> PyMem_Malloc ( 3 * sizeof(double))
        cdef double perim = getLength(self.vectMat._m, self.indPer._m, self.indPer._npoint)
        getCDG(cdg, self.vectMat._m, self.indPer._m, self.indPer._npoint)
        try:
            return "Wire(Perimeter=%.2f; CDG=(%.2f,%.2f,%.2f))" % (perim, cdg[0], cdg[1], cdg[2])
        finally:
            PyMem_Free(cdg)

    def __str__(self):
        cdef double * cdg = <double *> PyMem_Malloc ( 3 * sizeof(double))
        cdef double perim = getLength(self.vectMat._m, self.indPer._m, self.indPer._npoint)
        getCDG(cdg, self.vectMat._m, self.indPer._m, self.indPer._npoint)
        try:
            return "Wire(Perimeter=%.2f; CDG=(%.2f,%.2f,%.2f))" % (perim, cdg[0], cdg[1], cdg[2])
        finally:
            PyMem_Free(cdg)
