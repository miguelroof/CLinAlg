import sys
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from libc.math cimport M_PI, sin, fabs
from .AlgVector cimport Vector, PhantomVector
from . cimport AlgVector, AlgSegment, AlgQuaternion, AlgMatrix, AlgLine
from .AlgTool cimport presition


cdef double getLength(double * wirepoint, unsigned int numvect):
    cdef double length = 0
    cdef unsigned int i
    for i in range(<unsigned int>(numvect - 1)):
        length += AlgVector.distance(&wirepoint[i * 3], &wirepoint[(i + 1) * 3])
    return length

cdef double getLengthByIndex(double * wirepoint, list wireIndex):
    cdef double length = 0
    cdef int i 
    for i in range(len(wireIndex)-1):
        length += AlgVector.distance(&wirepoint[<int>wireIndex[i] * 3], &wirepoint[<int>wireIndex[i+1]*3])
    return length
        

cdef bint isInside(double * wirepoint, unsigned int numvect, double * punto, bint incEdges):
    cdef unsigned int i
    for i in range(<unsigned int> (numvect - 1)):
        if AlgSegment.isInside(&wirepoint[i * 3], &wirepoint[(i + 1) * 3], punto, incEdges):
            return True
    return False

cdef bint isClosed(double * wirePoint, int numpoints):
    return all([fabs(wirePoint[i] - wirePoint[i + (numpoints - 1) * 3]) < presition for i in range(3)])

cdef bint isInsideByIndex(double * wirepoint, list wireIndex, double * punto, bint incEdges):
    cdef int i 
    for i in range(len(wireIndex)-1):
        if AlgSegment.isInside(&wirepoint[<int>wireIndex[i]*3], &wirepoint[<int>wireIndex[i+1]*3], punto, incEdges):
            return True
    return False

cdef void getRotate(double * toPointer, double * wirePoint, unsigned int numpoints, double * center, double * axis,
                            double angle):
    cdef double * quat = <double *> PyMem_Malloc (4 * sizeof(double))
    AlgQuaternion.quatFromAxisAngle(quat, axis, angle)
    AlgQuaternion.rotateVectorMatrix(toPointer, quat, center, wirePoint, numpoints)
    PyMem_Free(quat)

cdef void getCDG(double * cdg, double * wirePoint, unsigned int numpoints):
    cdef unsigned int i, j
    cdef unsigned int numrows = numpoints
    if isClosed(wirePoint, numpoints):
        numrows = numrows -1
    for i in range(3):
        cdg[i] = 0
        for j in range(numrows):
            cdg[i] += wirePoint[i + j * 3]
        cdg[i] = cdg[i] / numrows

cdef void getCDGByIndex(double * cdg, double * wirePoint, list wireInd):
    cdef int i, j, numrows
    if AlgVector.isEqual(&wirePoint[<int>wireInd[0] * 3], &wirePoint[<int>wireInd[-1] * 3]):
        numrows = <unsigned int> (len(wireInd) - 1)
    else:
        numrows = <unsigned int> (len(wireInd))
    for i in range(3):
        cdg[i] = 0
        for j in range(numrows):
            cdg[i] += wirePoint[i + <int>wireInd[j] * 3]
        cdg[i] = cdg[i] / numrows

cdef bint getPointFromOrigin(double * toPoint, double * wirePoint, unsigned int numpoints, double length):
    cdef unsigned int i, j
    cdef double tlen = 0
    cdef double prop, slen
    cdef AlgSegment.PhantomSegment seg
    if length <= 0:
        for i in range(3):
            toPoint[i] = wirePoint[i]
        return True
    elif length > getLength(wirePoint, numpoints):
        for i in range(3):
            toPoint[i] = wirePoint[<unsigned int> (numpoints - 1) * 3 + i]
        return True
    for i in range(<unsigned int>(numpoints-1)):
        slen = AlgVector.distance(&wirePoint[i*3], &wirePoint[(i+1)*3])
        if slen and tlen + slen >= length:
            prop = (length - tlen) / slen
            for j in range(3):
                toPoint[j] = wirePoint[i*3+j] + prop * (wirePoint[(i+1)*3+j] - wirePoint[i*3+j])
            return True
        tlen += slen
    return False

cdef double getLengthFromOrigin(double * wirePoint, unsigned int numpoints, double * point):
    cdef double dist = 0
    cdef double ldist
    cdef unsigned int i
    for i in range(<unsigned int> (numpoints - 1)):
        if AlgSegment.isInside(&wirePoint[i * 3],
                                    &wirePoint[(i + 1) * 3],
                                    point, <bint>True):
            dist += AlgVector.distance(&wirePoint[i * 3], point)
            return dist
        else:
            dist += AlgVector.distance(&wirePoint[i * 3],
                                 &wirePoint[(i + 1) * 3])
    return sys.float_info.max

cdef bint getDirectionAtPoint(double * toDir, double * wirePoint, unsigned int numpoints, double * point):
    cdef int i
    for i in range(numpoints-1):
        if AlgVector.isEqual(&wirePoint[i*3], point) or AlgSegment.isInside(&wirePoint[i*3], &wirePoint[(i+1)*3], point, <bint>False):
            AlgVector.vdir(toDir, &wirePoint[i*3], &wirePoint[(i+1)*3])
            return True
    return False

cdef void getNearestPoint(double * toPoint, double * wirePoint, unsigned int numpoints, double * point):
    cdef unsigned int i, j
    cdef double val, curval
    cdef double * distpoint = <double *> PyMem_Malloc(3 * sizeof(double))
    val = sys.float_info.max
    for i in range(<unsigned int> (numpoints - 1)):
        AlgSegment.getNearestPoint(distpoint, &wirePoint[i * 3],
                                                    &wirePoint[(i + 1) * 3],
                                                    point)
        curval = AlgVector.distance(distpoint, point)
        if curval < val:
            val = curval
            AlgVector.copy(toPoint, distpoint)
    PyMem_Free(distpoint)

cdef bint isClockWise(double * wirePoint, unsigned int numpoints, double * obsPoint):
    "Determine if points are clockwise in reference to an observer point"
    cdef double * cdg = <double *> PyMem_Malloc (3 * sizeof(double))
    cdef double * v = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * vt1 = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * vt2 = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * vt3 = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef unsigned int i, j
    cdef double modulo, angle
    getCDG(cdg, wirePoint, numpoints)
    try:
        for j in range(3):
            v[j] = 0
        for i in range(<unsigned int> (numpoints - 1)):
            AlgVector.sub(vt1, &wirePoint[(i + 1) * 3], &wirePoint[i * 3])
            AlgVector.sub(vt2, &wirePoint[i * 3], cdg)
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


cdef list getShifted(double * wirePoint, unsigned int numpoints, double value):
    cdef unsigned int i, nump, iprev, inext, j
    cdef double * vt1 = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * vt2 = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * vt3 = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * vt4 = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * obs = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * cdg = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double ang, shp
    getCDG(cdg, wirePoint, numpoints)
    try:
        vlist = []
        if numpoints == 2:
            for i in range(2):
                vlist.append(Vector(wirePoint[i * 3], wirePoint[i * 3 + 1], wirePoint[i * 3 + 2]))
            return vlist
        for j in range(3):
            obs[j] = 0
        for i in range(<unsigned int> (numpoints - 1)):
            AlgVector.sub(vt1, &wirePoint[(i+1)*3], &wirePoint[i*3])
            AlgVector.sub(vt2, &wirePoint[i*3], cdg)
            AlgVector.cross(vt3, vt1, vt2)
            AlgVector.add(obs, vt3, obs)
        shp = AlgVector.module(obs)
        for j in range(3):
            obs[j] = obs[j] / shp
        if isClosed(wirePoint, numpoints):
            nump = numpoints - 1
            iprev = <unsigned int> ((-1) % nump)
            for i in range(nump):
                inext = <unsigned int> ((i + 1) % nump)
                AlgVector.vdir(vt1, &wirePoint[iprev*3], &wirePoint[i*3])
                AlgVector.vdir(vt2, &wirePoint[i*3], &wirePoint[inext*3])
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
                    vlist.append(Vector(wirePoint[i * 3] + vt3[0], wirePoint[i * 3 + 1] + vt3[1],
                                        wirePoint[i * 3 + 2] + vt3[2]))
                else:
                    vlist.append(Vector(wirePoint[i * 3] - vt3[0], wirePoint[i * 3 + 1] - vt3[1],
                                        wirePoint[i * 3 + 2] - vt3[2]))
                iprev = i
            vlist.append(vlist[0].copy())
        else:
            nump = numpoints
            iprev = 0
            for i in range(1, <unsigned int> (nump - 1)):
                inext = <unsigned int> ((i + 1) % nump)
                AlgVector.vdir(vt1, &wirePoint[3*iprev], &wirePoint[3*i])
                AlgVector.vdir(vt2, &wirePoint[3*i], &wirePoint[3*inext])
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
                        vlist.append(Vector(wirePoint[0] + vt2[0], wirePoint[1] + vt2[1],
                                            wirePoint[2] + vt2[2]))
                    else:
                        vlist.append(Vector(wirePoint[0] - vt2[0], wirePoint[1] - vt2[1],
                                            wirePoint[2] - vt2[2]))
                if ang <= M_PI / 2:
                    vlist.append(Vector(wirePoint[i * 3] + vt3[0], wirePoint[i * 3 + 1] + vt3[1],
                                        wirePoint[i * 3 + 2] + vt3[2]))
                else:
                    vlist.append(Vector(wirePoint[i * 3] - vt3[0], wirePoint[i * 3 + 1] - vt3[1],
                                        wirePoint[i * 3 + 2] - vt3[2]))
                iprev = i
            #agrego el ultimo
            AlgVector.cross(vt1, obs, vt2)
            shp = AlgVector.module(vt1)
            for j in range(3):
                vt1[j] = vt1[j] * value / shp
            iprev = <unsigned int> (numpoints- 1)
            if ang > M_PI / 2:
                vlist.append(Vector(wirePoint[iprev * 3] + vt1[0], wirePoint[iprev * 3 + 1] + vt1[1],
                                    wirePoint[iprev * 3 + 2] + vt1[2]))
            else:
                vlist.append(Vector(wirePoint[iprev * 3] - vt1[0], wirePoint[iprev * 3 + 1] - vt1[1],
                                    wirePoint[iprev * 3 + 2] - vt1[2]))
        return vlist
    finally:
        PyMem_Free(vt1)
        PyMem_Free(vt2)
        PyMem_Free(vt3)
        PyMem_Free(vt4)
        PyMem_Free(obs)
        PyMem_Free(cdg)

cdef list getIntersectionPointWithWire(double * wirePoint, unsigned int numpoints, 
                                       double * otherWirePoint, unsigned int otherNumpoints, 
                                       bint incEdge):
    "Return the intersection point between this wire and another wire"
    cdef int i, j
    cdef double * ipoint
    cdef double * tpoint = <double *> PyMem_Malloc (3 * sizeof(double))
    cdef bint sclosed = isClosed(wirePoint, numpoints)
    cdef bint oclosed = isClosed(otherWirePoint, otherNumpoints)
    cdef double * dirMIn = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * dirMOut = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * dirOIn = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * dirOOut = <double *> PyMem_Malloc(3 * sizeof(double))
    try:
        cortes = []
        for i in range(numpoints - 1):
            for j in range(otherNumpoints - 1):
                ipoint = <double *> PyMem_Malloc (3 * sizeof(double))
                if not AlgSegment.getIntersectionPointWithSegment(ipoint, &wirePoint[i * 3],
                                                                         &wirePoint[(i + 1) * 3],
                                                                         &otherWirePoint[j * 3],
                                                                         &otherWirePoint[(j + 1) * 3],
                                                                         <bint>True) or any([AlgVector.isEqual(ipoint, (<Vector> v)._v) for v in cortes]):

                    PyMem_Free(ipoint)
                    continue
                if (incEdge or (not AlgVector.isEqual(ipoint, &wirePoint[i * 3])
                                and not AlgVector.isEqual(ipoint, &wirePoint[(i + 1) * 3])
                                and not AlgVector.isEqual(ipoint, &otherWirePoint[j * 3])
                                and not AlgVector.isEqual(ipoint, &otherWirePoint[(j + 1) * 3]))):
                    newVect = Vector(None)
                    newVect._v = ipoint
                    cortes.append(newVect)
                    continue

                #vectors iniciales
                if AlgVector.isEqual(ipoint, &wirePoint[i * 3]):  # punto inicial
                    if i == 0:
                        if not sclosed:
                            PyMem_Free(ipoint)
                            continue
                        AlgVector.vdir(dirMIn, &wirePoint[(numpoints - 1) * 3],
                              &wirePoint[(numpoints - 2) * 3])
                    else:
                        AlgVector.vdir(dirMIn, &wirePoint[i * 3], &wirePoint[(i - 1) * 3])
                    AlgVector.vdir(dirMOut, &wirePoint[i * 3], &wirePoint[(i + 1) * 3])
                elif AlgVector.isEqual(ipoint, &wirePoint[(i + 1) * 3]):  #punto final
                    if i == (numpoints - 2):
                        if not sclosed:
                            PyMem_Free(ipoint)
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
                            PyMem_Free(ipoint)
                            continue
                        AlgVector.vdir(dirOIn, &otherWirePoint[(otherNumpoints - 1) * 3],
                              &otherWirePoint[(otherNumpoints - 2) * 3])
                    else:
                        AlgVector.vdir(dirOIn, &otherWirePoint[j * 3], &otherWirePoint[(j - 1) * 3])
                    AlgVector.vdir(dirOOut, &otherWirePoint[j * 3], &otherWirePoint[(j + 1) * 3])
                elif AlgVector.isEqual(ipoint, &otherWirePoint[(j + 1) * 3]):  #punto final
                    if j == (otherNumpoints - 2):
                        if not oclosed:
                            PyMem_Free(ipoint)
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
                else:
                    PyMem_Free(ipoint)
        return cortes
    finally:
        PyMem_Free(dirMIn)
        PyMem_Free(dirMOut)
        PyMem_Free(dirOIn)
        PyMem_Free(dirOOut)
        PyMem_Free(tpoint)

cdef list getIntersectionPointWithLine(double * wirePoint, unsigned int numpoints, double * linePoint, double * lineDir, bint incEdge):
    "Return the intersection point between this wire and one line"
    puntos = list(range(numpoints))
    return getIntersectionPointWithLineByIndex(wirePoint, <list>puntos, linePoint, lineDir, incEdge)

cdef list getIntersectionPointWithLineByIndex(double * basePoints, list wireInd, double * linePoint, double * lineDir,
                                           bint incEdge):
    # Funcion que genera el mismo algoritmo que BaseWire.c_getINtersecitonPointWidthLine pero usando un secuenciador de puntos
    cdef int i, n, pp, pc, pn
    cdef double * ipoint
    cdef double * tpoint = <double *> PyMem_Malloc (3 * sizeof(double))
    cdef bint sclosed
    cdef double * dirMIn = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * dirMOut = <double *> PyMem_Malloc(3 * sizeof(double))
    n = <int> (len(wireInd))
    if all([fabs(basePoints[<int> (wireInd[i])] - basePoints[i + (<int> (wireInd[-1])) * 3]) < presition for i in
            range(3)]):
        sclosed = True
    else:
        sclosed = False
    try:
        cortes = []
        for i in range(n - 1):
            pc = <int> (wireInd[i])
            pn = <int> (wireInd[i + 1])
            ipoint = <double *> PyMem_Malloc(3 * sizeof(double))
            if not AlgSegment.getIntersectionPointWithLine(ipoint, &basePoints[pc * 3],
                                                                  &basePoints[pn * 3],
                                                                  linePoint,
                                                                  lineDir,
                                                                  incEdge) or any([AlgVector.isEqual(ipoint, (<Vector> v)._v) for v in cortes]):

                PyMem_Free(ipoint)
                continue
            if (incEdge or (not AlgVector.isEqual(ipoint, &basePoints[pc * 3])
                            and not AlgVector.isEqual(ipoint, &basePoints[pn * 3]))):
                newVect = Vector(None)
                newVect._v = ipoint
                cortes.append(newVect)
                continue

            #vectors iniciales. Condicion a la que salto porque me ha coincidido sobre un vector
            if AlgVector.isEqual(ipoint, &basePoints[pc * 3]):  # punto inicial
                if i == 0:
                    if not sclosed:
                        PyMem_Free(ipoint)
                        continue
                    pp = <int> (wireInd[-2])
                    AlgVector.vdir(dirMIn, &basePoints[pc * 3], &basePoints[pp * 3])
                else:
                    pp = <int> (wireInd[i - 1])
                    AlgVector.vdir(dirMIn, &basePoints[pc * 3], &basePoints[pp * 3])
                AlgVector.vdir(dirMOut, &basePoints[pc * 3], &basePoints[pn * 3])
                AlgVector.add(dirMIn, &basePoints[pc * 3], dirMIn)
                AlgVector.add(dirMOut, &basePoints[pc * 3], dirMOut)
            elif AlgVector.isEqual(ipoint, &basePoints[pn * 3]):  #punto final
                if i == n - 2:
                    if not sclosed:
                        PyMem_Free(ipoint)
                        continue
                    pp = <int> (wireInd[1])
                    AlgVector.vdir(dirMOut, &basePoints[pn * 3], &basePoints[pp * 3])
                else:
                    pp = <int> (wireInd[i + 2])
                    AlgVector.vdir(dirMOut, &basePoints[pn * 3], &basePoints[pp * 3])
                AlgVector.vdir(dirMIn, &basePoints[pn * 3], &basePoints[pc * 3])
                AlgVector.add(dirMIn, &basePoints[pn * 3], dirMIn)
                AlgVector.add(dirMOut, &basePoints[pn * 3], dirMOut)
            else:  #ninguna de las anteriores
                PyMem_Free(ipoint)
                raise ValueError('No es posible esta combinacion')
            if AlgSegment.getIntersectionPointWithLine(tpoint, dirMIn, dirMOut, linePoint,
                                                                  lineDir, <bint>False):
                newVect = Vector(None)
                newVect._v = ipoint
                cortes.append(newVect)
            else:
                PyMem_Free(ipoint)
        return cortes
    finally:
        PyMem_Free(dirMIn)
        PyMem_Free(dirMOut)
        PyMem_Free(tpoint)


cdef list getIntersectionPointWithSegment(double * wirePoint, unsigned int numpoints, 
                                          double * vini, double * vfin,
                                          bint incEdge):
    "Return the intersection point between this wire and another wire"
    cdef double * segDir = <double *> PyMem_Malloc (3 * sizeof(double))
    try:
        AlgVector.vdir(segDir, vini, vfin)
        linecortes = getIntersectionPointWithLine(wirePoint, numpoints, vini, segDir, incEdge)
        cortes = []
        for i in range(len(linecortes)):
            if AlgSegment.isInside(vini, vfin, (<Vector>linecortes[i])._v, incEdge):
                cortes.append(linecortes[i])
        return cortes
    finally:
        PyMem_Free(segDir)


cdef Wire getSubWire(double * wirePoint, unsigned int numpoints, double length1, double length2):
    "Funcion que devuelve un subwire desde los puntos inicial y final. Destruye ciclos cerrados"
    cdef bint revert = False
    cdef double _l = 0
    cdef double elength, prop
    cdef unsigned int i, j
    if length2 < length1:
        revert = True
        length1, length2 = length2, length1
    length1 = min(<double>0, length1)
    length2 = max(getLength(wirePoint, numpoints), length2)
    puntos = []
    for i in range(0, <unsigned int> (numpoints - 1)):
        elength = AlgVector.distance(&wirePoint[i * 3], &wirePoint[(i + 1) * 3])
        if elength < 2 * presition:
            continue
        if _l > length2:
            break
        if _l - presition <= length1 < _l + elength - presition:  # estamos en el punto inicial
            prop = max((length1 - _l) / elength, 0)
            puntos.append(
                [wirePoint[i * 3 + j] + prop * (wirePoint[(i + 1) * 3 + j] - wirePoint[i * 3 + j])
                 for j in range(3)])
        if length1 + presition < _l < length2 - presition:
            puntos.append([wirePoint[i * 3 + j] for j in range(3)])
        if _l + presition < length2 <= _l + elength + presition:  # estamos en el punto final
            propfin = min(1, <int>((length2 - _l) / elength))
            puntos.append(
                [wirePoint[i * 3 + j] + prop * (wirePoint[(i + 1) * 3 + j] - wirePoint[i * 3 + j])
                 for j in range(3)])
            break
        _l += elength
    if not puntos:
        return
    if revert:
        return Wire(list(reversed(puntos)))
    else:
        return Wire(puntos)

cdef void getNormal(double * normal, double * wirePoint, unsigned int numpoints):
    "Funcion que devuelve la normal de un wire. Siempre define la direccion contrarreloj"
    cdef double * cdg = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * vtemp1 = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * vtemp2 = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * vtemp3 = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double modulo
    cdef int i, j
    try:
        getCDG(cdg, wirePoint, numpoints)
        if numpoints <= 2:
            for i in range(3):
                normal[i] = 0
            return
        for i in range(3):
            normal[i] = 0
        # eje principal
        for i in range(numpoints - 1):
            AlgVector.sub(vtemp1, &wirePoint[(i+1) * 3], &wirePoint[i * 3])
            AlgVector.sub(vtemp2, &wirePoint[i * 3], cdg)
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

cdef void getNormalByIndex(double * normal, double * wirePoint, list wireInd):
    "Funcion que devuelve la normal de un wire. Siempre define la direccion contrarreloj"
    cdef double * cdg = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * vtemp1 = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * vtemp2 = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * vtemp3 = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double modulo
    cdef int i, j, numrows
    try:
        getCDGByIndex(cdg, wirePoint, wireInd)
        if AlgVector.isEqual(&wirePoint[<int> (wireInd[0]) * 3], &wirePoint[<int> (wireInd[-1]) * 3]):
            # numrows = <unsigned int> (len(wireInd) - 1)
            wlist = [<int> (p) for p in <list>wireInd]
        else:
            wlist = [<int> (p) for p in <list>wireInd] + [<int> (wireInd[0])]
            # numrows = <unsigned int> (len(wireInd))
        numrows = <int> (len(wlist) - 1)
        for i in range(3):
            normal[i] = 0
        # eje principal
        for i in range(numrows - 1):
            AlgVector.sub(vtemp1, &wirePoint[wlist[i + 1] * 3], &wirePoint[wlist[i] * 3])
            AlgVector.sub(vtemp2, &wirePoint[wlist[i] * 3], cdg)
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

cdef class PointIterator():

    def __cinit__(self, mat:AlgMatrix.Matrix):
        self.mat = mat

    def __iter__(PointIterator self):
        self.n = 0
        return self

    def __next__(PointIterator self):
        cdef PhantomVector v
        if self.n < <int>self.mat._rows:
            v = PhantomVector()
            v._v = &self.mat._m[self.n * 3]
            self.n += 1
            return v
        else:
            raise StopIteration

    def __getitem__(self, index):
        cdef unsigned int i
        cdef int start, stop, step
        if isinstance(index, slice):
            start = 0 if index.start is None else index.start if index.start >= 0 else self.mat._rows + index.start
            stop = self.mat._rows if index.stop is None else index.stop if index.stop >= 0 else self.mat._rows + index.stop
            step = 1 if index.step is None else index.step
            phlist = []
            for i in range(start, stop, step):
                phv = PhantomVector()
                phv._v = &self.mat._m[i*3]
                phlist.append(phv)
            return phlist
        elif isinstance(index, int):
            i = (index % self.mat._rows) * 3
            phv = PhantomVector()
            phv._v = &self.mat._m[i]
            return phv

    def __len__(self):
        return self.mat._rows


cdef class SegmentIterator():

    def __cinit__(self, mat:AlgMatrix.Matrix):
        self.mat = mat

    def __iter__(SegmentIterator self):
        self.n = 0
        return self

    def __next__(SegmentIterator self):
        cdef AlgSegment.PhantomSegment seg
        if self.n < <int>(self.mat._rows - 1):
            seg = AlgSegment.PhantomSegment()
            seg._vini = &self.mat._m[self.n * 3]
            seg._vfin = &self.mat._m[(self.n + 1) * 3]
            self.n += 1
            return seg
        else:
            raise StopIteration

    def __getitem__(SegmentIterator self, index):
        cdef unsigned int i
        cdef int start, stop, step
        if isinstance(index, slice):
            start = 0 if index.start is None else index.start if index.start >= 0 else self.mat._rows + index.start - 1
            stop = self.mat._rows if index.stop is None else index.stop if index.stop >= 0 else self.mat._rows + index.stop - 1
            step = 1 if index.step is None else index.step
            phslist = []
            for i in range(start, stop, step):
                phs = AlgSegment.PhantomSegment()
                (<AlgSegment.PhantomSegment>phs)._vini = &self.mat._m[i*3]
                (<AlgSegment.PhantomSegment>phs)._vfin = &self.mat._m[(i+1)*3]
                phslist.append(phs)
            return phslist
        elif isinstance(index, int):
            i = (index % (self.mat._rows - 1)) * 3
            phs = AlgSegment.PhantomSegment()
            (<AlgSegment.PhantomSegment>phs)._vini = &self.mat._m[i]
            (<AlgSegment.PhantomSegment>phs)._vfin = &self.mat._m[i+3]
            return phs

    def __len__(PointIterator self):
        return self.mat._rows - 1




cdef class Wire:

    def __cinit__(self, pointList=[]):
        cdef unsigned int i, j
        self.vectMat = AlgMatrix.Matrix()
        if pointList:
            if isinstance(pointList, (list, tuple)):
                self.vectMat._m = <double *> PyMem_Malloc(len(pointList) * 3 * sizeof(double))
                self.vectMat._rows = <unsigned int> len(pointList)
                self.vectMat._cols = 3
                for i in range(len(pointList)):
                    for j in range(3):
                        self.vectMat._m[i*3 + j] = pointList[i][j]
            elif isinstance(pointList, AlgMatrix.Matrix):
                self.vectMat._m = <double *> PyMem_Malloc((<AlgMatrix.Matrix> pointList)._rows * (<AlgMatrix.Matrix> pointList)._cols * sizeof(double))
                self.vectMat._rows = (<AlgMatrix.Matrix> pointList)._rows
                self.vectMat._cols = (<AlgMatrix.Matrix> pointList)._cols
                for i in range((<AlgMatrix.Matrix> pointList)._rows * (<AlgMatrix.Matrix> pointList)._cols):
                    self.vectMat._m[i] = (<AlgMatrix.Matrix> pointList)._m[i]

    #.........................C-METHODS.................................

    def __json__(self):
        return {'__jsoncls__': 'CLinAlg.AlgWire:Wire.from_JSON', 'mat': self.vectMat}

    @classmethod
    def from_JSON(cls, jsondict):
        obj = cls(jsondict['mat'])
        return obj

    #........................PYTHON METHODS.............................

    @property
    def mat(self):
        return self.vectMat

    @property
    def Vertexes(self):
        return PointIterator(self.vectMat)

    @property
    def Edges(self):
        return SegmentIterator(self.vectMat)

    @property
    def Perimeter(self):
        return getLength(self.vectMat._m, self.vectMat._rows)

    cpdef double length(self):
        return getLength(self.vectMat._m, self.vectMat._rows)

    def isInside(self, point:Vector, bint incEdges=True):
        return isInside(self.vectMat._m, self.vectMat._rows, point._v, <bint>incEdges)

    @property
    def CDG(self):
        cdef Vector v = Vector(None)
        v._v = <double *> PyMem_Malloc (3 * sizeof(double))
        getCDG(v._v, self.vectMat._m, self.vectMat._rows)
        return v

    def isClosed(self):
        return isClosed(self.vectMat._m, self.vectMat._rows)

    @property
    def closed(self) -> bint:
        return self.isClosed()

    def getRotate(self, center:Vector, axis:Vector, double angle):
        cdef double * toPointer = <double *> PyMem_Malloc (self.vectMat._rows*self.vectMat._cols*sizeof(double))
        cdef Wire newWire = Wire()
        getRotate(toPointer, self.vectMat._m, self.vectMat._rows, center._v, axis._v, angle)
        newWire.vectMat._m = toPointer
        newWire.vectMat._rows = self.vectMat._rows
        newWire.vectMat._cols = 3
        return newWire

    def getPointFromOrigin(self, length):
        cdef Vector newvect = Vector(None)
        newvect._v = <double *> PyMem_Malloc (3 * sizeof(double))
        if getPointFromOrigin(newvect._v, self.vectMat._m, self.vectMat._rows, length):
            return newvect


    def getLengthFromOrigin(self, point:Vector):
        return getLengthFromOrigin(self.vectMat._m, self.vectMat._rows, point._v)

    def getDirectionAtPoint(self, point:Vector):
        cdef Vector newVect = Vector(None)
        newVect._v = <double *> PyMem_Malloc (3 * sizeof(double))
        if getDirectionAtPoint(newVect._v, self.vectMat._m, self.vectMat._rows, point._v):
            return newVect

    def getNearestPoint(Wire self, Vector point):
        cdef Vector newVect = Vector(None)
        newVect._v = <double *> PyMem_Malloc (3 * sizeof(double))
        getNearestPoint(newVect._v, self.vectMat._m, self.vectMat._rows, point._v)
        return newVect

    def isClockWise(self, obsPoint:Vector):
        return isClockWise(self.vectMat._m, self.vectMat._rows, obsPoint._v)

    def getShifted(self, value):
        "Shift a wire by the given value (+) to inner (-) to outer"
        cdef Wire newWire = Wire(getShifted(self.vectMat._m, self.vectMat._rows, value))
        return newWire

    def getIntersectionPoint(self, other, incEdge=False):
        if isinstance(other, Wire):
            return getIntersectionPointWithWire(self.vectMat._m, self.vectMat._rows,
                                                (<Wire>other).vectMat._m, (<Wire>other).vectMat._rows,
                                                <bint>incEdge)
        elif isinstance(other, AlgLine.Line):
            return getIntersectionPointWithLine(self.vectMat._m, self.vectMat._rows,
                                                (<AlgLine.Line>other)._pnt, (<AlgLine.Line>other)._dir,
                                                <bint>incEdge)
        elif isinstance(other, AlgSegment.Segment):
            return getIntersectionPointWithSegment(self.vectMat._m, self.vectMat._rows,
                                                   (<AlgSegment.Segment>other)._vini, (<AlgSegment.Segment>other)._vfin,
                                                   <bint>incEdge)
        else:
            raise ValueError("Intersection point arguments are wrong")

    def simplify(self):
        "Simplifica el modelo"
        cdef unsigned int i, nexti, previ, nump
        if self.isClosed():
            nump = <unsigned int> (self.vectMat._rows - 1)
        else:
            nump = self.vectMat._rows
        i = 1
        while i < nump:
            previ = <unsigned int> ((i - 1) % nump)
            nexti = <unsigned int> ((i + 1) % nump)
            if AlgSegment.isInside(&self.vectMat._m[previ * 3], &self.vectMat._m[nexti * 3],
                                        &self.vectMat._m[i * 3], <bint>True):
                self.vectMat.deleteRow(i)
                nump -= 1
                continue
            i += 1

    def getSubWire(self, length1, length2):
        return getSubWire(self.vectMat._m,self.vectMat._rows, length1, length2)

    def copy(Wire self) -> Wire:
        return Wire(self.vectMat)

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memodict={}):
        return self.copy()

    def __repr__(self):
        cdef double * cdg = <double *> PyMem_Malloc ( 3 * sizeof(double))
        cdef double perim = getLength(self.vectMat._m, self.vectMat._rows)
        getCDG(cdg, self.vectMat._m, self.vectMat._rows)
        try:
            return "Wire(Perimeter=%.2f; CDG=(%.2f,%.2f,%.2f))" % (perim, cdg[0], cdg[1], cdg[2])
        finally:
            PyMem_Free(cdg)

    def __str__(self):
        cdef double * cdg = <double *> PyMem_Malloc ( 3 * sizeof(double))
        cdef double perim = getLength(self.vectMat._m, self.vectMat._rows)
        getCDG(cdg, self.vectMat._m, self.vectMat._rows)
        try:
            return "Wire(Perimeter=%.2f; CDG=(%.2f,%.2f,%.2f))" % (perim, cdg[0], cdg[1], cdg[2])
        finally:
            PyMem_Free(cdg)