import sys
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from libc.math cimport M_PI, fabs
from .AlgTool cimport presition
from . cimport AlgVector, AlgWire, AlgSegment, AlgMatrix, AlgLine, AlgQuaternion


cdef list tessellate(double * surfPoint, unsigned int numpoints):
    """
    Tessellate a surface into triangles
    :param surfPoint: pointer to array that keeps points in rows
    :param numpoints: num of points in the array
    :return: list of index of points for triangles.
    """
    # toda la algebra de la funcion raiz se definia a partir de un contorno CERRADO
    cdef int i, j, ci, pi, ni
    cdef double modulo
    cdef double * ndir = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * pdir = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * axis = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * vtemp = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef list triangles = []
    try:
        indices = list(range(numpoints - 1))
        AlgWire.getNormal(axis, surfPoint, numpoints)
        i = 0
        while len(indices) > 3 and i < len(indices):
            pi = indices[-1] if i == 0 else indices[i - 1]
            ci = indices[i]
            ni = indices[i + 1] if i < len(indices) - 1 else indices[0]
            AlgVector.vdir(pdir, &surfPoint[pi * 3], &surfPoint[ci * 3])
            AlgVector.vdir(ndir, &surfPoint[ci * 3], &surfPoint[ni * 3])
            AlgVector.cross(vtemp, pdir, ndir)
            modulo = AlgVector.module(vtemp)
            if modulo > presition:
                for j in range(3):
                    vtemp[j] /= modulo
            else:
                i = i + 1
                continue
            AlgVector.add(vtemp, vtemp, axis)
            modulo = AlgVector.module(vtemp)
            if modulo < 0.5:
                i = i + 1
                continue
            # ahora detecto colisiones
            AlgVector.vdir(vtemp, &surfPoint[pi * 3], &surfPoint[ni * 3])
            linecortes = AlgWire.getIntersectionPointWithLineByIndex(surfPoint, indices[:] + [indices[0]], &surfPoint[ci * 3], vtemp, <bint>False)
            cortes = []
            valid = True
            for j in range(len(linecortes)):
                if AlgSegment.isInside(&surfPoint[pi * 3], &surfPoint[ni * 3],
                                            (<AlgVector.Vector>(linecortes[i]))._v, <bint>False):
                    # se ha producido un corte con otro tramo del wire. No lo doy por falido
                    valid = False
                    break
            if not valid:
                i = i + 1
                continue
            triangles.extend([pi, ci, ni])
            indices.pop(i)
        triangles.extend([indices[0], indices[1], indices[2]])
        return triangles
    finally:
        PyMem_Free(ndir)
        PyMem_Free(pdir)
        PyMem_Free(axis)
        PyMem_Free(vtemp)

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
        for i in range(0, len(surf.triMat), 3):
            AlgVector.sub(vAB, &surf.vectMat._m[surf.triMat._m[i + 1] * 3], &surf.vectMat._m[surf.triMat._m[i] * 3])
            AlgVector.sub(vAC, &surf.vectMat._m[surf.triMat._m[i + 2] * 3], &surf.vectMat._m[surf.triMat._m[i] * 3])
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
        for i in range(0, len(surf.triMat), 3):
            ox = sum([surf.vectMat._m[surf.triMat._m[i + k] * 3 + 0] for k in range(3)]) / 3.0
            oy = sum([surf.vectMat._m[surf.triMat._m[i + k] * 3 + 1] for k in range(3)]) / 3.0
            oz = sum([surf.vectMat._m[surf.triMat._m[i + k] * 3 + 2] for k in range(3)]) / 3.0
            AlgVector.sub(vAB, &surf.vectMat._m[surf.triMat._m[i + 1] * 3], &surf.vectMat._m[surf.triMat._m[i] * 3])
            AlgVector.sub(vAC, &surf.vectMat._m[surf.triMat._m[i + 2] * 3], &surf.vectMat._m[surf.triMat._m[i] * 3])
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
    cdef AlgMatrix.Matrix inertia = AlgMatrix.zeros(3,3)
    cdef AlgMatrix.Matrix V = AlgMatrix.zeros(3, 3)
    cdef AlgMatrix.Matrix S = AlgMatrix.Matrix(
        [[2.0 / 24, 1.0 / 24, 1.0 / 24], [1.0 / 24, 2.0 / 24, 1.0 / 24], [1.0 / 24, 1.0 / 24, 2.0 / 24]])
    cdef AlgMatrix.Matrix J = AlgMatrix.zeros(3, 3)
    cdef AlgMatrix.Matrix IDD = AlgMatrix.identity(3)
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
        for i in range(0, len(surf.triMat), 3):
            v0 = &surf.vectMat._m[surf.triMat._m[i] * 3]
            v1 = &surf.vectMat._m[surf.triMat._m[i + 1] * 3]
            v2 = &surf.vectMat._m[surf.triMat._m[i + 2] * 3]
            for j in range(3):
                for k in range(3):
                    V._m[j * 3 + k] = surf.vectMat._m[surf.triMat._m[i + j] * 3 + k]
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
            J = alfa * IDD
            J = J - C
            inertia = inertia + J
        inertia = stenierInGlobalAxis_p(inertia, surf.Area, pzero, (<AlgVector.Vector> surf.CDG)._v)
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
        for i in range(0, len(surf.triMat), 3):
            for j in range(3):
                cdg[j] = sum([surf.vectMat._m[surf.triMat._m[i + k] * 3 + j] for k in range(3)]) / 3.0
            AlgVector.sub(vAB, &surf.vectMat._m[surf.triMat._m[i + 1] * 3], &surf.vectMat._m[surf.triMat._m[i] * 3])
            AlgVector.sub(vAC, &surf.vectMat._m[surf.triMat._m[i + 2] * 3], &surf.vectMat._m[surf.triMat._m[i] * 3])
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
        for i in range(0, len(surf.triMat), 3):
            for j in range(3):
                AlgVector.sub(&mm[j * 3], &surf.vectMat._m[surf.triMat._m[i + j] * 3], point)
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
                            if abs(surf.triMat._m[i + indices[j][0]] - surf.triMat._m[
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
    cdef double * ipoint
    cdef double * tpoint = <double *> PyMem_Malloc (3 * sizeof(double))
    cdef bint sclosed
    cdef double * dirMIn = <double *> PyMem_Malloc(3 * sizeof(double))
    cdef double * dirMOut = <double *> PyMem_Malloc(3 * sizeof(double))
    wireInd = list(range(surf.vectMat._rows))
    n = <int> (len(wireInd))
    sclosed = True
    try:
        cortes = []
        indices = []
        for i in range(n - 1):
            pc = i
            pn = i+1
            ipoint = <double *> PyMem_Malloc (3 * sizeof(double))
            if not AlgSegment.getIntersectionPointWithLine(ipoint, &surf.vectMat._m[pc * 3],
                                                                  &surf.vectMat._m[pn * 3],
                                                                  point,
                                                                  dir,
                                                                  incEdge) or any([AlgVector.isEqual(ipoint, (<AlgVector.Vector> v)._v) for v in cortes]):
                PyMem_Free(ipoint)
                continue
            if (incEdge or (not AlgVector.isEqual(ipoint, &surf.vectMat._m[pc * 3])
                            and not AlgVector.isEqual(ipoint, &surf.vectMat._m[pn * 3]))):
                newVect = AlgVector.Vector(None)
                (<AlgVector.Vector>newVect)._v = ipoint
                cortes.append(newVect)
                indices.append(i)
                continue
            if AlgVector.isEqual(ipoint, &surf.vectMat._m[pc * 3]):  # punto inicial
                if i == 0:
                    if not sclosed:
                        PyMem_Free(ipoint)
                        continue
                    pp = n-2
                    AlgVector.vdir(dirMIn, &surf.vectMat._m[pc * 3], &surf.vectMat._m[pp * 3])
                else:
                    pp = <int>(i-1)
                    AlgVector.vdir(dirMIn, &surf.vectMat._m[pc * 3], &surf.vectMat._m[pp * 3])
                AlgVector.vdir(dirMOut, &surf.vectMat._m[pc * 3], &surf.vectMat._m[pn * 3])
                AlgVector.add(dirMIn, &surf.vectMat._m[pc * 3], dirMIn)
                AlgVector.add(dirMOut, &surf.vectMat._m[pc * 3], dirMOut)
            elif AlgVector.isEqual(ipoint, &surf.vectMat._m[pn * 3]):  #punto final
                if i == n - 2:
                    if not sclosed:
                        PyMem_Free(ipoint)
                        continue
                    pp = 1
                    AlgVector.vdir(dirMOut, &surf.vectMat._m[pn * 3], &surf.vectMat._m[pp * 3])
                else:
                    pp = i+2
                    AlgVector.vdir(dirMOut, &surf.vectMat._m[pn * 3], &surf.vectMat._m[pp * 3])
                AlgVector.vdir(dirMIn, &surf.vectMat._m[pn * 3], &surf.vectMat._m[pc * 3])
                AlgVector.add(dirMIn, &surf.vectMat._m[pn * 3], dirMIn)
                AlgVector.add(dirMOut, &surf.vectMat._m[pn * 3], dirMOut)
            else:  #ninguna de las anteriores
                PyMem_Free(ipoint)
                raise ValueError('No es posible esta combinacion')
            if AlgSegment.getIntersectionPointWithLine(tpoint, dirMIn, dirMOut, point,
                                                                  dir, <bint>False):
                newVect = AlgVector.Vector(None)
                (<AlgVector.Vector>newVect)._v = ipoint
                cortes.append(newVect)
                indices.append(i)
            else:
                PyMem_Free(ipoint)
        return indices, cortes
    finally:
        PyMem_Free(dirMIn)
        PyMem_Free(dirMOut)
        PyMem_Free(tpoint)

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
    for i in range(surf.vectMat._rows-1):
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

cdef class MatrixTriangle():
    "Objeto simple para representar un triangulo de una superficie"
    def __cinit__(MatrixTriangle self, data=None):
        if data:
            self.setList(data)

    def __dealloc__(self):
        PyMem_Free(self._m)

    cpdef void setList(MatrixTriangle self, list data):
        cdef int i
        if self._m:
            PyMem_Free(self._m)
        self._m = <int *> PyMem_Malloc(len(data) * sizeof(int))
        self._rows = <unsigned int> (len(data) // 3)
        # print("numero de filas", self._rows)
        for i in range(len(data)):
            self._m[i] = <int> (data[i])

    def copy(self) -> MatrixTriangle:
        return MatrixTriangle(self.__serialize__())

    def __serialize__(self) -> list:
        cdef int i
        return [self._m[i] for i in range(self._rows * 3)]

    def __len__(MatrixTriangle self):
        return self._rows * 3

    def __json__(self):
        return {'__jsoncls__': 'CLinAlg.AlgSurface:MatrixTriangle.from_JSON', 'data': self.__serialize__()}

    @classmethod
    def from_JSON(cls, jsondict):
        obj = cls(jsondict['data'])
        return obj
    

    
    
cdef class Surface():
    
    def __cinit__(self, pointList=[]):
        cdef unsigned int i,j
        self.vectMat = Matrix()
        self.triMat = MatrixTriangle()
        self.surfaceProp = {}
        if pointList:
            if isinstance(pointList, (list, tuple)):
                if pointList[0] != pointList[-1]:
                    pointList.append(pointList[0])
                self.vectMat._m = <double *> PyMem_Malloc(len(pointList) * 3 * sizeof(double))
                self.vectMat._rows = <unsigned int> len(pointList)
                self.vectMat._cols = 3
                for i in range(len(pointList)):
                    for j in range(3):
                        self.vectMat._m[i * 3 + j] = pointList[i][j]
            elif isinstance(pointList, Matrix):
                self.vectMat._m = <double *> PyMem_Malloc((<Matrix> pointList)._rows * (<Matrix> pointList)._cols * sizeof(double))
                self.vectMat._rows = (<Matrix> pointList)._rows
                self.vectMat._cols = (<Matrix> pointList)._cols
                for i in range((<Matrix> pointList)._rows * (<Matrix> pointList)._cols):
                    self.vectMat._m[i] = (<Matrix> pointList)._m[i]
            triangles = tessellate(self.vectMat._m, self.vectMat._rows)
            self.triMat.setList(triangles)

    def __json__(self):
        return {'__jsoncls__': 'CLinAlg:Surface.from_JSON', 'mat': self.vectMat, 'triMat': self.triMat,
                'surfaceProp': self.surfaceProp}

    @classmethod
    def from_JSON(cls, jsondict):
        obj = cls()
        obj.vectMat = jsondict['mat']
        obj.triMat = jsondict['triMat']
        obj.surfaceProp = jsondict['surfaceProp']
        return obj
        

    @property
    def Contour(self):
        """
        Return wire of contour
        """
        cdef AlgWire.Wire wire = AlgWire.Wire()
        wire.vectMat = self.vectMat
        return wire

    @property
    def Edges(self):
        "Return edges as iterator"
        return AlgWire.SegmentIterator(self.vectMat)

    @property
    def Points(self):
        "Return points as iterator"
        return AlgWire.PointIterator(self.vectMat)

    @property
    def Area(self):
        if not self.surfaceProp.get('area'):
            self.surfaceProp['area'] = getArea(self)
        return self.surfaceProp['area']

    @property
    def CDG(self):
        if not self.surfaceProp.get('cdg'):
            cdg = AlgVector.Vector(None)
            (<AlgVector.Vector>cdg)._v = <double *> PyMem_Malloc (3 * sizeof(double))
            getCDG((<AlgVector.Vector>cdg)._v, self)
            self.surfaceProp['cdg'] = cdg
        return self.surfaceProp['cdg']

    @property
    def Inertia(self):
        if not self.surfaceProp.get('inertia'):
            self.surfaceProp['inertia'] = getInertia(self)
        return self.surfaceProp['inertia']

    def staticMoment(self, line:AlgLine.Line):
        "Returns static moment to a line (point and direction)"
        return staticMoment(self, line._pnt, line._dir)

    def modulus(self):
        "Return section modulus at main axis of inertia"
        cdef list inertia
        cdef Matrix axis
        inertia, axis = mainAxisInertia(self.Inertia)  #lista de autovalores y matriz de autovectores
        cdg = <AlgVector.Vector>self.CDG
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
        inertia, axis = mainAxisInertia(self.Inertia)
        area = self.Area
        radius = []
        for i in range(3):
            radius.append((inertia[i] / area) ** 0.5)
        return radius, axis

    def isInside(self, point:AlgVector.Vector, incEdge=True):
        "Funcion que calcula si un punto esta dentro de la superficie"
        return isInside(self, point._v, <bint>incEdge)

    def getRotate(self, center:AlgVector.Vector, axis:AlgVector.Vector, angle:float):
        cdef double * quat = <double *> PyMem_Malloc (4 * sizeof(double))
        cdef Surface newSurface = Surface()
        newSurface.vectMat._m = <double *> PyMem_Malloc (3 * self.vectMat._rows*sizeof(double))
        newSurface.vectMat._rows = self.vectMat._rows
        newSurface.vectMat._cols = 3
        AlgQuaternion.quatFromAxisAngle(quat, axis._v, angle)
        AlgQuaternion.rotateVectorMatrix(newSurface.vectMat._m, quat, center._v, self.vectMat._m, self.vectMat._rows)
        newSurface.triMat.setList(tessellate(newSurface.vectMat._m, newSurface.vectMat._rows))
        PyMem_Free(quat)
        return newSurface

    def getTranslate(self, translateVector:AlgVector.Vector):
        """
        Returns a copy of surface translate by the given vector
        :param translateVector: Vector with the translation
        """
        cdef Surface newSurface = Surface()
        cdef int i, j
        newSurface.vectMat._m = <double *> PyMem_Malloc (self.vectMat._rows*3*sizeof(double))
        newSurface.vectMat._rows = self.vectMat._rows
        newSurface.vectMat._cols = self.vectMat._cols
        for i in range(<int>self.vectMat._rows):
            for j in range(3):
                newSurface.vectMat._m[j+i*3] = self.vectMat._m[j+i*3]+translateVector._v[j]
        newSurface.triMat._m = <int *> PyMem_Malloc (3*self.triMat._rows*3*sizeof(int))
        newSurface.triMat._rows = self.triMat._rows
        for i in range(3*self.triMat._rows):
            newSurface.triMat._m[i] = self.triMat._m[i]


    def splitByLine(self, line:AlgLine.Line):
        """
        Return a list with news surfaces created by cutting self surface by the line
        :param line: AlgLine.Line object
        :return: list of surfaces
        """
        print("SPLITBYLINE of SURFACE: Esta funcion todavia no esta completamente depurada")
        indcortes, puntocortes = dictIntersectByLine(self, line._pnt, line._dir, <bint>True)
        if (not indcortes) or (len(indcortes) <= 1):
            return [self.copy()]
        numvert = self.vectMat._rows-1
        dictcortes = {indcortes[i]: puntocortes[i] for i in range(len(indcortes))}
        sortcortes = sorted(dictcortes.keys(), key=lambda v: line.getParameter(v))
        sortedges = [dictcortes[v] for v in sortcortes]
        sortindex = [self.Vertexes.index(seg.pfin) for seg in sortedges]
        i = sortindex[0]
        newSect = [sortcortes[0]]
        vertexes = self.Vertexes
        if sortcortes[0] != vertexes[i % numvert]:
            newSect.append(vertexes[i % numvert])
        i = i + 1
        while (i % numvert) not in sortindex:
            newSect.append(vertexes[i % numvert])
            i += 1
        index = sortindex.index(i % numvert)
        if index > 1:  # hay que cambiar el giro
            newSect = [sortcortes[0]]
            i = i - 1
            if sortcortes[0] != vertexes[i % numvert]:
                newSect.append(vertexes[i % numvert])
            i = i - 1
            while (i % numvert) not in sortindex:
                newSect.append(vertexes[i % numvert])
                i = i - 1
        newSect.append(sortcortes[1])
        oldSect = []
        for v in range(sortindex[0], sortindex[0] + numvert):
            if vertexes[v % numvert] not in newSect:
                oldSect.append(vertexes[v % numvert])
        oldSect.append(sortcortes[0])
        oldSect.append(sortcortes[1])
        newSurf = [Surface(newSect)]
        oldSurf = Surface(oldSect)
        newSurf.extend(oldSurf.splitByLine(line))
        return newSurf

    def intersectByLine(self, line:AlgLine.Line, incEdge=True):
        "Returns a list of intersection points with a line"
        puntos = list(range(self.vectMat._rows))
        cortes = AlgWire.getIntersectionPointWithLine(self.vectMat._m, self.vectMat._rows, line._pnt, line._dir, <bint>incEdge)
        if len(cortes) < 2:
            return []
        return cortes

    def intersects(self, other:Surface) -> bool:
        "Returns True if surfaces intersect with each other"
        p: AlgVector.Vector
        for p in other.Points:
            if isInside(self, p._v, <bint>True):
                return True
        for p in self.Points:
            if isInside(other, p._v, <bint>True):
                return True
        intPoints = AlgWire.getIntersectionPointWithWire(self.vectMat._m, self.vectMat._rows, other.vectMat._m, other.vectMat._rows, <bint>True)
        return intPoints != []

    def intersected(self, surface:Surface):
        "Returns a list of surfaces that intersect with this surface"
        print("INTERSECTED of Surface: Esta funcion todavia no esta completamente depurada")
        myWire = self.Contour
        otherWire = surface.Contour
        intPoints = AlgWire.getIntersectionPointWithWire(self.vectMat._m, self.vectMat._rows, surface.vectMat._m, surface.vectMat._rows, <bint>True)
        if not intPoints:
            if self.isInside(surface.Points[0]):
                newSurface = Surface(self.vectMat)
                return newSurface
            elif surface.isInside(self.Points[0]):
                newSurface = Surface(surface.vectMat)
                return newSurface
            else:
                return None
        if len(intPoints) <= 1:
            return None
        fromMine = True
        intLenghts = [AlgWire.getLengthFromOrigin(self.vectMat._m, self.vectMat._rows, (<AlgVector.Vector>x)._v) for x in intPoints]
        if surface.isInside(self.Points[0], True):
            startInside = True
            subWire = AlgWire.getSubWire(self.vectMat._m, self.vectMat._rows, 0, intLenghts[0])
            index = 0
        else:
            startInside = False
            subWire = AlgWire.getSubWire(self.vectMat._m, self.vectMat._rows, intLenghts[0], intLenghts[1])
            index = 1
        points = subWire.points[:]
        while index < len(intPoints) - 1:
            if fromMine:
                subWire = otherWire.getSubWire(intPoints[index], intPoints[index + 1])
            else:
                subWire = myWire.getSubWire(intPoints[index], intPoints[index + 1])
            points.extend(subWire.points[1:])
            index += 1
            fromMine = not fromMine
        if startInside:
            subWire = myWire.getSubWire(intPoints[-1], 0)
        else:
            subWire = otherWire.getSubWire(intPoints[-1], intPoints[0])
        points.extend(subWire.points[1:-1])
        return Surface(points)

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
        newSurf.triMat = self.triMat.copy()
        for k, v in self.surfaceProp.items():
            if isinstance(v, (AlgMatrix.Matrix, AlgVector.Vector)):
                newSurf.surfaceProp[k] = v.copy()
            else:
                newSurf.surfaceProp[k] = v
        return newSurf

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memodict=None):
        return self.copy()

    def __repr__(self):
        return "Surface(Area=%.2f; Perimeter=%.2f; CDG=(%.2f,%.2f,%.2f))" % (
            self.Area, self.Contour.Perimeter, self.CDG.x, self.CDG.y, self.CDG.z)

    def __str__(self):
        return "Surface(Area=%.2f; Perimeter=%.2f; CDG=(%.2f,%.2f,%.2f))" % (
            self.Area, self.Contour.Perimeter, self.CDG.x, self.CDG.y, self.CDG.z)






    
    