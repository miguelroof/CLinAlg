from cpython.mem cimport PyMem_Malloc, PyMem_Free
from libc.math cimport sin, cos, fabs
from .cimport AlgVector, AlgWire, AlgSurface, AlgLine, AlgVector2D
from .AlgTool cimport presition

cdef class NeutralFiber():
    "Class that contents Neutral fiber based on strain vector e0, psi_y, psi_z. Everithing based on 2D YZ plane"
    def __cinit__(self, *coord):
        self._sv = <double *> PyMem_Malloc(3 * sizeof(double))

    def __init__(self, *coord):
        if len(coord) == 3:
            self._sv[0] = coord[0]
            self._sv[1] = coord[1]
            self._sv[2] = coord[2]
        elif len(coord) == 1:
            self._sv[0] = coord[0][0]
            self._sv[1] = coord[0][1]
            self._sv[2] = coord[0][2]

    def __dealloc__(self):
        if self._sv != NULL:
            PyMem_Free(self._sv)
            self._sv = NULL

    @property
    def strain_vector(self):
        cdef AlgVector.PhantomVector sv = AlgVector.PhantomVector()
        sv._v = self._sv
        sv._ref_object = self
        return sv

    @property
    def curvature(self):
        return AlgVector2D.module(&self._sv[1])

    def z(self, point):
        "Get the distance from one point to the fiber neutre"
        cdef double curvature = self.curvature
        if curvature < 1e-19:
            return 0.0  # no hay fibra neutra. a infinito
        if isinstance(point, (int, float)):
            # point is a deformation e
            return (point - self._sv[0]) / curvature
        else:  # asumire un punto en 3d
            return (self.e(point) / curvature)

    def e(self, point):
        cdef double curvature = self.curvature
        if curvature < 1e-19:
            return self._sv[0]
        if isinstance(point, (int, float)):
            return self._sv[0] + point * curvature
        else:
            return self._sv[0] - self._sv[2] * point[1] + self._sv[1] * point[2]

    def get_far_point(self, values, to_positif_z=True):
        # busca el punto mas lejano de entre los pasados
        cdef double curvature = self.curvature
        if isinstance(values, AlgSurface.Surface):
            if curvature < 1e-19:
                return (<AlgSurface.Surface> values).Point[0], 0.0
            dist_list = [self.e(pv) / curvature for pv in (<AlgSurface.Surface> values).Point]
            my_dist = max(dist_list) if to_positif_z else min(dist_list)
            pv = values.Point[dist_list.index(my_dist)]
            return pv, my_dist
        elif isinstance(values, AlgWire.Wire):
            if curvature < 1e-19:
                return (<AlgWire.Wire> values).Point[0], 0.0
            dist_list = [self.e(pv) / curvature for pv in (<AlgSurface.Surface> values).Point]
            my_dist = max(dist_list) if to_positif_z else min(dist_list)
            pv = values.Point[dist_list.index(my_dist)]
            return pv, my_dist
        elif isinstance(values, (list, tuple)):
            if curvature < 1e-19:
                return values[0], 0.0
            dist_list = [self.e(pv) / curvature for pv in values.Point]
            my_dist = max(dist_list) if to_positif_z else min(dist_list)
            pv = values[dist_list.index(my_dist)]
            return pv, my_dist

    def get_line(self, z: float = 0.0):
        "Returns a line at z distnace"
        # cdef double curvature = self.curvature
        cdef int signo
        if self._sv[1] == 0.0 and self._sv[2] == 0.0:
            # tengo que devolver una linea en la cachimbamba
            signo = -1 if self._sv[0] >= 0 else 1
            return AlgLine.Line(AlgVector.Vector(0, 0, 1e19 * signo), AlgVector.Vector(0, 1, 1e19 * signo))
        if self._sv[1] != 0:
            p0 = AlgVector.Vector(0, 0, -self._sv[0] / self._sv[1])
            p1 = AlgVector.Vector(0, self._sv[1], self._sv[2])
        else:
            p0 = AlgVector.Vector(0, self._sv[0] / self._sv[2], 0)
            p1 = AlgVector.Vector(0, self._sv[1], self._sv[2])
        p1.normalize()
        return AlgLine.Line(p0, p0 + p1)

    def copy(self):
        return NeutralFiber(self.strain_vector)

    # la algebra
    def translate_by_e(self, delta_e: float):
        self._sv[0] = self._sv[0] + delta_e

    def translate_by_z(self, z):
        self.translate_by_e(-self.e(z))

    def roll(self, center, angle):
        "Funcion que rota en torno al eje x"
        center_strain = self.e(center)
        psi_y = self._sv[1] * cos(angle) - self._sv[2] * sin(angle)
        psi_z = self._sv[1] * sin(angle) + self._sv[2] * cos(angle)
        e0 = center_strain - center[2] * psi_y + center[1] * psi_z
        self._sv[0] = e0
        self._sv[1] = psi_y
        self._sv[2] = psi_z

    def pitch(self, pivot_point, ref_point, ref_e):
        "Funcion que gira la deformada, manteniendo la inclinacion y el valor en pivot_point, para que en ref_point tenga ref_value"
        cdef double zx, zc, ex, ec, curve, ncurve
        cdef AlgVector.Vector direction
        curve = self.curvature
        if curve:  # tenia curvatura previa
            zx = self.z(ref_point)
            ex = self.e(ref_point)
            zc = self.z(pivot_point)
            ec = self.e(pivot_point)
            z0 = self._sv[0] / curve
            if fabs(zx - zc) < presition:
                return
            ncurve = (ref_e - ec) / (zx - zc)
            self._sv[0] = ref_e - ncurve * (zx - z0)
            self._sv[1] = self._sv[1] * ncurve / curve
            self._sv[2] = self._sv[2] * ncurve / curve
        else:
            ec = self._sv[0]
            ncurve = (ref_e - ec) / pivot_point.distance(ref_point)
            direction = (pivot_point - ref_point).normalize()
            if ref_e < ec:
                direction = - direction
            self._sv[2] = ncurve * direction[1]
            self._sv[1] = -ncurve * direction[2]
            self._sv[0] = ec + self._sv[2] * pivot_point[1] - self._sv[1] * pivot_point[2]

    def set_curvature(self, point, strain):
        "Funcion que ajusta la curvatura para que cumpla con las condiciones"
        z_point = self.z(point)
        if z_point < presition:
            raise RuntimeError("POint for set_curvature should be far from NF line")
        old_curv = self.curvature
        if old_curv < presition:  # antes no habia curvatura, por lo que no tengo referencia
            self._sv[0] = strain
            return
        new_curv = strain / z_point
        self._sv[0] = self._sv[0] * new_curv / old_curv
        self._sv[1] = self._sv[1] * new_curv / old_curv
        self._sv[2] = self._sv[2] * new_curv / old_curv

    def __eq__(NeutralFiber self, NeutralFiber other) -> bool:
        cdef int k
        return all([fabs(self._sv[k] - other._sv[k]) < presition for k in range(3)])

    def __ne__(NeutralFiber self, NeutralFiber other):
        return not self.__eq__(other)

    def __add__(NeutralFiber self, NeutralFiber other):
        return NeutralFiber(self.strain_vector + other.strain_vector)

    def __sub__(NeutralFiber self, NeutralFiber other):
        return NeutralFiber(self.strain_vector - other.strain_vector)

    def __mul__(NeutralFiber self, other):
        return NeutralFiber(other * self.strain_vector)

    def __rmul__(NeutralFiber self, other):
        return NeutralFiber(other * self.strain_vector)

    def __bool__(NeutralFiber self):
        cdef int k
        return any([fabs(self._sv[k]) > presition for k in range(3)])

    def __serialize__(NeutralFiber self) -> list:
        cdef unsigned int i
        return [self._sv[i] for i in range(3)]

    def toList(NeutralFiber self) -> list:
        return self.__serialize__()

    def pythonized(NeutralFiber self) -> list:
        return self.__serialize__()

    def toTuple(NeutralFiber self) -> tuple:
        return (self._sv[0], self._sv[1], self._sv[2])

    def __hash__(NeutralFiber self):
        return hash(self.toTuple())

    def __repr__(NeutralFiber self) -> str:
        return 'NF(%.2f,%.2f,%.2f)' % (self._sv[0], self._sv[1], self._sv[2])

    def __str__(NeutralFiber self) -> str:
        return 'NF(%.2f,%.2f,%.2f)' % (self._sv[0], self._sv[1], self._sv[2])

    def __getitem__(NeutralFiber self, unsigned int key) -> float:
        return self._sv[key]

    def __setitem__(NeutralFiber self, unsigned int key, float value):
        self._sv[key] = value

    def __copy__(NeutralFiber self):
        return self.copy()

    def __deepcopy__(NeutralFiber self, memo=None):
        return self.copy()


    def __json__(self):
        return {'__jsoncls__': 'CLinAlg.ToolFiber:NeutralFiber.from_JSON', 'vector': self.toTuple()}

    @classmethod
    def from_JSON(cls, jsondict):
        cdef NeutralFiber nf = cls(jsondict['vector'])
        return nf