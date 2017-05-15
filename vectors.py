from math import pi, cos, sin, tan, sqrt, acos, asin, atan

def radians(degrees):
    return pi*degrees/180

def degrees(radians):
    return 180*radians/pi
    
    
class Angle:
    
    def __init__(self, radians=0.0):
        self.radians = radians
        self.degrees = 180*radians/pi
        self.sin = sin(radians)
        self.cos = cos(radians)
        self.tan = self.sin/self.cos
        self.sec = 1/self.cos
        self.csc = 1/self.sin
        self.cot = self.cos/self.sin
        
    def __add__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            return Angle(self.radians + other)
        else:
            return Angle(self.radians + other.radians)
            
    def __radd__(self, other):
        return self + other
            
    def __sub__(self, other):
        if isinstance(other, int) or isinstance(other, float):
            return Angle(self.radians - other)
        else:
            return Angle(self.radians - other.radians)
            
    def __mul__(self, other):
        return Angle(self.radians*other)
        
    def __rmul__(self, other):
        return self * other
        
    def __truediv__(self, other):
        return Angle(self.radians/other)
        
    def __neg__(self):
        return Angle(self.radians * -1)
        
    def __repr__(self):
        return 'Angle({})'.format(self.radians)
        
    def __str__(self):
        return '{}\N{DEGREE SIGN}'.format(round(self.degrees, 2))
    

class Vector3:
    
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = x
        self.y = y
        self.z = z                                          
        
    def magnitude(self):
        return sqrt(self.x**2 + self.y**2 + self.z**2)
    
    def normalize(self):
        ux = self.x/self.magnitude()
        uy = self.y/self.magnitude()
        uz = self.z/self.magnitude()
        return Vector3(ux, uy, uz)
    
    def angles(self):
        alpha = Angle(acos(self.x/self.magnitude()))
        beta = Angle(acos(self.y/self.magnitude()))
        gamma = Angle(acos(self.z/self.magnitude()))
        return {'alpha': alpha, 'beta': beta, 'gamma': gamma}
        
    def scale(self, scalar):
        return Vector3(self.x*scalar, self.y*scalar, self.z*scalar)
    
    def __mul__(self, scalar):
        return self.scale(scalar)
        
    def __rmul__(self, scalar):
        return self.scale(scalar)
    
    def __truediv__(self, scalar):
        return self.scale(1/scalar)
        
    def __neg__(self):
        return self.scale(-1)
    
    def __repr__(self):
        return 'Vector3({}, {}, {})'.format(self.x, self.y, self.z)
    
    def __str__(self):
        output, plus, coords, hats = '', False, self.components(), 'ijk'
        def evaluate(index):
            component = ''
            sign = coords[index] > 0
            if round(abs(coords[index]), 12) > 0:
                component = hats[index]
                if round(abs(coords[index]), 12) != 1:
                    component = '{:.2f}{}'.format(abs(coords[index]), component)
                if not sign:
                    component = ' - ' + component
            return {'component': component, 'sign': sign}
                
        output, plus = evaluate(2)['component'], evaluate(2)['sign']
        for n in reversed(range(2)):
            current, previous = evaluate(n), evaluate(n + 1)
            if current['component'] != '':
                if plus:
                    output = ' + ' + output
                plus = current['sign']
            output = current['component'] + output 
        return output.strip()
    
    def translate(self, other):
        return Vector3(self.x + other.x, self.y + other.y, self.z + other.z)
    
    def __add__(self, other):
        return self.translate(other)
    
    def __sub__(self, other):
        return self.translate(-other)
    
    def dot(self, other):
        return self.x*other.x + self.y*other.y + self.z*other.z
    
    def cross(self, other):
        cx = self.y*other.z - self.z*other.y
        cy = self.z*other.x - self.x*other.z
        cz = self.x*other.y - self.y*other.x
        return Vector3(cx, cy, cz)
    
    def angle(self, other=None):
        if isinstance(other, Vector3):
            return Angle(acos(self.dot(other)/(self.magnitude()*other.magnitude())))
        elif isinstance(other, Line3):
            return Angle(acos(abs(self.dot(other.unit)/(self.magnitude()))))
        elif isinstance(other, Plane):
            return Angle((pi/2) - acos(abs(self.dot(other.normal)/(self.magnitude()))))
        else:
            return self.angles()
    
    def components(self, other=None):
        if other == None:
            return (self.x, self.y, self.z)
        else:
            unit = other.normalize()
            projection = unit.scale(self.dot(unit))
            return {'parallel': projection, 'perpendicular': self - projection}
    
    def project(self, other):
        if isinstance(other, Vector3):
            return self.components(other)['parallel']
        elif isinstance(other, Line3):
            return self.components(other.vector)['parallel'] + other.point
        elif isinstance(other, Plane):
            return self.components(other.normal)['perpendicular'] - other.d*other.normal
            
    def to_plane(self, point=None):
        if point == None:
            point = Vector3(0.0, 0.0, 0.0)
        unit = self.normalize()
        pa, pb, pc = unit.x, unit.y, unit.z
        pd = -unit.dot(point)
        return Plane(pa, pb, pc, pd)

    def to_line(self, point=None):
        if point == None:
            point = Vector3(0.0, 0.0, 0.0)
        return Line3(self, point) 
    
    def intersect(self, plane):
        scalar = -plane.d/self.dot(plane.normal)
        return self*scalar
    
    def distance(self, other):
        if isinstance(other, Vector3):
            return (self - other).magnitude()
        elif isinstance(other, Line3):
            return (self - self.project(other)).magnitude()
        elif isinstance(other, Plane):
            return (self - self.project(other)).magnitude()
       
 
class Line3:
    def __init__(self, a, b=None, c=None, x=None, y=None, z=None):
        if b == None:
            unit, point = a.normalize(), vector(0.0, 0.0, 0.0)
        elif c == None:
            unit, point = a.normalize(), b.components(a)['perpendicular']
        elif y == None:
            if isinstance(a, Vector3):
                unit = a.normalize()
                point = vector(b, c, x).components(unit)['perpendicular']
            else:
                unit = vector(a, b, c).normalize()
                point = x.components(unit)['perpendicular']
        else:
            unit = vector(a, b, c).normalize()
            point = vector(x, y, z).components(unit)['perpendicular']
        self.a = unit.x
        self.b = unit.y
        self.c = unit.z
        self.x = point.x
        self.y = point.y
        self.z = point.z
        self.unit = unit
        self.point = point

    def __repr__(self):
        return 'Line3({}, {}, {}, {}, {}, {})'.format(self.a, self.b, self.c, self.x, self.y, self.z)

    def __str__(self):
        return 'Line (Unit: {}, Point: {})'.format(self.unit, self.point)
    
    def translate(self, vector):
        return self.unit.to_line(self.point + vector)
        
    def __add__(self, vector):
        return self.translate(vector)
    
    def __radd__(self, vector):
        return self.translate(vector)
    
    def __sub__(self, other):
        return self.unit.to_plane(self.point - vector)
        
    def components(self):
        return (self.unit, self.point)
    
    def project(self, plane):
        return Line3(self.unit.components(plane.normal)['perpendicular'], self.point.project(plane))
    
    def intersect(self, plane):
        plane = plane - self.point
        scalar = -plane.d/self.unit.dot(plane.normal)
        return self.unit*scalar + self.point
    


class Plane:
    def __init__(self, a=0.0, b=0.0, c=0.0, d=0.0):
        if isinstance(a, Vector3):
            unit = a.normalize()
            if isinstance(b, Vector3):
                d = -unit.dot(b)
        else:
            unit = vector(a, b, c).normalize()
        self.a = unit.x
        self.b = unit.y
        self.c = unit.z 
        self.d = d
        self.normal = unit
        self.point = -d*unit

    def __repr__(self):
        return 'Plane({}, {}, {}, {})'.format(self.a, self.b, self.c, self.d)

    def __str__(self):
        return 'Plane (Normal: {}, Point: {})'.format(self.normal, self.point)
    
    def translate(self, vector):
        return self.normal.to_plane(self.point + vector)
        
    def __add__(self, vector):
        return self.translate(vector)
    
    def __radd__(self, vector):
        return self.translate(vector)
    
    def __sub__(self, vector):
        return self.translate(-vector)
    
    def components(self):
        return (self.a, self.b, self.c, self.d)
    
    def intersect(self, other):
        if isinstance(other, Plane):
            i_vector = self.normal.cross(other.normal)
            y_vector = i_vector.cross(self.normal).normalize()
            angle = self.normal.angle(other.normal).radians
            y_magnitude = (abs(other.d) - cos(angle)*abs(self.d))/sin(angle)
            i_point = self.point + y_magnitude*y_vector
            return Line3(i_vector, i_point)
        else:
            return other.intersect(self)
    

def angle(radians=0.0):
    return Angle(radians)
def unit_vector(x=0.0, y=0.0, z=0.0):
    return Vector3(x, y, z).normalize()
def vector(x=0.0, y=0.0, z=0.0):
    return Vector3(x, y, z)    
def line(a, b=None, c=None, x=None, y=None, z=None):
    return Line3(a, b, c, x, y, z)    
def plane(a=0.0, b=0.0, c=0.0, d=0.0):
    return Plane(a, b, c, d)

def vector_sum(*vectors):
    if hasattr(vectors[0], '__iter__'):
        vectors = vectors[0]
    resultant = Vector3(0, 0, 0)
    for vector in vectors:
        resultant += vector
    return resultant
    
def print_angles(vector):
    alpha = vector.angles()['alpha']
    beta = vector.angles()['beta']
    gamma = vector.angles()['gamma']
    print('Alpha: {}\nBeta: {}\nGamma: {}'.format(alpha, beta, gamma))

origin = vector(0, 0, 0)
i_hat, j_hat, k_hat = vector(1, 0, 0), vector(0, 1, 0), vector(0, 0, 1)
xy_plane, yz_plane, zx_plane = plane(0, 0, 1, 0), plane(1, 0, 0, 0), plane(0, 1, 0, 0)
