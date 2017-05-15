#from vectors import *


class Force3:
    def __init__(self, magnitude, direction=None, lever_arm=None):
        if isinstance(magnitude, Vector3):
            if isinstance(direction, Vector3):
                lever_arm = direction
            magnitude, direction = magnitude.magnitude(), magnitude
        if direction.magnitude() > 0:
            vector = abs(magnitude)*direction.normalize()
        else: vector = Vector3()

        if lever_arm == None:
            lever_arm = Vector3(0, 0, 0)
        self.magnitude = abs(magnitude)
        self.vector = vector
        self.lever_arm = lever_arm
        self.moment = lever_arm.cross(vector)
        
    def __repr__(self):
        return 'Force3({}, {}, {})'.format(self.magnitude, self.vector, self.lever_arm)
    def __str__(self):
        return 'Force (Magnitude: {}, Vector: {}, Lever Arm: {}, Moment: {})'.format(self.magnitude, 
                                                                                     self.vector, 
                                                                                     self.lever_arm, 
                                                                                     self.moment)
    
    def __mul__(self, scalar):
        return Force3(self.magnitude*abs(scalar), self.vector*scalar, self.lever_arm)
    
    def __rmul__(self, scalar):
        return self*scalar
    
    def __truediv__(self, scalar):
        return self*(1/scalar)
    
    def __neg__(self):
        return self*(-1)
    
    def translate(self, other):
        return Force(self.magnitude, self.vector, self.lever_arm + other)
    
    def __add__(self, other):
        if isinstance(other, Vector3):
            return self.translate(other)
        else:
            vector = self.vector + other.vector
            magnitude = vector.magnitude()
            moment = self.moment + other.moment
            arms = [0, 0, 0]
            for axis in range(3):
                coords = [0, 0, 0]; coords[axis] = 1; unit = Vector3(*coords)
                sum_vector = vector.components(unit)['perpendicular']
                if sum_vector.magnitude() > 0:
                    self_vector, other_vector = self.vector.project(sum_vector), other.vector.project(sum_vector)
                else:
                    self_vector, other_vector = Vector3(), Vector3()
                self_arm, other_arm = self.lever_arm.project(unit), other.lever_arm.project(unit)
                self_moment, other_moment = Vector3(), Vector3()
                if self_vector.magnitude()*self_arm.magnitude() > 0:
                    self_moment = self_arm.cross(self_vector)
                if other_vector.magnitude()*other_arm.magnitude() > 0:
                    other_moment = other_arm.cross(other_vector)
                try:
                    
                    arm_unit = sum_vector.cross(self_moment + other_moment).normalize()
                    arms[axis] = (arm_unit*round(((self_moment + other_moment).magnitude()/sum_vector.magnitude()), 12)).components()[axis]
                except:
                    continue
                
            lever_arm = Vector3(*arms)
            return Force3(magnitude, vector, lever_arm)
        
    def __radd__(self, other):
        return self + other
    
    def __sub__(self, other):
        return self + Force3(other.magnitude, -other.vector, other.lever_arm)
        
    def project(other):
        vector = self.vector.project(other) - other.point
        magnitude = self.vector.magnitude()
        lever_arm = self.lever_arm.project(other)
        return Force3(magnitude, vector, lever_arm)
        
    
    def add_moment(self, moment, constraint=None):
        if constraint == None:
            vector = self.vector
        else:
            vector = self.vector.components(constraint)['perpendicular']
        component = vector.components(moment)['perpendicular']
        direction = component.cross(moment).normalize()
        lever_arm = (moment.magnitude()/component.magnitude())*direction
        return Force3(self.magnitude, self.vector, self.lever_arm + lever_arm)

def force(magnitude, vector, lever_arm):
    return Force3(magnitude, vector, lever_arm)
    
def force_sum(*forces):
    if hasattr(forces[0], '__iter__'):
        forces = forces[0]
    resultant = forces[0]
    for f in range(1, len(forces)):
        resultant += forces[f]
    return resultant
