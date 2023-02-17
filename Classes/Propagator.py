from math import sqrt, sin, sinh, cos, cosh
from numpy import cross, dot, array

class Propagator:

    #   alp     : alpha = 1/a
    #   dt      : time delta
    #   epoch   : Starting time
    #   SFC      : Stumpff Function output C
    #   SFS      : Stumpff Function output S
    #   smja    : Semi-Major Axis (a)
    #   X       : Universal Anomaly
    #   z       : "z" = alpha * X ** 2

    def __init__(self, sys=None, smja=0, timedelta=0.001):
        if sys is None:
            print("No system initialized")
        else:
            self.SFC = 0
            self.SFS = 0

            self.dt = timedelta
            self.thres = 10**-8
            self.sys = sys

            self.propagate_system()

    def propagate_system(self):
        for body in self.sys.Bodies.values():
            r_mag = sqrt(dot(body.r, body.r))
            v_mag = sqrt(dot(body.v, body.v))

            if body.smja is None:
                body.alp = 2/r_mag - (v_mag**2)/self.sys.mu

            self.find_universal_anomaly(body)
            body.r, body.v = self.solve_pos_vel(r_mag, v_mag, body)

    def solve_pos_vel(self, r_0_mag, v_0_mag, body):
        """
        :variables
            r_0     :   Initial Position
            r_0_mag :   Initial position magnitude
            v_0     :   Initial velocity vector
            v_0_mag :   Initial velocity magnitude
            f       :   Lagrange multiplier f
            g       :   Lagrange multiplier g
            f_dot   :   Lagrange multiplier f_dot
            g_got   :   Lagrange multiplier f_dot

        :param body:
        :return: f, g, g_dot, f_dot
        """

        r_0 = body.r
        v_0 = body.v

        v_r = body.v_r

        #   Calculate f and g
        f = 1 - body.X**2/r_0_mag*self.SFC
        g = self.dt - 1/sqrt(self.sys.mu)*body.X**3 *self.SFS

        #   Calculate r and r_mg
        r = f*r_0 + g*v_0
        r_mag = sqrt(dot(r, r))

        #   Calculate f_dot and g_dot
        g_dot = 1 - body.X**2/r_mag*self.SFC
        f_dot = sqrt(self.sys.mu) / (r_0_mag * r_mag) * (
                    body.alp * body.X ** 3 * self.SFS - body.X)
        # Calculate v
        v = f_dot*r_0 + g_dot*v_0

        return r, v

    def calc_SFS(self, body):

        #   Calculate Stumpff function S(z)
        if body.z > 0:
            SFS = (sqrt(body.z) - sin(sqrt(body.z))) / sqrt(body.z) ** 3
        elif body.z < 0:
            SFS = (sinh(sqrt(-body.z)) - sqrt(-body.z)) / sqrt(-body.z) ** 3
        elif body.z == 0:
            SFS = 1 / 6

        return SFS

    def calc_SFC(self, body):

        #   Calculate Stumpff function C(z)
        if body.z > 0:
            SFC = (1 - cos(sqrt(body.z))) / body.z
        elif body.z < 0:
            SFC = (cosh(sqrt(-body.z)) - 1) / -body.z
        elif body.z == 0:
            SFC = 1 / 2

        return SFC

    def find_universal_anomaly(self, body):
        ratio = 1 + self.thres
        mu = self.sys.mu
        r = body.r
        r_mag = sqrt(dot(r, r))
        v = body.v
        X = sqrt(mu) * abs(body.alp) * self.dt
        v_r = dot(v, r)/r_mag
        h = cross(r, v)
        e = 1/mu * (cross(v, h) - mu * r/r_mag)
        e = sqrt(dot(e, e))
        while abs(ratio) > self.thres:

            body.z = body.alp * X ** 2
            SFC = self.calc_SFC(body)
            SFS = self.calc_SFS(body)

            f = r_mag * v_r / sqrt(mu) * X ** 2 * SFC + (1 - body.alp * r_mag) * X ** 3 * SFS + r_mag * X - sqrt(
                mu) * self.dt
            df = r_mag * v_r / sqrt(mu) * X * (1 - body.alp * X ** 2 * SFS) + (1 - body.alp * r_mag) * X \
                 ** 2 * SFC + r_mag

            ratio = f/df
            X = X - ratio


        body.X = X
        self.SFC = SFC
        self.SFS = SFS