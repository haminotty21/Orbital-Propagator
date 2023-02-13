from math import sqrt, sin, sinh, cos, cosh
from numpy import cross, dot

class Propagator:

    #   alp     : alpha = 1/a
    #   dt      : time delta
    #   epoch   : Starting time
    #   SFC      : Stumpff Function output C
    #   SFS      : Stumpff Function output S
    #   smja    : Semi-Major Axis (a)
    #   X       : Universal Anomaly
    #   z       : "z" = alpha * X ^ 2

    def __init__(self, sys=None, smja=0, timedelta=0.001):
        if sys is None:
            print("No system initialized")
        else:
            self.SFC = 0
            self.SFS = 0

            self.dt = timedelta
            self.thres = 10^-8
            self.sys = sys

            self.propagate_system

    def propagate_system(self):
        for body in self.sys.Bodies:
            r_mag = sqrt(dot(body.r, body.r))
            v_mag = sqrt(dot(body.v, body.v))

            if body.smja is None:
                body.alp = 2/r_mag - v_mag ^ 2/self.sys.mu

            self.find_universal_anomaly(body)
            self.solve_pos_vel(r_mag, v_mag, body)

    def solve_pos_vel(self,r_0_mag, v_0_mag, body):
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
        alpha = body.alp

        #   Calculate f and g
        f = 1 - self.X^2/r_0_mag*self.SFC*self.alp*self.X^2
        g = self.dt - 1/sqrt(self.sys.mu)*self.X^3 *self.SFS*self.alp*self.X^2

        #   Calculate r and r_mg
        r = f*r_0 + g*v_0
        r_mag = sqrt(r * r)

        #   Calculate f_dot and g_dot
        g_dot = 1 - self.X^2/r_mag*self.SFC
        f_dot = sqrt(self.sys.mu) / (r_0_mag * r_mag) * (
                    self.alp * self.X ^ 3 * self.SFS * self.alp * self.X ^ 2 - self.X ^ 2)
        # Calculate v
        v = f_dot*r_0 + g_dot*v_0

        return r, v

    def calc_SFS(self):

        #   Calculate Stumpff function S(z)
        if self.z > 0:
            self.SFS = (sqrt(self.z) - sin(sqrt(self.z))) / sqrt(self.z ^ 3)
        elif self.z < 0:
            self.SFS = (sinh(sqrt(-self.z)) - sqrt(-self.z)) / sqrt(self.z) ^ 3
        elif self.z == 0:
            self.SFS = 1 / 6

    def calc_SFC(self):

        #   Calculate Stumpff function C(z)
        if self.z > 0:
            self.SFC = (1 - cos(sqrt(self.z))) / self.z
        elif self.z < 0:
            self.SFC = (cosh(sqrt(-self.z)) - 1) / -self.ze
        elif self.z == 0:
            self.SFC = 1 / 2

    def solve_keplar_ratio(self, r, v, X, mu):

       f_xi =  r*v/sqrt(mu)*X^2*self.SFC+(1-self.alp*r)*X^3*self.SFS + r*X - sqrt(mu)*self.dt
       df_xi = r*v/sqrt(mu)*X*(1-self.alp*r)*X^2*self.SFC+r
       return f_xi/df_xi

    def find_universal_anomaly(self, body):
        mu = self.sys.mu
        r = body.r
        r_mag = sqrt(dot(r, r))
        v = body.v
        X = sqrt(mu) * abs(body.alp) * self.dt

        h = cross(r, v)
        e = 1/mu * (cross(v, h) - mu * r/r_mag)

        while ratio > self.thres:

            self.z = self.alp * X ^ 2
            self.calc_SFC()
            self.calc_SFS()

            ratio = self.solve_keplar_ratio(self.sys, body)

            X = X - ratio
            self.z = self.alp * X ^ 2

        self.X = X