from datetime import datetime, timedelta
from math import sqrt, sin, sinh, cos, cosh
from numpy import cross, dot, array, zeros, ceil
from numpy.linalg import norm


class Propagator:

    #   alp     : alpha = 1/a
    #   dt      : time delta
    #   epoch   : Starting time
    #   SFC      : Stumpff Function output C
    #   SFS      : Stumpff Function output S
    #   smja    : Semi-Major Axis (a)
    #   X       : Universal Anomaly
    #   z       : "z" = alpha * X ** 2

    def __init__(self, sys=None, smja=0, starttime=datetime(2017, 1, 1), stoptime=datetime(2019, 1,  31), dt=1):
        if sys is None:
            print("No system initialized")
        else:

            self.starttime = starttime
            self.stoptime = stoptime
            self.dt = dt
            self.thres = 10**-8
            self.sys = sys

            self.propagate_system()

    def propagate_system(self):

        timetoprop = int(ceil((self.stoptime - self.starttime)/self.dt))
        for body in self.sys.Bodies.values():
            body.r_array = zeros([timetoprop, 3])
            body.v_array = zeros([timetoprop, 3])
            body.time_array = zeros([timetoprop, 3])
            body.r_array[0] = body.r
            body.v_array[0] = body.v
            #body.time_array[0] = self.starttime
            iter = 1
            cur_time = self.starttime
            body.time_array = [self.starttime + self.dt*x for x in range(timetoprop)]

            for day in range(timetoprop-1):
                # print(iter)
                tof = iter * self.dt
                body.r_array[iter], body.v_array[iter] = self.propagate(body, tof)
                # body.r_array[iter] = body.r
                # body.v_array[iter] = body.v
                iter = iter + 1


    def calc_SFS(self, z):

        #   Calculate Stumpff function S(z)
        if z > 0:
            SFS = (sqrt(z) - sin(sqrt(z))) / sqrt(z) ** 3
        elif z < 0:
            SFS = (sinh(sqrt(-z)) - sqrt(-z)) / sqrt(-z) ** 3
        elif z == 0:
            SFS = 1 / 6

        return SFS

    def calc_SFC(self, z):

        #   Calculate Stumpff function C(z)
        if z > 0:
            SFC = (1 - cos(sqrt(z))) / z
        elif z < 0:
            SFC = (cosh(sqrt(-z)) - 1) / -z
        elif z == 0:
            SFC = 1 / 2

        return SFC

    def propagate(self, body, dt):
        mu = self.sys.mu

        r = body.r
        r_mag = norm(r)

        v = body.v
        v_mag = norm(body.v)
        v_r = dot(v, r)/r_mag

        h = cross(r, v)

        e = r/r_mag - cross(v, h)/mu
        e = norm(e)

        alp = 2 / r_mag - (v_mag ** 2) / mu
        X = sqrt(mu) * abs(alp) * dt
        ratio = 1 + self.thres
        while abs(ratio) > self.thres:

            z = alp * X ** 2
            if z < -100000:
                test = 1
            SFC = self.calc_SFC(z)
            SFS = self.calc_SFS(z)

            f = r_mag * v_r / sqrt(mu) * X ** 2 * SFC + (1 - alp * r_mag) * X ** 3 * SFS + r_mag * X - sqrt(
                mu) * dt
            df = r_mag * v_r / sqrt(mu) * X * (1 - alp * X ** 2 * SFS) + (1 - alp * r_mag) * X \
                 ** 2 * SFC + r_mag

            ratio = f/df
            X = X - ratio

        #   Calculate f and g
        f = 1 - X**2/r_mag*SFC
        g = dt - 1/sqrt(self.sys.mu)*X**3 * SFS

        #   Calculate r and r_mg
        r_0 = r
        r_0_mag = r_mag

        r = f*r + g*v
        r_mag = norm(r)

        #   Calculate f_dot and g_dot
        g_dot = 1 - X**2/r_mag*SFC
        f_dot = sqrt(self.sys.mu) / (r_0_mag * r_mag) * (
                    alp * X ** 3 * SFS - X)
        # Calculate v
        v = f_dot*r_0 + g_dot*v

        return r, v
