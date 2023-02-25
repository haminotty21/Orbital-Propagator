from datetime import datetime

from Classes.Body import Body
from Classes.Plots import Pork_Chop_Plot
from Classes.System import System
from Classes.Propagator import Propagator
from numpy import array, dot, sqrt, linspace, pi, sin, arctan, tan, cos, arange, empty, NaN
from Classes.Lambert_Solver_Izzio import Izzio_Lambert_Solver
from Conversions import AU_to_SI, SI_to_AU, orbital_elements, get_julian_datetime
from numpy.linalg import norm

GRAV_CONST_SUN_AU = 0.000295912
GRAV_CONST_SUN_SI = 132712 * 10 ** 6
GRAV_CONST_EARTH = 398600.0


class Test:

    def __init__(self):
        ID = 1

    def orbital_elements_test(self):
        mu = 398600
        r = array([-6045, -3490, 2500])
        v = array([-3.457, 6.618, 2.533])

        return orbital_elements(r, v, mu)

    def percent_diff(self, x, y):
        return (norm(x) - norm(y))/norm(y) * 100

    def propapagator(self):
        POS_TEST = array([20000.0, -105000.0, -19000.0])
        VEL_TEST = array([0.9, -3.4, -1.5])

        expected_pos = array([26338.0, -128750.0, -29656.0])       #km
        expected_vel = array([0.86280, -3.2116, -1.4613])   #km/s

        sys_test = System(mu=GRAV_CONST_EARTH)
        test = Body(name="test", r0=POS_TEST, v0=VEL_TEST)

        sys_test.addbody(test, "test")
        Propagator(sys_test, 0, 3600 * 2)

        result_pos = sys_test.Bodies["test"].r
        result_vel = sys_test.Bodies["test"].v

        accur_pos = (sqrt(dot(expected_pos, expected_pos)) - sqrt(dot(result_pos, result_pos)))/sqrt(dot(expected_pos, expected_pos)) * 100
        accur_vel = (sqrt(dot(expected_vel, expected_vel)) - sqrt(dot(result_vel, result_vel)))/sqrt(dot(expected_vel, expected_vel)) * 100
        t = 1

    def Lamberts_Problem(self):
        r1 = array([5644.0, 2830.0, 4170.0])
        r2 = array([-2240.0, 7320.0, 4980.0])
        expected_vel1 = 10.84
        expected_vel2 = 9.970

        v1, v2 = Izzio_Lambert_Solver(r1, r2, mu=398600, t=20*60)
        v1_m = norm(v1)
        v2_m = norm(v2)
        v1_acc = self.percent_diff(expected_vel1, v1)
        v2_acc = self.percent_diff(expected_vel2, v2)
        return v1, v2


    def Porkchop_Plots(self):

        pos_earth = array([-2.292663414632956E-2, -1.015900691571000, 1.443136695093059E-5])
        vel_earth = array([1.691315265663616E-2, -4.490327244856349E-4, 4.905338943924240E-7 ])

        pos_mars = array([1.040530727385831, -9.156827817917115E-1, 1.040530727385831 ])
        vel_mars = array([ 9.775555888077589E-3 , 1.170416979098932E-2, 5.089187517767667E-6])

        t = SI_to_AU(mu=GRAV_CONST_SUN_SI)

        sys = System(GRAV_CONST_SUN_AU)
        delta = 1

        Earth = Body(name="Earth", _mass=0, r0=pos_earth, v0=vel_earth)

        Mars = Body(name="Mars", _mass=0, r0=pos_mars, v0=vel_mars)

        sys.addbody(Earth, "Earth")
        sys.addbody(Mars, "Mars")

        d_start = get_julian_datetime(datetime(2005, 6, 20))
        d_end = get_julian_datetime(datetime(2005, 11, 7))

        a_start = get_julian_datetime(datetime(2005, 12, 1))
        a_end = get_julian_datetime(datetime(2007, 2, 24))

        Propagator(sys, 0, starttime=d_start,
                   stoptime=a_end, dt=delta)
        sys.plot_pos()
        departure_dates = arange(d_start, d_end, delta)
        arrival_dates = arange(a_start, a_end, delta)
        departure_dates = departure_dates.round(1)
        arrival_dates = arrival_dates.round(1)

        v1 = empty([len(arrival_dates), len(departure_dates)])
        v2 = empty([len(arrival_dates), len(departure_dates)])
        v1[:] = NaN
        v2[:] = NaN
        i = -1
        z = 0
        text = open("test.txt", "w")
        for d_date in departure_dates:
            i = i + 1
            print(i)
            j = 0
            for a_date in arrival_dates:
                tof = float(a_date - d_date)
                if tof <= 0:
                    continue
                if i == 106:
                    print(j)
                r1 = sys.Bodies['Earth'].r_array[sys.Bodies['Earth'].time_array.index(d_date)]
                r2 = sys.Bodies['Mars'].r_array[sys.Bodies['Mars'].time_array.index(a_date)]
                if i == 50 and j == 6:
                    print(j)
                v1_vec, v2_vec = Izzio_Lambert_Solver(r1=r1, r2=r2, mu=sys.mu, t=tof, M=0)
                p, v1_vec, q = AU_to_SI(0, v1_vec)
                p, v2_vec, q = AU_to_SI(0, v2_vec)
                v1_m= norm(v1_vec)
                v2_m = norm(v1_vec)
                vt_m = v1 + v2
                if type(v1_vec) is not float and type(v2_vec) is not float:
                    text.write(f"{r1[0]} , {r1[1]}, {r1[2]}, {r2[0]}, {r2[1]}, {r2[2]}, {tof}, {v1_vec[0]} , {v1_vec[1]}, {v1_vec[2]}, {v2_vec[0]} , {v2_vec[1]}, {v2_vec[2]}\n")
                v1[j, i] = norm(v1_vec)
                v2[j, i] = norm(v2_vec)
                j = j + 1

        text.close()
        tot_v_sys1 = v1 + v2
        # p, tot_v_sys1, q = AU_to_SI(v=tot_v_sys1)
        # p, v1_SI, q = AU_to_SI(v=v1)


        fig1, ax1 = Pork_Chop_Plot(departure_dates, arrival_dates, tot_v_sys1, arange(1, 52, .1))
        fig2, ax2 = Pork_Chop_Plot(departure_dates, arrival_dates, v1, arange(1, 50, .1))
        ax1.plot.show()
        ax2.plot.show()