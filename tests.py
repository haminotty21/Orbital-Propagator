
from Classes.Body import Body
from Classes.System import System
from Classes.Propagator import Propagator
from numpy import array, dot, sqrt
from Classes.Lambert_Solver_Izzio import Izzio_Lambert_Solver
from Conversions import AU_to_SI, SI_to_AU


GRAV_CONST_EARTH = 398600.0


class Test:

    def __init__(self):
        ID = 1


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
        r1 = array([5644, 2830, 4170])
        r2 = array([-2240, 7320, 4980])
        expected_vel1 = 10.84
        expected_vel2 = 9.970

        Izzio_Lambert_Solver(r1, r2, 398600, 20*60)

