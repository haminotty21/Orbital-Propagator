# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
from matplotlib.pyplot import axes, figure
from numpy import array, zeros, meshgrid, transpose, empty, NaN, arange
from datetime import datetime, timedelta

from numpy.linalg import norm

from Classes.Lambert_Solver_Izzio import Izzio_Lambert_Solver
from Classes.Plots import Pork_Chop_Plot
from Classes.System import System
from Classes.Body import Body
from Classes.Propagator import Propagator
from tests import Test
from matplotlib import pyplot as plot, figure, axes
from Conversions import AU_to_SI, SI_to_AU, get_julian_datetime, orbital_elements

"""
    This is the script that drives the creation of the system, it's bodies, and then
    propagates
"""
MASS_SUN = 1  # Units of mass sun
MASS_EARTH = 3300
GRAV_CONST_SUN = 0.000295912
GAUSSIAN_GRAV_CONST = 0.01720209895
"""
    Epoch
    UTC: 2017-Jan-01 00:00:00:.0000
    JD : 2457754.5                      au/day
"""
EPOCH = datetime(year=2017, month=1, day=1)

"""
    Initial State Vectors (AU)
"""
POS_EARTH = array([-1.796136509111975 * 10 ** -1, 9.667949206859814 * 10 ** -1, -3.668681017942158 * 10 ** -5])
VEL_EARTH = array([-1.720038360888334 * 10 ** -2, -3.211186197806460 * 10 ** -3, 7.927736735960840 * 10 ** -7])

"""
    Notes
    _____
        Body1 : ’Oumouamoua (Figure 1) was the first interstellar object to be discovered by Robert Weryk on 
        October 19        , 2017
        
        Body2 : In 2019, Comet 2I/Borisov was the second object discovered to be on a hyperbolic trajectory through our 
        Solar System
    -----
"""

POS_B1 = array([3.515868886595499 * 10 ** -2, -3.162046390773074, 4.493983111703389])
VEL_B1 = array([-2.317577766980901 * 10 ** -3, 9.843360903693031 * 10 ** -3, -1.541856855538041 * 10 ** -2])

POS_B2 = array([7.249472033259724, 14.61063037906177, 14.24274452216359])
VEL_B2 = array([-8.241709369476881 * 10 ** -3, -1.156219024581502 * 10 ** -2, -1.317135977481448 * 10 ** -2])


def drive():
    delta = 1
    # Testing against values given in Curtis Chapter 3 problem 3.20
    test = Test()
    # test.propapagator()
    t_v1, t_v2 = test.Lamberts_Problem()
    # h, i, raan, e, aop, theta = test.orbital_elements_test()
    # test.Porkchop_Plots()

    h1, i1, raan1, e1, aop1, theta1 = orbital_elements(POS_B1, VEL_B1, GRAV_CONST_SUN)
    h2, i2, raan2, e2, aop2, theta2 = orbital_elements(POS_B2, VEL_B2, GRAV_CONST_SUN)

    sys1 = System(GRAV_CONST_SUN)
    sys2 = System(GRAV_CONST_SUN)

    Earth = Body(name="Earth", _mass=MASS_EARTH, r0=POS_EARTH, v0=VEL_EARTH)
    Oumouamoua = Body(name="Oumouamoua", _mass=0, r0=POS_B1, v0=VEL_B1)
    sys1.addbody(Earth, "Earth")
    sys1.addbody(Oumouamoua, "Oumouamoua")

    Propagator(sys1, 0, starttime=get_julian_datetime(datetime(2017, 1, 1)), stoptime=get_julian_datetime(datetime(2019, 2,  1)), dt=delta)


    d_start = get_julian_datetime(datetime(2017, 6, 1))
    d_end = get_julian_datetime(datetime(2017, 12,  31))
    a_start = get_julian_datetime(datetime(2017, 8,  1))
    a_end = get_julian_datetime(datetime(2019, 2,  1))

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

    for d_date in departure_dates:
        i = i + 1
        print(i)
        j = 0
        for a_date in arrival_dates:
            tof = float(a_date - d_date)
            if tof <= 0:
                continue

            r1 = sys1.Bodies['Earth'].r_array[sys1.Bodies['Earth'].time_array.index(d_date)]
            r2 = sys1.Bodies['Oumouamoua'].r_array[sys1.Bodies['Oumouamoua'].time_array.index(a_date)]
            v1_vec, v2_vec = Izzio_Lambert_Solver(r1=r1, r2=r2, mu=sys1.mu, t=tof, M=0)
            v1[j, i] = norm(v1_vec)
            v2[j, i] = norm(v2_vec)
            j = j + 1

    tot_v_sys1 = v1 + v2
    p, tot_v_sys1, q = AU_to_SI(v=tot_v_sys1)
    p, v_sys1, q = AU_to_SI(v=tot_v_sys1)

    X = departure_dates
    Y = arrival_dates
    Z = tot_v_sys1
    fig1, ax1 = Pork_Chop_Plot(departure_dates, arrival_dates, tot_v_sys1, arange(1, 50, .1))
    fig2, ax2 = Pork_Chop_Plot(departure_dates, arrival_dates, v_sys1, arange(1, 20, .1))


    Borisov = Body(name="Borisov", _mass=0, r0=POS_B2, v0=VEL_B2)
    Earth = Body(name="Earth", _mass=MASS_EARTH, r0=POS_EARTH, v0=VEL_EARTH)
    sys2.addbody(Earth, "Earth")
    sys2.addbody(Borisov, "Borisov")

    Propagator(sys2, 0, starttime=get_julian_datetime(datetime(2017, 1, 1)), stoptime=get_julian_datetime(datetime(2022, 2,  1)), dt=delta)
    d_start = get_julian_datetime(datetime(2017, 1,  1))
    d_end = get_julian_datetime(datetime(2020, 8,  1))

    a_start = get_julian_datetime(datetime(2019, 6,  1))
    a_end = get_julian_datetime(datetime(2022, 2,  1))
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
    for d_date in departure_dates:
        i = i + 1
        print(i)
        j = 0
        for a_date in arrival_dates:
            tof = a_date - d_date
            if tof <= 0:
                continue
            r1 = sys2.Bodies['Earth'].r_array[sys2.Bodies['Earth'].time_array.index(d_date)]
            r2 = sys2.Bodies['Borisov'].r_array[sys2.Bodies['Borisov'].time_array.index(a_date)]
            v1_vec, v2_vec = Izzio_Lambert_Solver(r1=r1, r2=r2, mu=sys2.mu, t=tof, M=0)
            v1[j, i] = norm(v1_vec)
            v2[j, i] = norm(v2_vec)
            j = j + 1

    tot_v_sys2 = v1 + v2
    p, tot_v_sys2, q = AU_to_SI(v=tot_v_sys2)
    p, v1_SI, q = AU_to_SI(v=v1)


    X = departure_dates
    Y = arrival_dates
    Z = tot_v_sys2
    fig3, ax3 = Pork_Chop_Plot(departure_dates, arrival_dates, tot_v_sys2, arange(1, 60, .1))
    fig4, ax4 = Pork_Chop_Plot(departure_dates, arrival_dates, v1_SI, arange(1, 20, .1))
    ax1.plot.show()
    ax2.plot.show()
    ax3.plot.show()
    ax4.plot.show()

if __name__ == '__main__':
    drive()


