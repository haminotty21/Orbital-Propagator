# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
from numpy import array
from datetime import datetime
from Classes.System import System
from Classes.Body import Body
from Classes.Propagator import Propagator
from tests import Test
from Conversions import AU_to_SI, SI_to_AU

"""
    This is the script that drives the creation of the system, it's bodies, and then
    propagates
"""
MASS_SUN = 1  # Units of mass sun
MASS_EARTH = 3300
GRAV_CONST_SUN = 0.000295912  # AU^3/km^3

"""
    Epoch
    UTC: 2017-Jan-01 00:00:00:.0000
    JD : 2457754.5                      au/day
"""
EPOCH = datetime(year=2017, month=1, day=1)

"""
    Initial State Vectors (AU)
"""
POS_EARTH = [-1.796136509111975 * 10 ** -1, 9.667949206859814 * 10 ** -1, -3.668681017942158 * 10 ** -5]
VEL_EARTH = [-1.720038360888334 * 10 ** -2, -3.211186197806460 * 10 ** -3, 7.927736735960840 * 10 ** -7]

"""
    Notes
    _____
        Body1 : ’Oumouamoua (Figure 1) was the first interstellar object to be discovered by Robert Weryk on 
        October 19        , 2017
        
        Body2 : In 2019, Comet 2I/Borisov was the second object discovered to be on a hyperbolic trajectory through our 
        Solar System
    -----
"""

POS_B1 = [3.515868886595499 * 10 ** -2, -3.162046390773074, 4.493983111703389]
VEL_B1 = [-2.317577766980901 * 10 ** -3, 9.843360903693031 * 10 ** -3, -1.541856855538041 * 10 ** -2]

POS_B2 = [7.249472033259724, 14.61063037906177, 14.24274452216359]
VEL_B2 = [-8.241709369476881 * 10 ** -3, -1.156219024581502 * 10 ** -2, -1.317135977481448 * 10 ** -2]


def drive():
    # Testing against values given in Curtis Chapter 3 problem 3.20
    # test = Test()
    # test.propapagator()

    sys1 = System(MASS_SUN)
    sys2 = System(MASS_SUN)

    Earth = Body(name="Earth", _mass=MASS_EARTH, r0=POS_EARTH, v0=VEL_EARTH)
    Oumouamoua = Body(name="Oumouamoua", _mass=0, r0=POS_B1, v0=VEL_B1)
    Borisov = Body(name="Borisov", _mass=0, r0=POS_B2, v0=VEL_B2)

    sys1.addbody(Earth, "Earth")
    sys1.addbody(Oumouamoua, "Oumouamoua")

    sys2.addbody(Earth, "Earth")
    sys2.addbody(Borisov, "Borisov")

    Propagator(sys1, 0, 3600 / 86400)
    Propagator(sys2, 0, 3600 / 86400)


if __name__ == '__main__':
    drive()
