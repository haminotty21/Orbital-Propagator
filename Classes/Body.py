from matplotlib.pyplot import axes, figure
from numpy import dot, array
from math import sqrt

class Body:
    """
        Variables
        ---------
        smja    :   Semi-major axis (a)
        r       :   Initial position vector
        v       :   Initial velocity vector
    """

    Name = ""
    mass = 0.0
    r = array([0.0, 0.0, 0.0])
    v = array([0.0, 0.0, 0.0])
    z = 0.0
    alp = 0.0
    v_array = 0
    r_array = 0
    time_array = 0


    def __init__(self, name="", _mass=0, r0=[0.0, 0.0, 0.0], v0=[0.0, 0.0, 0.0], ellipticity=0, semimajoraxis=None):
        self.Name = name
        self.mass = _mass
        r0 = array(r0)
        v0 = array(v0)
        self.r = r0
        r_mag = sqrt(dot(r0, r0))
        self.v = v0
        self.e = ellipticity
        self.smja = semimajoraxis
        if semimajoraxis is not None:
            self.alp = 1/semimajoraxis

        self.v_r = dot(self.v, self.r)/r_mag

    def plot_pos(self):
        fig = figure()

        ax = axes(projection='3d')
        ax.plot3D(self.r_array[:, 0], self.r_array[:, 1], self.r_array[:, 2])

