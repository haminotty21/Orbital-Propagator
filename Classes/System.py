

class System:
    """
    Parameters
    ----------
    G       : Universal Constant of Gravitation (N*m^2/kg^2)
    MU      : Standard Gravitional Parameter
    Bodies  : Dictionary of all bodies in the System

    """

    G = 6.673 * 10^-11
    MU = None

    def __init__(self, mass=0):
        self.MU = mass * self.G
        self.Bodies = dict()

    def addbody(self, body, name=""):
        if name == "":
            str = "Body " + len(self.Bodies)
            self.Bodies[str] = body
        else:
            self.Bodies[name] = body
