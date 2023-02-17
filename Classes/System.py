

class System:
    mu = None
    """
    Parameters
    ----------
    G       : Universal Constant of Gravitation (N*m^2/kg^2)
    MU      : Standard Gravitional Parameter
    Bodies  : Dictionary of all bodies in the System

    """
    def __init__(self, mu=0, mass=0):
        self.mu = mu
        self.Bodies = dict()

    def addbody(self, body, name=""):
        if name == "":
            str = "Body " + len(self.Bodies)
            self.Bodies[str] = body
        else:
            self.Bodies[name] = body
