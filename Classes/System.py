from matplotlib.pyplot import figure, axes, xlim, ylim




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

    def plot_pos(self):
        fig = figure()
        mx = 0
        my = 0
        mz = 0
        ax = axes(projection='3d')
        ax.plot(0, 0, 0, color='yellow', marker='o', markersize=50)
        for body in self.Bodies.values():
            if mx < max(body.r_array[:, 0]):
                mx = max(body.r_array[:, 0])
            if my < max(body.r_array[:, 1]):
                my = max(body.r_array[:, 1])
            if mz < max(body.r_array[:, 2]):
                mz = max(body.r_array[:, 2])
            ax.plot3D(body.r_array[:, 0], body.r_array[:, 1], body.r_array[:, 2])
        xlim([-mx, mx])
        ylim([-my, my])
        ax.set_zlim(-mz, mz)
        ax.set_aspect('equal', adjustable='box')
        fig.show()
        test = 1