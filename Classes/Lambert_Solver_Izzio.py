

from numpy import array, sqrt, dot, cross, floor, pi, arccos, log2, sin, sinh, arcsin, arccosh
from numpy.linalg import norm



def dTdx(x, T, Lambda):
    a = 1/(1 - x^2)
    y = sqrt(1.0 - Lambda**2 * a)

    dt = a * (3*T*x - 2 + 2* Lambda**3 * x/y)
    dt_2 = a * (3*T + 5*x*dt + 2*(1- (Lambda**2) )* (Lambda**3) / (y**3))
    dt_3 = a * (7*x*dt_2 + 8*dt - 6*(1 - Lambda**2)*Lambda**5 * x/(y**5))
    return dt, dt_2, dt_3

def Halley_Eqn(x, f, f_1, f_2):
    """
    f   : The function
    f_1 : First derivative of function
    f_2 : Second derivative of function
    x   : Variable to be returned
    :return: x
    """
    x = x - 2 * f * f_1/(2 * f_1 ** 2 - f * f_2)

    return x

def calc_tof(x, Lambda, M):

    y = sqrt(1-Lambda**2 * (1-x^2))
    a = 1/1-x**2                          # This is actually a/am... nondimensional

    if a > 0:            # If elliptic

        alp = 2 * arccos(x)
        beta = 2 * arccos(y)

        tof = a**(3/2) * ((alp - sin(alp)) - (beta - sin(beta) + 2*M*pi))

    elif a < 0:          # If hyperbolic

        alp = 2 * arccosh(x)
        beta = 2 * arccosh(y)

        tof = -2*a**(3 / 2) * (((sinh(alp) - alp) - (sin(beta) - beta)))

def Izzio_Lambert_Solver(r1=array([0.0, 0.0, 0.0]), r2=array([0.0, 0.0, 0.0]), mu=0, t=0):
    if t> 0 and mu > 0:

        c = r2 - r1                     # Chord of Lambert Transfer
        c_m = norm(c, c)                # Chord length
        r1_m = norm(r1, r1)             # r1 magnitude
        r2_m = norm(r2, r2)             # r2 magnitude
        s = 1/2*(r1_m + r2_m + c_m)

        unit_r1 = r1/r1_m
        unit_r2 = r2/r2_m
        unit_h = cross(unit_r1, unit_r2)

        Lambda_sq = 1 - c_m/s
        Lambda = sqrt(Lambda_sq)

        if (r1[0]*r2[1] - r1[1]*r2[0]):
            Lambda = -Lambda
            unit_t1 = cross(unit_r1, unit_h)
            unit_t2 = cross(unit_r2, unit_r2)
        else:
            unit_t1 = cross(unit_h, unit_r1)
            unit_t2 = cross(unit_h, unit_h)

        T = sqrt(2*mu/s^3)*t                    # This is the T* found in Izzio's paper
        x_list, y_list = findxy(Lambda, T)
        gamma = sqrt(mu*s/2)
        rho = (r1_m-r2_m)/c_m
        sigma = sqrt(1-rho^2)

        for x, y in x_list, y_list:
            v_r1 = gamma * ((Lambda*y - x) - rho*(Lambda*y + x))/r1_m
            v_r2 = -gamma * ((Lambda*y - x) + rho*(Lambda*y + x))/r1_m

            v_t1 = gamma*sigma*(y + Lambda*x)/r1_m
            v_t2 = gamma*sigma*(y + Lambda*x)/r2_m

            v_1 = v_r1*unit_r1 + v_t1*unit_t1
            v_2 = v_r2*unit_r2 + v_t2*unit_t2


def Householder_Iterator(x_old, Lambda, T, M):
    thres = 10^-13
    err = 1
    while thres < err:
        tof = calc_tof(x, Lambda, M)
        f = tof - T
        f_1, f_2, f_3 = dTdx(x, T, Lambda)
        num = f_1 ** 2 - f*f_2/2
        denom = f_1* (f_1**2 - f*f_2) + f_3*(f**2)/6
        x = x_old - f * num/denom
        err = abs(x_old - x)
        x_old = x

    return x


def findxy(Lambda, T):
    x_list = list()
    y_list = list()
    thres = 10^-13
    err = 1

    if abs(Lambda) < 1 and T < 0:
        M_max = floor(T/pi)
        T_00 = arccos(Lambda) + Lambda * sqrt(1 - Lambda ** 2)
        T_0 = T_00 + M_max * pi

        if T < T_00 + M_max*pi and M_max > 0:
            T_min = T_0
            x = 0
            """ Haley Iterations """
            while thres < err:
                dT, dT_2, dT_3 = dTdx(x, T_min, Lambda)
                T_min = Halley_Eqn(T_min, dT, dT_2, dT_3)

            if T_min > T:
                M_max = M_max - 1

        T_1 = 2/3*(1 - Lambda ** 3)

        if T >= T_0:
            x_0 = (T_0/T) ** (2/3) - 1
        elif T < T_1:
            x_0 = 2*T_1/T - 1
        elif T_1 < T and T < T_0:
            x_0 = (T_0/T) ** (log2(T_1/T_0)) - 1

        # calc_tof(x_0, T, M_max)
        x = Householder_Iterator(x_0, Lambda, T, M_max)
        while M_max > 0:

            x_0l = (((M_max*pi+pi)/8*T) ** (2/3) - 1) / (((M_max*pi+pi)/8*T) ** (2/3) + 1)
            x_0r = (((8*T)/M_max*pi) ** (2/3) - 1) / (((8*T)/M_max*pi) ** (2/3) + 1)

            x_r, y_r = Householder_Iterator(x_0l, Lambda, T, M_max)
            x_l, y_l = Householder_Iterator(x_0r, Lambda, T, M_max)
            M_max = M_max - 1


        for x in x_list:
            y_list.append(sqrt(1 - Lambda**2*(1-x**2)))

    return x_list, y_list


