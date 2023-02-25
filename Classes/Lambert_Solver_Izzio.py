

from numpy import array, sqrt, dot, cross, floor, pi, arccos, log2, sin, sinh, arcsin, arccosh, log, linspace, cos, \
    arctan, tan, arcsinh, arctan2
from numpy.array_api import atan2
from numpy.linalg import norm










def Izzio_Lambert_Solver(r1=array([0.0, 0.0, 0.0]), r2=array([0.0, 0.0, 0.0]), mu=0, t=0, M =-1):
    if t > 0 and mu > 0:

        r1_m = norm(r1)             # r1 magnitude
        r2_m = norm(r2)             # r2 magnitude

        c = r2 - r1                 # Chord of Lambert Transfer
        c_m = norm(c)               # Chord length

        s = 1/2*(r1_m + r2_m + c_m)

        u_r1 = r1/r1_m
        u_r2 = r2/r2_m
        u_h = cross(u_r1, u_r2)

        # Taking the norm again to ensure the magnitude of these vectors are actually one.
        u_r1 = u_r1 / norm(u_r1)
        u_r2 = u_r2 / norm(u_r2)
        u_h = u_h / norm(u_h)

        Lambda_sq = 1 - c_m/s
        Lambda = sqrt(Lambda_sq)

        if (r1[0]*r2[1] - r1[1]*r2[0]) < 0:
            Lambda = -Lambda
            u_t1 = cross(u_r1, u_h)
            u_t2 = cross(u_r2, u_h)
        else:
            u_t1 = cross(u_h, u_r1)
            u_t2 = cross(u_h, u_r2)

        u_t1 = u_t1 / norm(u_t1)
        u_t2 = u_t2 / norm(u_t2)

        T = sqrt(2*mu/s**3)*t

        x_list, y_list = findxy(Lambda, T, M)

        if x_list is None and y_list is None:
            return 0, 0
        gamma = sqrt(mu*s/2)
        rho = (r1_m-r2_m)/c_m
        sigma = sqrt(1-rho ** 2)

        for x, y in zip(x_list, y_list):
            v_r1 = gamma * ((Lambda*y - x) - rho*(Lambda*y + x))/r1_m
            v_r2 = -gamma * ((Lambda*y - x) + rho*(Lambda*y + x))/r2_m

            v_t1 = gamma*sigma*(y + Lambda*x)/r1_m
            v_t2 = gamma*sigma*(y + Lambda*x)/r2_m

            v_1 = v_r1 * u_r1 + v_t1 * u_t1
            v_2 = v_r2 * u_r2 + v_t2 * u_t2

        return v_1, v_2


def findxy(Lambda, T, M):
    x_list = list()
    y_list = list()
    thres = 10**-13
    err = 1

    if M == -1:
        M_max = floor(T/pi)
    else:
        M_max = M

    T_00 = arccos(Lambda) + Lambda * sqrt(1 - Lambda ** 2)
    T_0 = T_00 + M_max * pi
    T_1 = 2 / 3 * (1 - Lambda ** 3)

    if T < T_0 and M_max > 0:
        T_min = T_0
        x_old = 0
        iter = 0
        """ Haley Iterations """
        while thres < err:
            dT, dT_2, dT_3 = dTdx(x, T_min, Lambda)
            x_new = x_old - 2 * dT * dT_2/(2 * dT_2 ** 2 - dT * dT_3)

            err = abs(x_old - x_new)
            if (err < thres) or iter > 15:
                break
            if T_min > T:
                M_max = M_max - 1

            calc_tof(x_new, Lambda, M_max)
            iter = iter + 1


    if T >= T_0:
        x_0 = (T_0/T) ** (2/3) - 1
    elif T < T_1:
        x_0 = 5/2 * T_1*(T_1 - T)/(T * (1-Lambda**5)) + 1
    elif T_1 < T and T < T_0:
        x_0 = (T_0/T) ** (log2(T_1/T_0)) - 1

    tmp = Householder_Iterator(x_0, Lambda, T, M_max)

    if tmp is None:
        return None, None

    x_list.append(tmp)
    while M_max > 0:

        x_0l = (((M_max*pi+pi)/8*T) ** (2/3) - 1) / (((M_max*pi+pi)/8*T) ** (2/3) + 1)
        x_0r = (((8*T)/M_max*pi) ** (2/3) - 1) / (((8*T)/M_max*pi) ** (2/3) + 1)

        x_r, y_r = Householder_Iterator(x_0l, Lambda, T, M_max)
        x_l, y_l = Householder_Iterator(x_0r, Lambda, T, M_max)
        M_max = M_max - 1


    for x in x_list:
        y_list.append(sqrt(1 - Lambda**2 * (1-x**2)))

    return x_list, y_list


def Householder_Iterator(x_old, Lambda, T, M):
    thres = 10**-8
    err = 1
    iter = 0
    while thres < err:
        tof = calc_tof(x_old, Lambda, M)
        f = tof - T
        f_1, f_2, f_3 = dTdx(x_old, tof, Lambda)
        num = f_1 ** 2 - f*f_2/2
        denom = f_1 * (f_1**2 - f * f_2) + f_3 * f**2/6
        x = x_old - f * num/denom
        err = abs(x_old - x)
        if x < -1:
            test = 1
        x_old = x
        iter = iter + 1
        if iter > 16:
            return None

    return x


def dTdx(x, T, Lambda):
    y = sqrt(1.0 - Lambda**2 * (1 - x**2))

    dt = 1/(1 - x**2) * (3*T*x - 2 + 2 * Lambda ** 3 * x/y)
    dt_2 = 1/(1 - x**2) * (3*T + 5*x*dt + 2 * (1 - Lambda**2) * (Lambda**3) / (y**3))
    dt_3 = 1/(1 - x**2) * (7*x*dt_2 + 8*dt - 6*(1 - Lambda**2)*Lambda**5 * x/(y**5))
    return dt, dt_2, dt_3


def calc_tof(x, Lambda, M):
    y = sqrt(1 - Lambda ** 2 * (1 - x ** 2))
    battin = 0.01
    lagrange = 0.2
    dist = abs(x-1)


    if dist < lagrange and dist > battin:
        a = 1/(1-x**2)
        if a > 0:            # If elliptic
            alp = 2 * arccos(x)
            beta = 2 * arccos(y)
            tof = a**(3/2) * ((alp - sin(alp)) - (beta - sin(beta) + 2*M*pi))
        elif a < 0:          # If hyperbolic
            alp = 2 * arccosh(x)
            beta = 2 * arccosh(y)
            tof = 2*(-a)**(3 / 2) * (((sinh(alp) - alp) - (sin(beta) - beta)))

        return tof

    elif dist < battin:           # Use battin
        return Battin_tof(Lambda, x, M)
    else:
        return Lancaster_tof(Lambda, x, M)


def Battin_tof(Lambda, x, M):
    y = sqrt(1 + Lambda**2 * (x**2 - 1))
    rho = abs(x ** 2 - 1)
    eta = y - Lambda * x
    S1 = 0.5 * (1 - Lambda - x * eta)
    Q = 4 / 3 * (Hypergeometric_Battin(S1))
    tof = 1 / 2 * (eta ** 3 * Q + 4 * Lambda * eta) + M*pi / rho ** 1.5

    return tof


def Hypergeometric_Battin(z):
    A_old = 1.0
    S_old = 1.0

    err = 1.0
    thres = 10**-13
    k = 0
    while abs(err) > thres:

        A_new = A_old * (3.0 + k) * (1.0 + k) / (2.5 + k) * z/(k+1)
        S_new = S_old + A_new

        A_old = A_new
        S_old = S_new
        err = abs(A_old)
        k  = k + 1

    return S_old




def Lancaster_tof(Lambda, x, M):
    y = sqrt(1-Lambda**2*(1-x**2))

    if x < 1:
        psi_sin = arcsin((y - x*Lambda) * sqrt(1 - x ** 2))
        psi_cos = arccos(x*y + Lambda * (1 - x**2))
        psi_tan = arctan2(sin(psi_sin), cos(psi_cos))
        psi = psi_tan
        # print(str(psi_tan) + "==" + str(psi_cos))

    if x >= 1:
        psi_sin = arcsinh((y - x * Lambda) * sqrt(x ** 2) - 1)
        psi_cos = arccosh(x * y - Lambda * (x ** 2 - 1))
        psi = psi_cos
      #  psi_tan = atanh2(sinh(psi_sin), cosh(psi_cos))


    tof = 1 / (1 - x ** 2) * ((psi + M * pi) / sqrt(abs(1 - x ** 2)) - x + Lambda * y)
    return tof











