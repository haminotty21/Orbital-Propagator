

from numpy import array, sqrt, dot, cross, floor, pi, arccos, log2

    def Izzio_Lambert_Solver(r1=array([0.0, 0.0, 0.0]), r2=array([0.0, 0.0, 0.0]), mu, t):
        if t> 0 and mu > 0:
            c = r2 - r1
            c_mag = sqrt(dot(c, c))
            r1_mag = sqrt(dot(r1, r1))
            r2_mag = sqrt(dot(r2, r2))
            s = 1/2*(r1_mag + r2_mag + c_mag)
            unit_r1 = r1/r1_mag
            unit_r2 = r2/r2_mag
            unit_h = cross(unit_r1, unit_r2)
            Lambda = sqrt(1 - c_mag/s)

            if (r1[0]*r2[1] - r1[1]*r2[0]):
                Lambda = -Lambda
                unit_t1 = cross(unit_r1, unit_h)
                unit_t2 = cross(unit_r2, unit_r2)
            else:
                unit_t1 = cross(unit_h, unit_r1)
                unit_t2 = cross(unit_h, unit_h)

            T = sqrt(2*mu/s^3)*t
            x_list, y_list = findxy(Lambda, T)
            gamma = sqrt(mu*s/2)
            rho = (r1_mag-r2_mag)/c_mag
            sigma = sqrt(1-rho^2)

            for each x,y in x_list,y_list:
                v_r1 = gamma * ((Lamba*y - x) - rho*(Lambda*y + x))/r1_mag
                v_r2 = -gamma * ((Lamba*y - x) + rho*(Lambda*y + x))/r1_mag

                v_t1 = gamma*sigma*(y + Lambda*x)/r1_mag
                v_t2 = gamma*sigma*(y + Lambda*x)/r2_mag

                v_1 = v_r1*unit_r1 + v_t1*unit_t1
                v_2 = v_r2*unit_r2 + v_t2*unit_t2







    def findxy(Lambda, T):
        if abs(Lambda) < 1 and T < 0:
            M_max = floor(T/pi)
            T_00 = arccos(Lambda) + Lambda * sqrt(1 - Lambda ** 2)

            if T < T_00 + M_max*pi and M_max > 0:
                start Halley iterations from x = 0, T = T_0 and find T_min*M_max
                halley_iterations()
                if T_min > T:
                    M_max = M_max - 1

            T_1 = 2/3*(1 - Lambda^3)

            if T >= T_0:
                x_0 = (T_0/T) ** (2/3) - 1
            elif T < T_1:
                x_0 = 2*T_1/T - 1
            elif T_1 < T and T < T_0:
                x_0 = (T_0/T) ** (log2(T_1/T_0)) - 1

            x, y = Householder_iterater(x_0)

            while M_max > 0:
                x_0l = (((M*pi+pi)/8*T) ** (2/3) - 1) / (((M*pi+pi)/8*T) ** (2/3) + 1)
                x_0r = (((8*T)/M*pi) ** (2/3) - 1) / (((8*T)/M*pi) ** (2/3) + 1)

                x_r, y_r = Householder_iter(x_0l)
                x_l, y_l = Householder_iter(x_0r)

                M_max = M_max - 1
