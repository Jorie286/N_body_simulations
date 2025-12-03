# file with the necessary classes to run the calculations for three, four, and five body orbits

# import the necessary packages
import numpy as np
from scipy.integrate import solve_ivp

# NEED TO MODIFY THE DY_DT AND SOLVE_ODE FUCNTIONS FOR THE PROBLEMS I AM SOLVING!!!


class three_body_orbit:
    """
    Equations of motion for a three body gravitational system with equal masses.
    """
    def __init__(self, m=1, G=1): # initialize the necessary constants
        self.m = m
        self.G = G

    def dy_dt(self, t, y):
        """
        This function returns the right-hand side of the diffeq:
        [dr/dt d^2r/dt^2 dphi/dt]

        Parameters
        ----------
        t : float
            time
        y : float
            3-component vector with y[0] = r(t), y[1] = dr/dt, y[2] = dphi/dt

        """
        return [ y[1],
                -1./self.mu * self.Ueff_deriv(y[0]),
                self.ang_mom / (self.mu * y[0]**2) ]

    def solve_ode(self, t_pts, r_0, r_dot_0, phi_0,
                  abserr=1.0e-8, relerr=1.0e-8):
        """
        Solve the ODE given initial conditions.
        Specify smaller abserr and relerr to get more precision.
        """
        y = [r_0, r_dot_0, phi_0]
        solution = solve_ivp(self.dy_dt, (t_pts[0], t_pts[-1]),
                             y, t_eval=t_pts,
                             atol=abserr, rtol=relerr)
        r, r_dot, phi = solution.y
        return r, r_dot, phi


class four_body_orbit:
    """
    Equations of motion for a four body gravitational system with equal masses.
    """
    def __init__(self, m=1, G=1): # initialize the necessary constants
        self.m = m
        self.G = G

    def dy_dt(self, t, y):
        """
        This function returns the right-hand side of the diffeq:
        [dr/dt d^2r/dt^2 dphi/dt]

        Parameters
        ----------
        t : float
            time
        y : float
            3-component vector with y[0] = r(t), y[1] = dr/dt, y[2] = dphi/dt

        """
        return [ y[1],
                -1./self.mu * self.Ueff_deriv(y[0]),
                self.ang_mom / (self.mu * y[0]**2) ]

    def solve_ode(self, t_pts, r_0, r_dot_0, phi_0,
                  abserr=1.0e-8, relerr=1.0e-8):
        """
        Solve the ODE given initial conditions.
        Specify smaller abserr and relerr to get more precision.
        """
        y = [r_0, r_dot_0, phi_0]
        solution = solve_ivp(self.dy_dt, (t_pts[0], t_pts[-1]),
                             y, t_eval=t_pts,
                             atol=abserr, rtol=relerr)
        r, r_dot, phi = solution.y
        return r, r_dot, phi


class five_body_orbit:
    """
    Equations of motion for a five body gravitational system with equal masses.
    """
    def __init__(self, m=1, G=1): # initialize the necessary constants
        self.m = m
        self.G = G

    def dy_dt(self, t, y):
        """
        This function returns the right-hand side of the diffeq:
        [dr/dt d^2r/dt^2 dphi/dt]

        Parameters
        ----------
        t : float
            time
        y : float
            3-component vector with y[0] = r(t), y[1] = dr/dt, y[2] = dphi/dt

        """
        return [ y[1],
                -1./self.mu * self.Ueff_deriv(y[0]),
                self.ang_mom / (self.mu * y[0]**2) ]

    def solve_ode(self, t_pts, r_0, r_dot_0, phi_0,
                  abserr=1.0e-8, relerr=1.0e-8):
        """
        Solve the ODE given initial conditions.
        Specify smaller abserr and relerr to get more precision.
        """
        y = [r_0, r_dot_0, phi_0]
        solution = solve_ivp(self.dy_dt, (t_pts[0], t_pts[-1]),
                             y, t_eval=t_pts,
                             atol=abserr, rtol=relerr)
        r, r_dot, phi = solution.y
        return r, r_dot, phi
