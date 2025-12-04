# file with the necessary classes to run the calculations for three, four, and five body orbits

# import the necessary packages
import numpy as np


class N_body_orbit:
    """
    Equations of motion for a three body gravitational system with equal masses.
    """
    def __init__(self, ang_mom, m, n=-1, k=1, mu=1):
        self.ang_mom = ang_mom
        self.m = m # mass of a body in the system (assume they are all the same to simplify the problem)
        self.n = n
        self.k = k
        self.mu = mu

    def U(self, r):
        """Potential energy of the form U = kr^n."""
        U_tot = 0
        for i in range(len(r)-1):
            U_i = self.k * np.abs(r[i] - r[i+1])**self.n
            U_tot += U_i
        return U_tot

    # get the central and effective energy
    def Ucf(self, r):
        """Centrifugal potential energy"""
        return self.ang_mom**2 / (2. * self.mu * (sum([self.m * r[i] for i in range(len(r)-1)])/(self.m**(len(r))))**2)

    def Ueff(self, r):
        """Effective potential energy"""
        return self.U(r) + self.Ucf(r)

    def U_deriv(self, r):
        """dU/dr"""
        U_deriv_tot = 0
        for i in range(len(r)-1):
            U_deriv_i = self.n * self.k * (r[i] - r[i+1]) * np.abs(r[i] - r[i+1])**(self.n - 1)
            U_deriv_tot += U_deriv_i
        return U_deriv_tot

    # get the derivative of the central and effective energy
    def Ucf_deriv(self, r):
        """dU_cf/dr"""
        return -2. * self.ang_mom**2 / (2. * self.mu * (sum([self.m * r[i] for i in range(len(r)-1)])/(self.m**(len(r))))**3)

    def Ueff_deriv(self, r):
        """dU_eff/dr"""
        return self.U_deriv(r) + self.Ucf_deriv(r)

    # the differential equation
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
        y=np.array(y)
        return [ y[:,1],
                -1./self.mu * self.Ueff_deriv(y[:,0]),
                self.ang_mom / (self.mu * y[:,0]**2) ]

    # integrate over the differential equation and return the final results
    def solve_leapfrog(self, t_pts, delta_t, r_0, r_dot_0, phi_0):
        """
        Solve ODE using leapfrog method.
        """
        # get the intial conditions and add them to the arrays that will be output
        r = np.zeros((len(r_0), len(t_pts)))
        r_dot = np.zeros((len(r_0), len(t_pts)))
        r_dot_half = np.zeros((len(r_0), len(t_pts)))
        phi = np.zeros((len(r_0), len(t_pts)))
        r[:,0] = r_0
        r_dot[:,0] = r_dot_0
        phi[:,0] = phi_0

        # step through the points and get the values
        for i in range(len(t_pts)-1):
            y = np.array([r[:,i], r_dot[:, i], phi[:, i]]).T
            r_dot_half[:, i] = r_dot[:, i] + self.dy_dt(t_pts[i], y)[1] * (delta_t/2)
            r[:, i+1] = r[:, i] + r_dot_half[:, i] * delta_t
            phi_dot_half = self.ang_mom / (self.mu * r[:, i]**2)
            phi[:, i+1] = phi[:, i] + phi_dot_half * delta_t

            yplus = np.array([r[:, i+1], r_dot[:, i+1], phi[:, i+1]]).T
            r_dot[:, i+1] = r_dot_half[:, i] + self.dy_dt(t_pts[i], yplus)[1] * (delta_t/2)
        return r, r_dot, phi

    # get the energy of the system
    def energy(self, t_pts, r, r_dot):
        """Evaluate the energy as a function of time"""
        return (self.mu/2.) * r_dot**2 + self.Ueff(r)
