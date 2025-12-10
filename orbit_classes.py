# file with the necessary classes to run the calculations for three, four, and five body orbits

# import the necessary packages
import numpy as np


class N_body_orbit_polar:
    """
    Equations of motion for a three body gravitational system with equal masses.
    """
    def __init__(self, ang_mom, m, n=-1, k=-1, mu=1):
        self.ang_mom = ang_mom
        self.m = m # mass of a body in the system (assume they are all the same to simplify the problem)
        self.n = n
        self.k = k
        self.mu = mu

    def U(self, r, phi):
        """Potential energy of the form U = kr^n."""
        U_tot = 0
        for i in range(len(r)-1):
            x_i, y_i = r[i]*np.cos(phi[i]), r[i]*np.sin(phi[i]) # find the distance between each set of bodies
            x_ip1, y_ip1 = r[i+1]*np.cos(phi[i+1]), r[i+1]*np.sin(phi[i+1])

            d = np.sqrt(((x_i-x_ip1)**2) + ((y_i-y_ip1)**2))

            U_i = self.k * np.abs(d)**self.n
            U_tot += U_i

        # find the looped potential
        x_0, y_0 = r[0]*np.cos(phi[0]), r[0]*np.sin(phi[0])
        x_fin, y_fin = r[-1]*np.cos(phi[-1]), r[-1]*np.sin(phi[-1])

        d = np.sqrt(((x_fin - x_0)**2) + ((y_fin - y_0)**2))

        U_tot += self.k * np.abs(d)**self.n
        return U_tot

    # get the central and effective energy
    def Ucf(self, r):
        """Centrifugal potential energy"""
        return self.ang_mom**2 / (2. * self.mu * r**2)

    def Ueff(self, r, phi):
        """Effective potential energy"""
        return self.U(r, phi) + self.Ucf(r)

    def U_deriv(self, r, phi):
        """dU/dr"""
        U_deriv_tot = 0
        for i in range(len(r)-1):
            x_i, y_i = r[i]*np.cos(phi[i]), r[i]*np.sin(phi[i]) # find the distance between each set of bodies
            x_ip1, y_ip1 = r[i+1]*np.cos(phi[i+1]), r[i+1]*np.sin(phi[i+1])

            d = np.sqrt(((x_i-x_ip1)**2) + ((y_i-y_ip1)**2))

            U_deriv_i = self.n * self.k * np.abs(d) * np.abs(d)**(self.n - 1)
            U_deriv_tot += U_deriv_i

        # find the looped potential
        x_0, y_0 = r[0]*np.cos(phi[0]), r[0]*np.sin(phi[0])
        x_fin, y_fin = r[-1]*np.cos(phi[-1]), r[-1]*np.sin(phi[-1])

        d = np.sqrt(((x_fin - x_0)**2) + ((y_fin - y_0)**2))

        U_deriv_tot += self.n * self.k * np.abs(d) * np.abs(d)**(self.n - 1)
        return U_deriv_tot

    # get the derivative of the central and effective energy
    def Ucf_deriv(self, r):
        """dU_cf/dr"""
        return -2. * self.ang_mom**2 / (2. * self.mu * r**3)

    def Ueff_deriv(self, r, phi):
        """dU_eff/dr"""
        return self.U_deriv(r, phi) + self.Ucf_deriv(r)

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
            3-component vector with y[0] = r(t), y[1] = dr/dt, y[2] = phi(t), y[3] = dphi/dt

        """
        y=np.array(y)
        return [ y[:,1],
                -1./self.mu * self.Ueff_deriv(y[:,0], y[:,2]),
                y[:, 3],
                -2. * self.ang_mom**2 / (2. * self.mu * y[:,0]**3) ]

    # integrate over the differential equation and return the final results
    def solve_leapfrog(self, t_pts, delta_t, r_0, r_dot_0, phi_0, phi_dot_0):
        """
        Solve ODE using leapfrog method.
        """
        # get the intial conditions and add them to the arrays that will be output
        r = np.zeros((len(r_0), len(t_pts)))
        r_dot = np.zeros((len(r_0), len(t_pts)))
        r_dot_half = np.zeros((len(r_0), len(t_pts)))

        phi = np.zeros((len(r_0), len(t_pts)))
        phi_dot = np.zeros((len(r_0), len(t_pts)))
        phi_dot_half = np.zeros((len(r_0), len(t_pts)))

        r[:,0] = r_0
        r_dot[:,0] = r_dot_0
        phi[:,0] = phi_0
        phi_dot[:,0] = phi_dot_0

        # step through the points and get the values
        for i in range(len(t_pts)-1):
            y = np.array([r[:,i], r_dot[:, i], phi[:, i], phi_dot[:, i]]).T

            r_dot_half[:, i] = r_dot[:, i] + self.dy_dt(t_pts[i], y)[1] * (delta_t/2)
            r[:, i+1] = r[:, i] + r_dot_half[:, i] * delta_t

            phi[:, i+1] = phi[:, i] + self.dy_dt(t_pts[i], y)[3] * delta_t

            a_new = -1./self.mu * self.Ueff_deriv(r[:,i+1], phi[:,i+1]) # get the new acceleration before updating r_dot
            r_dot[:, i+1] = r_dot_half[:, i] + a_new * (delta_t/2)

            a_phi_new = -2. * self.ang_mom**2 / (2. * self.mu * r[:,i+1]**3) # get the new angular acceleration before updating phi_dot
            phi_dot[:, i+1] = self.dy_dt(t_pts[i], y)[3] + a_phi_new * (delta_t/2)

        return r, r_dot, phi, phi_dot

    # get the energy of the system
    def energy(self, t_pts, r, r_dot, phi):
        """Evaluate the energy as a function of time"""
        return (self.mu/2.) * r_dot**2 + self.Ueff(r, phi)

class N_body_orbit_cartesian:
    """
    Equations of motion for a three body gravitational system with equal masses.
    """
    def __init__(self, ang_mom, m, n=-1, k=-1, mu=1):
        self.ang_mom = ang_mom
        self.m = m # mass of a body in the system (assume they are all the same to simplify the problem)
        self.n = n
        self.k = k
        self.mu = mu

    def U(self, r, phi):
        """Potential energy of the form U = kr^n."""
        U_tot = 0
        for i in range(len(r)-1):
            x_i, y_i = r[i]*np.cos(phi[i]), r[i]*np.sin(phi[i]) # find the distance between each set of bodies
            x_ip1, y_ip1 = r[i+1]*np.cos(phi[i+1]), r[i+1]*np.sin(phi[i+1])

            d = np.sqrt(((x[i]-x[i+1])**2) + ((y[i]-y[i+1])**2))
           

            U_i = self.k * np.abs(d)**self.n
            U_tot += U_i

        # find the looped potential
        x_0, y_0 = r[0]*np.cos(phi[0]), r[0]*np.sin(phi[0])
        x_fin, y_fin = r[-1]*np.cos(phi[-1]), r[-1]*np.sin(phi[-1])

        d = np.sqrt(((x_fin - x_0)**2) + ((y_fin - y_0)**2))

        U_tot += self.k * np.abs(d)**self.n
        return U_tot

    # get the central and effective energy
    def Ucf(self, r):
        """Centrifugal potential energy"""
        return self.ang_mom**2 / (2. * self.mu * r**2)

    def Ueff(self, r, phi):
        """Effective potential energy"""
        return self.U(r, phi) + self.Ucf(r)

    def U_deriv(self, r, phi):
        """dU/dr"""
        U_deriv_tot = 0
        for i in range(len(r)-1):
            x_i, y_i = r[i]*np.cos(phi[i]), r[i]*np.sin(phi[i]) # find the distance between each set of bodies
            x_ip1, y_ip1 = r[i+1]*np.cos(phi[i+1]), r[i+1]*np.sin(phi[i+1])

            d = np.sqrt(((x_i-x_ip1)**2) + ((y_i-y_ip1)**2))
            print("distance: ", d)

            U_deriv_i = self.n * self.k * np.abs(d) * np.abs(d)**(self.n - 1)
            U_deriv_tot += U_deriv_i

        # find the looped potential
        x_0, y_0 = r[0]*np.cos(phi[0]), r[0]*np.sin(phi[0])
        x_fin, y_fin = r[-1]*np.cos(phi[-1]), r[-1]*np.sin(phi[-1])

        d = np.sqrt(((x_fin - x_0)**2) + ((y_fin - y_0)**2))
        print("distance2: ", d)

        U_deriv_tot += self.n * self.k * np.abs(d) * np.abs(d)**(self.n - 1)
        return U_deriv_tot

    # get the derivative of the central and effective energy
    def Ucf_deriv(self, r):
        """dU_cf/dr"""
        return -2. * self.ang_mom**2 / (2. * self.mu * r**3)

    def Ueff_deriv(self, r, phi):
        """dU_eff/dr"""
        return self.U_deriv(r, phi) + self.Ucf_deriv(r)

    # the differential equation
    def dy_dt(self, t, y):
        """
        This function returns the right-hand side of the diffeq:
        [dr/dt d^2r/dt^2 dphi/dt d^2phi/dt^2]

        Parameters
        ----------
        t : float
            time
        y : float
            3-component vector with y[0] = x(t), y[1] = dx/dt, y[2] = y(t), y[3] = dy/dt

        """
        y=np.array(y)        
        for i in range(len(y[:,0])):
            y[i, 0] = (y[i, 0]**2) + (y[i, 2]**2)  # convert to r for potential calculations
            y[i, 1] = np.sqrt((y[i, 1]**2) + (y[i, 3]**2))  # convert to r_dot for potential calculations
            y[i, 2] = np.arctan2(y[i, 2], y[i, 0])  # convert to phi for potential calculations
            y[i, 3] = np.arctan2(y[i, 3], y[i, 1])  # convert to phi_dot for potential calculations

        return [ y[:,1],
                -1./self.mu * self.Ueff_deriv(y[:,0], y[:,2]),
                y[:, 3],#self.ang_mom**2 / (2. * self.mu * y[:, 0]**2),
                -2. * self.ang_mom**2 / (2. * self.mu * y[:, 0]**3) ]

    # integrate over the differential equation and return the final results
    def solve_leapfrog(self, t_pts, delta_t, x_0, x_dot_0, y_0, y_dot_0):
        """
        Solve ODE using leapfrog method.
        """
        # get the intial conditions and add them to the arrays that will be output
        x = np.zeros((len(x_0), len(t_pts)))
        x_dot = np.zeros((len(x_0), len(t_pts)))
        x_dot_half = np.zeros((len(x_0), len(t_pts)))
        y = np.zeros((len(x_0), len(t_pts)))
        y_dot = np.zeros((len(x_0), len(t_pts)))
        y_dot_half = np.zeros((len(x_0), len(t_pts)))
        x[:,0] = x_0
        x_dot[:,0] = x_dot_0
        y[:,0] = y_0
        y_dot[:,0] = y_dot_0

        # step through the points and get the values
        for i in range(len(t_pts)-1):
            p = np.array([x[:,i], x_dot[:, i], y[:, i], y_dot[:, i]]).T
            x_dot_half[:, i] = x_dot[:, i] + self.dy_dt(t_pts[i], p)[1] * np.cos(self.dy_dt(t_pts[i], p)[3]) * (delta_t/2)
            x[:, i+1] = x[:, i] + x_dot_half[:, i] * delta_t

            y_dot_half[:, i] = y_dot[:, i] + self.dy_dt(t_pts[i], p)[1] * np.sin(self.dy_dt(t_pts[i], p)[3]) * (delta_t/2)
            y[:, i+1] = y[:, i] + y_dot_half[:, i] * delta_t

            x_dot[:, i+1] = x_dot_half[:, i] + self.dy_dt(t_pts[i], p)[1] * np.cos(self.dy_dt(t_pts[i], p)[3]) * (delta_t/2)
            y_dot[:, i+1] = y_dot_half[:, i] + self.dy_dt(t_pts[i], p)[1] * np.sin(self.dy_dt(t_pts[i], p)[3]) * (delta_t/2)

        return x, x_dot, y, y_dot

    # get the energy of the system
    def energy(self, t_pts, x, y, x_dot, y_dot):
        """Evaluate the energy as a function of time"""
        r_dot = np.sqrt(x_dot**2 + y_dot**2)
        r = np.sqrt(x**2 + y**2)
        phi = np.arctan2(y, x)
        return (self.mu/2.) * r_dot**2 + self.Ueff(r, phi)