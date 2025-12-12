import numpy as np
import sys
from dataclasses import dataclass
from orbit_classes import N_body_orbit_polar
import pandas as pd

# instantiate the orbit integrator
n = -1
k = -3.
ang_mom = 1.
m = 1.
o1_polar = N_body_orbit_polar(ang_mom, m, n=n, k=k, mu=0.5)

# Plotting time
t_start = 0.
t_end = 10.
delta_t = 0.001
deg_2_rad = np.pi/180

t_pts = np.arange(t_start, t_end+delta_t, delta_t)

# for each stable 3, 4, and 5 body orbit, we define the initial conditions necessary within dataclasses
@dataclass
class ThreeBodyButterflyI:
    x_0 = np.array([1.0, -1.0, 0.0])
    y_0 = np.array([0.0, 0.0, 0.0])
    xdot_1: float = 0.30689
    xdot_2: float = xdot_1
    xdot_3: float = -2*xdot_1
    ydot_1: float = 0.12551
    ydot_2: float = ydot_1
    ydot_3: float = -2*ydot_1
    xdot_0 = np.array([xdot_1, xdot_2, xdot_3])
    ydot_0 = np.array([ydot_1, ydot_2, ydot_3])

    r_0, r_dot_0, phi_0, phi_dot_0 = np.zeros(len(x_0)), np.zeros(len(x_0)), np.zeros(len(x_0)), np.zeros(len(x_0))
    for i in range(len(x_0)):
        r_0[i] = np.sqrt(x_0[i]**2 + y_0[i]**2)
        r_dot_0[i] = np.sqrt(xdot_0[i]**2 + ydot_0[i]**2)
        phi_0[i] = np.arctan2(y_0[i], x_0[i])
        phi_dot_0[i] = np.arctan2(ydot_0[i], xdot_0[i])

@dataclass
class ThreeBodyButterflyII:
    x_0 = np.array([1.0, -1.0, 0.0])
    y_0 = np.array([0.0, 0.0, 0.0])
    xdot_1: float = 0.39295
    xdot_2: float = xdot_1
    xdot_3: float = -2*xdot_1
    ydot_1: float = 0.09758
    ydot_2: float = ydot_1
    ydot_3: float = -2*ydot_1
    xdot_0 = np.array([xdot_1, xdot_2, xdot_3])
    ydot_0 = np.array([ydot_1, ydot_2, ydot_3])
    r_0, r_dot_0, phi_0, phi_dot_0 = np.zeros(len(x_0)), np.zeros(len(x_0)), np.zeros(len(x_0)), np.zeros(len(x_0))
    for i in range(len(x_0)):
        r_0[i] = np.sqrt(x_0[i]**2 + y_0[i]**2)
        r_dot_0[i] = np.sqrt(xdot_0[i]**2 + ydot_0[i]**2)
        phi_0[i] = np.arctan2(y_0[i], x_0[i])
        phi_dot_0[i] = np.arctan2(ydot_0[i], xdot_0[i])

@dataclass
class ThreeBodyBumblebee:
    x_0 = np.array([2.0, -2.0, 0.0])
    y_0 = np.array([0.0, 0.0, 0.0])
    xdot_1: float = 0.18428
    xdot_2: float = xdot_1
    xdot_3: float = -2*xdot_1
    ydot_1: float = 0.58719
    ydot_2: float = ydot_1
    ydot_3: float = -2*ydot_1
    xdot_0 = np.array([xdot_1, xdot_2, xdot_3])
    ydot_0 = np.array([ydot_1, ydot_2, ydot_3])

    r_0, r_dot_0, phi_0, phi_dot_0 = np.zeros(len(x_0)), np.zeros(len(x_0)), np.zeros(len(x_0)), np.zeros(len(x_0))
    for i in range(len(x_0)):
        r_0[i] = np.sqrt(x_0[i]**2 + y_0[i]**2)
        r_dot_0[i] = np.sqrt(xdot_0[i]**2 + ydot_0[i]**2)
        phi_0[i] = np.arctan2(y_0[i], x_0[i])
        phi_dot_0[i] = np.arctan2(ydot_0[i], xdot_0[i])

@dataclass
class threebody_flower:
    r_0 = np.array([1.0, 1.0, 1.0])
    r_dot_0 = np.array([0.0, 0.0, 0.0])
    phi_0 = np.array([0.0, 120.0*deg_2_rad, -120.0*deg_2_rad])
    phi_dot_0 = np.array([20.0, 20.0, 20.0]) #ang_mom / (0.5 * r_0**2)

@dataclass
class fourbody_flower:
    r_0 = np.array([1.0, 1.0, 1.0, 1.0])
    r_dot_0 = np.array([0.0, 0.0, 0.0, 0.0])
    phi_0 = np.array([0.0, 90.0*deg_2_rad, -90.0*deg_2_rad, 180.0*deg_2_rad])
    phi_dot_0 = ang_mom / (0.5 * r_0**2)

@dataclass
class fourbody_limax:
    x_0 = np.array([1.000001, -0.500001, 0.0, -0.500001])
    y_0 = np.array([0.000001, 0.500001, 0.0, -0.500001])
    xdot_1: float = 0.0
    xdot_2: float = -0.5
    xdot_3: float = 0.0
    xdot_4: float = 0.5
    ydot_1: float = 1.5
    ydot_2: float = -1.0
    ydot_3: float = 0.5
    ydot_4: float = -1.0
    xdot_0 = np.array([xdot_1, xdot_2, xdot_3, xdot_4])
    ydot_0 = np.array([ydot_1, ydot_2, ydot_3, ydot_4])

    r_0, r_dot_0, phi_0, phi_dot_0 = np.zeros(len(x_0)), np.zeros(len(x_0)), np.zeros(len(x_0)), np.zeros(len(x_0))
    for i in range(len(x_0)):
        r_0[i] = np.sqrt(x_0[i]**2 + y_0[i]**2)
        r_dot_0[i] = np.sqrt(xdot_0[i]**2 + ydot_0[i]**2)
        phi_0[i] = np.arctan2(y_0[i], x_0[i])
        phi_dot_0[i] = np.arctan2(ydot_0[i], xdot_0[i])

@dataclass
class fivebody_flower:
    r_0 = np.array([1.0, 1.0, 1.0, 1.0, 1.0])
    r_dot_0 = np.array([0.0, 0.0, 0.0, 0.0, 0.0])
    phi_0 = np.array([0.0, 72.0*deg_2_rad, -72.0*deg_2_rad, 144.0*deg_2_rad, -144.0*deg_2_rad])
    phi_dot_0 = ang_mom / (0.5 * r_0**2)

def get_orbit_polar(r_0, r_dot_0, phi_0, phi_dot_0):
    """
    Integrate the motion of the bodies for the type of equillibrium system given.

    Inputs:
        r_0, initial positions (list)
        r_dot_0, initial velocities (list)
        phi_0, initial angles if not in a linear configuration (list)

    Returns:
        r_pts, the positions of the bodies over the orbital integration time (array)
        r_dot_pts, the velocities of the bodies over the integration time (array)
        phi_pts, angular position over the integration time (array)
    """
    r_pts, r_dot_pts, phi_pts, phi_dot_pts = o1_polar.solve_leapfrog(t_pts, delta_t, np.array(r_0), np.array(r_dot_0), np.array(phi_0), np.array(phi_dot_0))
    return r_pts, r_dot_pts, phi_pts, phi_dot_pts

# get the names of all stable orbits defined above
#orbits = [ThreeBodyButterflyI, ThreeBodyButterflyII, ThreeBodyBumblebee, threebody_flower, fourbody_flower, fourbody_limax,
         #fivebody_flower]
orbits = [ThreeBodyBumblebee, threebody_flower, ThreeBodyButterflyI]
for orbit_name in orbits:
    # run the orbit integrator
    r_pts, r_dot_pts, phi_pts, phi_dot_pts = get_orbit_polar(orbit_name.r_0, orbit_name.r_dot_0, orbit_name.phi_0, orbit_name.phi_dot_0)
    # convert the output to a dataframe for convenience; a new set of columns is added for each body in the system
    df = pd.DataFrame()  # Create an empty DataFrame
    for row in range(len(r_pts[:,0])):
        df[f'body_{row}_r'] = r_pts[row,:]
        df[f'body_{row}_rdot'] = r_dot_pts[row,:]
        df[f'body_{row}_phi'] = phi_pts[row,:]
        df[f'body_{row}_phidot'] = phi_dot_pts[row,:]

    # return the data as a dat file for easier use when read in when creating an animation
    df.to_csv(f'{orbit_name.__name__}_output.dat', sep=";")