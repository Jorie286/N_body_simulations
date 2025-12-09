import numpy as np
import sys
from dataclasses import dataclass
from orbit_classes import N_body_orbit
import pandas as pd

# instantiate the orbit integrator
n = -1
k = -3.
ang_mom = 1.
m = 1.
o1 = N_body_orbit(ang_mom, m, n=n, k=k, mu=0.5)

# Plotting time
t_start = 0.
t_end = 10.
delta_t = 0.01

t_pts = np.arange(t_start, t_end+delta_t, delta_t)

# for each stable 3, 4, and 5 body orbit, we define the initial conditions necessary within dataclasses
@dataclass
class ThreeBodyButterflyI:
    r_0 = [1.000001, -1.000002, 0.000003]
    x_1: float = 0.30689
    x_2: float = x_1
    x_3: float = -2*x_1
    y_1: float = 0.12551
    y_2: float = y_1
    y_3: float = -2*y_1
    r_dot_0 = [np.sqrt((x_1**2)+(y_1**2)), np.sqrt((x_2**2)+(y_2**2)), np.sqrt((x_3**2)+(y_3**2))]
    phi_0 = [0.0, 0.0, 0.0]

@dataclass
class ThreeBodyButterflyII:
    r_0 = [1.000001, -1.000002, 0.000003]
    x_1: float = 0.39295
    x_2: float = x_1
    x_3: float = -2*x_1
    y_1: float = 0.09758
    y_2: float = y_1
    y_3: float = -2*y_1
    r_dot_0 = [np.sqrt((x_1**2)+(y_1**2)), np.sqrt((x_2**2)+(y_2**2)), np.sqrt((x_3**2)+(y_3**2))]
    phi_0 = [0.0, 0.0, 0.0]

@dataclass
class ThreeBodyBumblebee:
    r_0 = [1.000001, -1.000002, 0.000003]
    x_1: float = 0.18428
    x_2: float = x_1
    x_3: float = -2*x_1
    y_1: float = 0.58719
    y_2: float = y_1
    y_3: float = -2*y_1
    r_dot_0 = [np.sqrt((x_1**2)+(y_1**2)), np.sqrt((x_2**2)+(y_2**2)), np.sqrt((x_3**2)+(y_3**2))]
    phi_0 = [0.0, 0.0, 0.0]

def get_orbit(r_0, r_dot_0, phi_0):
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
    r_pts, r_dot_pts, phi_pts = o1.solve_leapfrog(t_pts, delta_t, np.array(r_0), np.array(r_dot_0), np.array(phi_0))
    return r_pts, r_dot_pts, phi_pts

# get the names of all stable orbits defined above
orbits = [ThreeBodyButterflyI, ThreeBodyButterflyII, ThreeBodyBumblebee]
for orbit_name in orbits:
    # run the orbit integrator
    r_pts, r_dot_pts, phi_pts = get_orbit(orbit_name.r_0, orbit_name.r_dot_0, orbit_name.phi_0)

    # convert the output to a dataframe for convenience; a new set of columns is added for each body in the system
    df = pd.DataFrame()  # Create an empty DataFrame
    for row in range(len(r_pts[:,0])):
        df[f'body_{row}_r'] = r_pts[row,:]
        df[f'body_{row}_rdot'] = r_dot_pts[row,:]
        df[f'body_{row}_phi'] = phi_pts[row,:]

    # return the data as a dat file for easier use when read in when creating an animation
    df.to_csv(f'{orbit_name.__name__}_output.dat', sep=";")
