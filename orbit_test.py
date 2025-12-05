import matplotlib.pyplot as plt
import numpy as np
from orbit_classes import N_body_orbit

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

# Initial conditions
r_0 = np.array([1.001, -1.002, -0.5])
r_dot_0 = np.array([0.0, 0.0, 0.0])

deg_2_rad = np.pi/180
phi_0 = np.array([0.0, 0.0, 0.0])
r_pts_l, r_dot_pts_l, phi_pts_l = o1.solve_leapfrog(t_pts, delta_t, r_0, r_dot_0, phi_0)

c = o1.ang_mom**2 / (np.abs(o1.k) * o1.mu)
epsilon = c / r_0 - 1.
energy_0 = o1.mu/2. * r_dot_0**2 + o1.Ueff(r_0)

# print(f'energy = {energy_0:.2f}')
# print(f'eccentricity = {epsilon:.2f}')

fig_4 = plt.figure(figsize=(8,8))

overall_title = 'Gravitational orbit:  ' + \
                rf' $n = {o1.n},$' + \
                rf' $k = {o1.k:.1f},$' + \
                rf' $l = {o1.ang_mom:.1f},$' + \
                rf' $r_0 = {r_0},$' + \
                rf' $\dot r_0 = {r_dot_0},$' + \
                rf' $\phi_0 = {phi_0}$' + \
                '\n'     # \n means a new line (adds some space here)
fig_4.suptitle(overall_title, va='baseline')

ax_4c = fig_4.add_subplot(1,1,1)
ax_4c.plot(r_pts_l[0,:]*np.cos(phi_pts_l[0,:]), r_pts_l[0,:]*np.sin(phi_pts_l[0,:]), color='blue', label = "Leapfrog")
ax_4c.plot(r_pts_l[1,:]*np.cos(phi_pts_l[1,:]), r_pts_l[1,:]*np.sin(phi_pts_l[1,:]), color='red', label = "Leapfrog")
ax_4c.plot(r_pts_l[2,:]*np.cos(phi_pts_l[2,:]), r_pts_l[2,:]*np.sin(phi_pts_l[2,:]), color='green', label = "Leapfrog")
#ax_4c.plot(r_pts_l[3,:]*np.cos(phi_pts_l[3,:]), r_pts_l[3,:]*np.sin(phi_pts_l[3,:]), color='green', label = "Leapfrog")
ax_4c.set_xlabel(r'$x$')
ax_4c.set_ylabel(r'$y$')
ax_4c.set_aspect(1)
ax_4c.set_title('Cartesian plot')
ax_4c.set_xlim(-10,10)
ax_4c.set_ylim(-10,10)

fig_4.tight_layout()

fig_4.savefig("Orbits_plots.png", bbox_inches = "tight")


fig_4 = plt.figure(figsize=(8,8))

overall_title = 'Gravitational orbit:  ' + \
                rf' $n = {o1.n},$' + \
                rf' $k = {o1.k:.1f},$' + \
                rf' $l = {o1.ang_mom:.1f},$' + \
                rf' $r_0 = {r_0},$' + \
                rf' $\dot r_0 = {r_dot_0},$' + \
                rf' $\phi_0 = {phi_0}$' + \
                '\n'     # \n means a new line (adds some space here)
fig_4.suptitle(overall_title, va='baseline')

ax_4d = fig_4.add_subplot(1,1,1, polar=True)

for i, l in enumerate(r_pts_l): # fix any negative radii so that they don't interfere with the polar plot
    for j, val in enumerate(l):
        if val<0:
            phi_pts_l[i, j] = phi_pts_l[i, j]-np.pi
            r_pts_l[i, j] = np.abs(r_pts_l[i, j])

ax_4d.plot(phi_pts_l[0,:], r_pts_l[0,:], color='blue', label = "Leapfrog")
ax_4d.plot(phi_pts_l[1,:], r_pts_l[1,:], color='red', label = "Leapfrog")
ax_4d.plot(phi_pts_l[2,:], r_pts_l[2,:], color='green', label = "Leapfrog")
ax_4d.set_title('Polar plot', pad=20.)
#ax_4d.set_ylim(0,4)

fig_4.tight_layout()

fig_4.savefig("Orbits_plots_polar.png", bbox_inches = "tight")
