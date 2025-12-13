# N_body_simulations
This repository contains the code from the Physics 5300 final project. It calculates and creates animations of the equilibrium orbits of N body gravitational systems.

To run the orbital integrator used by this project is an adapted leapfrog integrator. The code for this integrator is all contained within the orbit_classes.py file. To run the integrator, the orbit_runs.py file should be used. This file currently contains the initial conditions for a number of equillibrium N body orbits (up to 8 bodies) and the code to run all the orbits. The data from these runs are saved to individual dat files that can be used to animate the orbits or plot them. To run a new orbit, a new dataclass must be added with the initial conditions and name of the orbit. This name must then be added to the orbits list for it to be calculated. To run this file and get the orbits, simply call

  python orbit_runs.py

We analyze the output of these orbits in the animations.ipynb notebook. The notebook demonstrates how to make matplotlib plots of the orbits as well as the code to create Manim animations. NOTE: The animation code must be adapted for each individual orbit as it does not recognize the number of masses in each data file. At the end of the the notebook we added a widget which allows users to view a variety of equillibrium orbits previously rendered as Manim videos.

Last updated: 12/13/2025
