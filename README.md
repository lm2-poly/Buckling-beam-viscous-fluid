# Buckling-beam-viscous-fluid


This code solves the equations of a beam extruding in a highly viscous fluid. This code is used to generate the numerical results published in the reference below :

Gosselin, F. P., Neetzow, P., & Paak, M. (2014). **Buckling of a beam extruded into highly viscous fluid**. *Physical Review E*, 90(5), 052718.

https://link.aps.org/doi/10.1103/PhysRevE.90.052718

https://arxiv.org/abs/1406.2365 **&larr; Preprint**

In a nutshell, the code uses a backward Euler algorithm to integrate in time the large-deflection dynamics of the rod extruded in stagnant viscous fluid. At each time step, the BVP4C Matlab function is called to solve implicitly the equations of the beam on the Lagrangian domain. The length of the beam is scaled to always be 1. 

To launch the code, download all files and run extrudingbeamviscousflow.m with Matlab.

THIS CODE WAS DEVELOPPED BY FREDERICK P. GOSSELIN AT POLYTECHNIQUE MONTREAL IN 2014.
