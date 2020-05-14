Tutorial: Particle gyration
##############################

In this tutorial we use Runko to run a small simulation of a particle orbiting in the plane perpendicular to a constant magnetic field.

.. image:: https://cdn.jsdelivr.net/gh/natj/pb-utilities@master/movies/shock.gif


Running a particle gyration simulation
========================================
You will find the related scripts to run the simulation and the analytical tools in

.. code-block:: bash

   /runko/projects/tests/prtcl_gyration


Configuring the simulation
==========================
To configure the simulation, you can edit the gyration.ini file. As the simulation is very simple, it should finish in seconds.

The following settings are given:



- [io]
   - `outdir`: Defines the output directory name
   - `interval`: No. of steps between output files written
   - `restart`: No. of steps between restart files written
   - `laprestart`: Simulation lap to restart from (-1 = no restart, 0 = automatic, X lap to restart)

Additionally, more advanced options include

   - `stride`: Output image reduce factor
   - `full_interval`: No. of steps between full restart snapshot written to file (-1 = disabled)

- [simulation]
   - `cfl`: The simulation time step in units of grid space
   - `Nt`: The total number of laps in the simulation
   
- [problem]
   - `Nspecies`: Number of species
   - `npasses`: Number of current filter passes, default 0, don't edit.
   
Additionally, more advanced options include

   - `bphi`: external B-field z angle
   - `btheta`: external B-field x-y angle


- [grid]
   - `Nx`, `Ny`, `Nz`: No. of mesh tiles in x,y,z direction, default 1, don't edit.
   - `NxMesh`/`NyMesh`/`NzMesh`: Size of tiles in x,y,z direction
   
Complete grid size is then `Ni*NiMesh`, however, this will be overwritten in the gyration test case with a mesh sized to contain the Larmor gyration of the particle, the radius

.. math::

   r_l = \frac{\gamma m c v}{qB}

where the orthogonal velocity :math:`v` will be calculated from the relativestic factor input :math:`\gamma` (see below).

- [particles]
   - `gamma`: Relativistic factor of the particle, determines speed, increasing the value increases the Larmor radius, and the simulation box must be sized up accordingly.
   - `ppc`: Particles per cell per species
   - `c_omp`: Skindepth resolution

Additionally, more advanced options include

   - `n_test_prtcls`: Number of test particles tracked in simulation, default 2 to track the two manually loaded particles, don't edit.
   


Running the simulation
======================
To run a gyration simulation in Runko using a given particle integrator <pusher>, use the following command:

.. code-block:: bash

   python3 <pusher>_pic.py --conf gyration.ini

To run a gyration simulation in Runko using multiple cores, use the following command:

.. code-block:: bash

   mpirun [-n no_of_cores] python3 <pusher>_pic.py --conf gyration.ini



Using the analysis tools
========================
There are four tools provided for use in studying the results:

1. Particle orbit
-------------------

"<pusher>_orbit.pdf" provides a plot of the particle orbit in XY perpendicular to the input magnetic field.

2. Particle velocity oscillations
-----------------

"<pusher>_vx.pdf" and "<pusher>_vy.pdf" provide plots of the particle x and y velocities vs. time and should reinforce the harmonically oscillating nature of the system.
  
3. Kinetic energy
-------------

"<pusher>_energy.pdf" provides a plot of the particle kinetic energy over time, the behaviour of which will depend on the chosen particle pusher.
