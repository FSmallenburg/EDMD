# Event-driven molecular dynamics code for hard spheres


This is an event-driven molecular dynamics (EDMD) code for hard spheres. It uses neighbor lists [1] and an efficient event calendar [2] for efficiency.

See http://arxiv.org/abs/2201.01100 for a complete description and benchmarking.


## Code variants

Four simulation codes are included, in their respective subfolders. Three of these are variations of a microcanonical simulation code, simulating a constant number $N$ particles in a cubic volume $V$ with periodic boundary conditions at constant energy. 
- **Cell** contains a simulation code where collision checks are based on a cell list.
- **Multi** contains a simualtion code where collision checks are based on a neighbor list, and is more typically more efficient than the **Cell** code at sufficiently high densities.
- **Single** contains a simualtion code where collision checks are based on a neighbor list, and only a single event is scheduled per particle. This is typically more efficient than the **Multi** code.

Additionally, a separate simulation code **Grow** is included in which the particles grow over time until a desired packing fraction is reached. This can be helpful for creating initial configurations.



## Compilation details

Each code consists of a single code file plus a header file. A makefile is included for each code, as well as sample initial configurations. All simulation parameters are defined as global variables near the top of the main code file.


## Simulation units and output

We simulate systems of $N$ hard spheres in a constant volume $V$ in three dimensions. Each particle $i$ has a position $\mathbf{r}_i$, a diameter $\sigma_i$, a mass $m_i$, and a radius $R_i$. 
The simulation code operates in the following units:
-  Lengths are measured in units of the maximum particle size $\sigma$.
-  Mass is measured in units of a reference mass $m$, which is typically chosen to be the mass of a particle with diameter $\sigma$.
-  Time is measured in units of $\tau = \sqrt{\beta m \sigma^2}$. Here, $\beta = 1/k_B T$ with $k_B$ Boltzmann's constant.

Note that the simulation assumes that no particles with a diameter greater than $1 \sigma$ exist in the simulation, and uses assumption in the creation of the cell list. Hence, all particles necessarily have a diameter $\sigma_i \leq \sigma$. By default, the simulation code sets the mass of all particles to be equal to $m$ when loading an initial configuration, but mass is taken into account when determining the effect of collisions. Hence, other choices for the mass can readily be implemented by adapting the initialization functions.

The simulation code measures the pressure $P$ during the simulation, and outputs it in the form of a reduced pressure $P^* = \beta P \sigma^3$. The average pressure in a given time interval $[t_a, t_b]$ is measured via the virial expression 
$P = \rho k_B T + \frac{1}{3V} \frac{\sum  \mathbf{\delta p}_i \cdot \mathbf{r}_{ij}}{t_b - t_a},$ 
where  $\rho = N/V$ is the number density with $V$ the system volume and $N$ the number of particles, and the sum in the last term is taken over all collisions in the time interval $[t_a, t_b]$. For each collision, $\mathbf{r}_{ij}$ is the center-to-center vector connecting the two colliding particles $i$ and $j$, and $\mathbf{\delta p}_i$ is the momentum change of particle $i$ due to the collision. 


Additionally, as a check on conservation of energy, the simulation measures the temperature of the system, using the equipartition theorem 
$\sum_i \frac{1}{2}m v_i^2 = \frac{3N}{2} k_B T.$ 
The temperature is reported in units of $m \sigma^2 / \tau^2 / k_B$, and should be constant during the simulation (up to numerical accuracy). A thermostat function is included in the simulation codes, but is disabled by default. Since the total kinetic energy is a conserved quantity, the temperature remains constant even without a thermostat. Hence, all simulations reported in the main text are performed in the microcanonical ensemble (i.e. at constant total energy).



    
## Snapshot file format

The configuration files from the simulation are written in a simple text-based format, which can contain multiple snapshots per file. For each frame, the format consists of $N+2$ lines (with $N$ the number of particles), as follows:
- One line containing just the number of particles
- One line containing the box size, specifying the box length $L_x$, $L_y$, and $L_z$ along the three axes, separated by whitespace.
- One line per particle containing: a letter indicating particle type, three numbers indicating the real-space particle coordinates, and one number indicating the particle radius. 

The movie files that are created by the simulation code include multiple of these frames consecutively in a single text file. Note that although the code assumes periodic boundary conditions, coordinates of particles that leave the box during the simulation will be printed as being outside of the simulation box, to allow for analysis of long-time dynamics. Hence, any structural analysis or visualization should apply periodic boundary conditions explicitly. 

The simulation codes can read in snapshots in this format as initial configurations. Periodic boundaries will be applied to the snapshot at the start of the simulation.  Note that the simulation code assumes that all box lengths, positions, and radii are given in units of $\sigma$, which is the largest possible particle diameter. Hence, the radius of a particle in the initial configuration should never be given as a number larger than 0.5.

Adaptation to different configuration file formats can be done via modification of the ``loadparticles``, ``write``, and ``outputsnapshot`` functions.

