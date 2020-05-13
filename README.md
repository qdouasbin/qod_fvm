# Quasi-1D Finite Volume solver

## Details:

 1. Equations: Euler equations
 2. Numerics:
    - 1st order in time (forward Euler)
    - 1st order in space (Rusanov or Lax-Friedrich flux reconstruction)
 3. Boundary conditions:
    - left: inlet_uty
    - right: outlet_p
    - types: Dirichlet or flux formulation
 4. Input file: TOML format
 5. Fluid properties
    - gaseous
    - calorically perfect gas
    - non-viscous
    
## Theory

A theory guide is available at the following path `docs/source/theory/Quasi-1D FVM solver.md`

## How to use it

The input files are in the TOML format. Here is an example:

```TOML
# Example input file for the Quasi-1D Finite Volume solver

[TimeMarching]
CFL = 0.8
time_init = 0.0
n_steps_max = 30000
time_end = 30.0

[Geometry]
simple_geom = true
simple_geom_type = "constant_area"

[Mesh]
x_min = -1
x_max = 2
n_cells = 5

[Numerics]
stencil = 1

[Fluid]
GAMMA = 1.4
R_GAS = 290

[InitialSolution]
init_type = "Uniform"
init_u = 3.0
init_T = 300.0
init_P = 101325.0

[BoundaryConditions]
    [BoundaryConditions.left]
    type = "inlet_UTY"
    method = "flux"
    T = 650
    u = 4.5
    idx_cell = 0

    [BoundaryConditions.right]
    type = "outlet_P"
    method = "Dirichlet"
    P = 101325
    idx_cell = -1

[IO]
	frequency = 500
	directory = 'solution'
```
