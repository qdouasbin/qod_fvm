# Domain
x_min = -1
x_max = 2
n_cell = 5

# Time
t_init = 0.
n_steps = 500 * 60
CFL = 0.8

# Simple geometry
simple_geom = True
simple_geom_type = "constant_area"

# Numercial scheme
stencil_left = 1
stencil_right = 1

# fluid
GAMMA = 1.4
R_GAS = 290

# init
init_u = 3
init_T = 300
init_P = 101325

# BC
bc_left_T = 2*init_T
# bc_left_u = 100.0
bc_left_u = 1.5*init_u

# BC right
bc_right_P = init_P

# Solution output
output_freq = 500*60
output_dir = 'solut'
