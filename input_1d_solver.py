# Domain
x_min = -1
x_max = 2
n_cell = 20000

# Time
t_init = 0.
n_steps = 2
CFL = 0.7

# Numercial scheme
stencil_left = 1
stencil_right = 1

# fluid
GAMMA = 1.4
R_GAS = 290

# init
init_u = 20
init_T = 300
init_P = 101325

# BC
bc_left_T = 1 * init_T
bc_left_u = 1

# BC right
bc_right_P = 0.99* init_P

# Solution output
output_freq = 1
output_dir = 'solut'
