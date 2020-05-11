# Domain
x_min = 0
x_max = 1
n_cell = 20

# Time
t_init = 0.
n_steps = 2000
CFL = 0.7

# Numercial scheme
stencil_left = 1
stencil_right = 1

# fluid
GAMMA = 1.4
R_GAS = 290

# init
init_u = 2.
init_T = 321
init_P = 101325

# BC
bc_left_T = 1.1* init_T
bc_left_P = init_P
bc_left_u = 250

# BC right
bc_right_P = init_P

# Solution output
output_freq = 10
output_dir = 'solut'

