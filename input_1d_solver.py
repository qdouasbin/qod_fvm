# Domain
x_min = 0
x_max = 1
n_cell = 5

# Time
t_init = 0.
n_steps = 100
CFL = 0.1

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
bc_left_T = init_T
bc_left_P = init_P*1.1
bc_left_u = 3

# BC right
bc_right_P = init_P

# Solution output
output_freq = 1
