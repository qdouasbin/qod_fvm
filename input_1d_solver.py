# Domain
x_min = 0
x_max = 1
n_cell = 10

# Time
t_init = 0.
n_steps = 1000
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
init_P = 1.01325e5

# BC
bc_left_T = 900.
bc_left_P = 1.1e5
bc_left_u = 5

# BC right
bc_right_P = 9e4

# Solution output
output_freq = 999
