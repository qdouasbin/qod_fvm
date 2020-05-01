# Domain
x_min = 0
x_max = 1
n_cell = 25

# Time
t_init = 0.
n_steps = 2
CFL = 0.1

# Numercial scheme
stencil_left = 1
stencil_right = 1

# fluid
GAMMA = 1.4
R_GAS = 290

# BC
bc_left_T = 500.
bc_left_P = 1.1e5
bc_left_u = 5

# BC right
bc_right_P = 1e5

# Solution outpur
output_freq = 1
