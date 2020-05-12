# Domain
x_min = 0
x_max = 1
n_cell = 200

# Time
t_init = 0.
n_steps = 100000
CFL = 0.7

# Numercial scheme
stencil_left = 1
stencil_right = 1

# fluid
GAMMA = 1.4
R_GAS = 290

# init
init_u = 1.
init_T = 300
init_P = 101325

# BC
bc_left_T = 1.2 * init_T
bc_left_P = 1.0*init_P
bc_left_u = 10. * init_u

# BC right
bc_right_P = 1.*init_P

# Solution output
output_freq = 2000
output_dir = 'solut'
