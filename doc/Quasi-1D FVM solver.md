# Quasi-1D FVM solver

## Euler equation

### Continuity 

The continuity equation reads:
$$
\frac{\partial}{\partial t} \iiint_V \rho dV + \iint_S \rho \mathbf{u}\cdot \mathbf{dS} = 0
$$
Considering a quasi-1D control volume we have:

![](C:\Users\quent\PycharmProjects\one_dimensional_fv_solver\doc\cv.svg)

First, we note that since $dAdx$ is small, we have the following simplification for the first term: 
$$
lim_{dx \to 0} \frac{\partial}{\partial t} \iiint_V \rho dV = \frac{\partial}{\partial t} \left( \rho A dx \right)
$$
The second term is:
$$
\iint_S \rho \mathbf{u}\cdot \mathbf{dS} = -\rho u A + (\rho +d\rho)(u + du)(A+dA)
$$
Keep only the first order terms, _i.e._ neglecting the cross products $dq_1dq_2$, we have:
$$
\iint_S \rho \mathbf{u}\cdot \mathbf{dS} = \rho u dA + \rho A du + u A d\rho  = d \left(  \rho A u\right)
$$
Combining Eqs. (2) and (4) and dividing by $dx$ leads to the quasi-1D continuity equation:
$$
\frac{\partial }{\partial t} \left( \rho A\right) + \frac{\partial }{\partial x} \left( \rho A u\right) = 0
$$
We note that this equation reduces to the 1D Euler equation for constant cross sections.

### Momentum

The momentum conservation over the control volume is is written as:
$$
\frac{\partial}{\partial t}\iiint_V \left( \rho u \right)dV + \iint_S \left( \rho u  \mathbf{u}\right) \cdot \mathbf{dS} = - \iint_S \frac{\partial }{\partial x} \left( p \mathbf{dS} \right)
$$
Similarly to the continuity equation, the first term can be simplified as:
$$
\frac{\partial}{\partial t}\iiint_V \left( \rho u \right)dV = \frac{\partial}{\partial t} \left( \rho u A dx \right)
$$
The second term:
$$
\iint_S \left( \rho u  \mathbf{u}\right) \cdot \mathbf{dS} = -\rho u^2 A + (\rho +d\rho)(u + du)^2(A+dA) = \rho u^2 A
$$
The right hand side is:
$$
- \iint_S \frac{\partial }{\partial x} \left( p \mathbf{dS} \right) = -pA + (p + dp)(A+dA) - 2 p \frac{dA}{2} = -A d p
$$
Combining this terms and dividing by $dx$ leads to the quasi-1D momentum equation:
$$
\frac{\partial}{\partial t} \left( \rho u A \right) + \frac{\partial}{\partial x} \left( \rho u^2 A\right) = - A \frac{\partial p}{\partial x}
$$

### Energy

Following the derivation of the quasi-1D continuity and energy equations, the quasi-1D energy equation reads:
$$
\frac{\partial}{\partial t} \left[ \rho \left( e + \frac{u^2}{2}\right) \right] + \frac{\partial}{\partial x} \left[ \rho u A \left( e + \frac{u^2}{2}\right) \right] 
=
- \frac{\partial}{\partial x} \left( p A u \right)
$$


## Euler equation: Matrix form

In the matrix form, the Quasi-1D Euler equations can be written as follows:
$$
\frac{\partial \mathbf{U}}{\partial t} + \frac{\partial \mathbf{F}}{\partial x} = - \mathbf{S}
$$
With
$$
\mathbf{U}
 = 
 \begin{bmatrix}
\rho A\\
\rho u A\\
\rho E A
\end{bmatrix}
\quad \text{with} \quad E = e + \frac{u_2}{2}
$$

$$
\mathbf{F}
 = 
 \begin{bmatrix}
\rho u A\\
\rho u^2 A\\
\rho E u A
\end{bmatrix}
$$

$$
\mathbf{S}
 = 
 \begin{bmatrix}
0\\
A \frac{\partial p}{\partial x}\\
\frac{\partial}{\partial x} \left( p u A \right)
\end{bmatrix}
$$

Moving the right hand sides to the fluxes we have:
$$
\frac{\partial \mathbf{U}}{\partial t} + \frac{\partial \mathbf{F}}{\partial x} = - \mathbf{S}
$$
With
$$
\mathbf{U} =  \begin{bmatrix}\rho A\\\rho u A\\\rho E A\end{bmatrix}\quad \text{with} \quad E = e + \frac{u_2}{2}
$$

$$
\mathbf{F} =  A \begin{bmatrix}
\rho u \\
\rho u^2 + p\\
u \left( \rho E + p\right)
\end{bmatrix}
$$

$$
\mathbf{S} =  
\begin{bmatrix}
0\\0\\0
\end{bmatrix}
$$

From this we have:
$$
\frac{\partial \mathbf{U}}{\partial t} = - \left[ \frac{\partial \mathbf{F}}{\partial x} + \mathbf{S} \right]
$$
