# Generic formulation for Collocation and Runge-Kutta methods

## Main topic

Let consider the generic non-linear ODE problem

$$
\frac{du}{dt} = f(u,t), \quad u(t_0) = u_0,
$$

solved for one given time step of size $\Delta{t}$
between $[t_0, t]$.

A generic collocation method consider the solution at nodes
$\tau_1, \tau_2, \dots, \tau_M$ used to discretize $[t_0, t]$,
and the representation of the solution on this point using the vector
${\bf u}=[u(t_1), u(t_2), \dots, u(t_M)]^T$,
with $t_m = t_0 + \Delta{t}\tau_m$.
To given nodes is associated a $Q$ matrice used to find the collocation
solution with the system :

$$
(I-\Delta{t}Q\otimes f){\bf u} = {\bf u}_0,
$$

where ${\bf u}_0 = [u_0, u_0, \dots, u_0]^T$ and
$Q\otimes f$ is the non-linear operator evaluating the matrix vector product between $Q$ and
$[f(u(t_1), t_1), [f(u(t_2), t_2),\dots,[f(u(t_M), t_M)]^T$.

Finally, the solution at time $t_1$ is computed using the collocation update (or prolongation) :

$$
u(t) = u_0 + \Delta{t}\sum_{m=1}^{M} \omega_m u(t_m),
$$

where $(\omega_m)_{1\leq m \leq M}$ are the weights associated to the collocation nodes.

There is a link between Runge-Kutta method and collocation method in the unconscious scientific culture of applied mathematics, but it seems that it's full description hase fadded over the years as the people who found it too easy to write it down died progressively. So here is the main question :

> **Can we write any Runge-Kutta method into the same $Q$ matrix formulation as collocation methods ? And if yes (which is very likely), how ?**

From this can be derived additional research topics that could merge technics used for collocation methods and Runge-Kutta ones, for instance :

- _Can we apply a Picard iteration for a Runge-Kutta method, and how stable is it in comparison to the same Runge-Kutta time-stepping ?_
- _Can we develop SDC based iterations using a Runge-Kutta method instead of collocation method for the Q matrix ? And can we use the idea of diagonal optimized SDC on it to develop parallel-in-time RK methods ?_
- _..._

## Repository structure

- [grib](./grib/README.md) : all the notes, reports, articles etc ... on the main and derived topics.
- [python](./python) : python code ([library](./python/collgen/README.md), [scripts](./python/scripts/README.md) and [notebooks](./python/notebook/README.md)).

## Starting points

1. [Runge-Kutta type collocation class](https://github.com/Parallel-in-Time/pySDC/blob/master/pySDC/implementations/sweeper_classes/Runge_Kutta.py) implemented by Thomas in pySDC.
   - Simplified implementation of the $Q$-form generation in [python](./python/collgen/thomas.py)
   - First investigation [script](./python/scripts/01_RK4.py) for RK4