# Generic formulation for Collocation and Runge-Kutta methods

## Main topic

Let consider the generic non-linear ODE problem

$$
\frac{du}{dt} = f(u,t), \quad u(0) = u_0,
$$

solved for one given time step between $[0, \Delta{t}]$.

A generic collocation method consider the solution at nodes
$\tau_1, \tau_2, \dots, \tau_M$ distributed on $[0, \Delta{t}]$,
and the representation of the solution on this point using the vector
${\bf u}=[u(\tau_1), u(\tau_2), \dots, u(\tau_M)]^T$.
To given nodes is associated a $Q$ matrice used to find the collocation
solution with the system :

$$
(I-\Delta{t}Q\otimes f){\bf u} = {\bf u}_0,
$$

where ${\bf u}_0 = [u_0, u_0, \dots, u_0]^T$ and
$Q\otimes f$ is the non-linear operator evaluating the matrix vector product between $Q$ and
$[f(u(\tau_1), \tau_1), [f(u(\tau_2), \tau_2),\dots,[f(u(\tau_M), \tau_M)]^T$.

There is a link between Runge-Kutta method and collocation method in the unconscious scientific culture of applied mathematics, but it seems to have fadded over the years as the people who found it too easy to write it down died progressively. So here is the main question :

> **Can we write any Runge-Kutta method into the same $Q$ matrix formulation as collocation methods ? And if yes (which is very likely), how ?**

From this can be derived additional research topics that could merge technics used for collocation methods and Runge-Kutta ones, for instance :

- _Can we apply a Picard iteration for a Runge-Kutta method, and how stable is it in comparison to the same Runge-Kutta time-stepping ?_
- _Can we develop SDC based iterations using a Runge-Kutta method instead of collocation method for the Q matrix ? And can we use the idea of diagonal optimized SDC on it to develop parallel-in-time RK methods ?_
- _[And probably many others ...]_

## Repository structure

- [grib](./grib/README.md) : all the notes, reports, articles etc ... on the main and derived topics.
- [python](./python) : python code ([library](./python/code/README.md), [scripts](./python/scripts/README.md) and [notebooks](./python/notebook/README.md)).

## Starting points

1.