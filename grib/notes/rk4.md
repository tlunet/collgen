# 4th order Runge Kutta method (the classical one)

## Generic formulation

RK4 is usually written using $k$ coefficients :

$$
u_{n+1} = u_{n} + \frac{\Delta{t}}{6}(k_1 + 2k_2 + 2k_3 + k_4)
$$

with

$$
\begin{align}
k_1 &= f(u_n,t_n) \\
k_2 &= f\left(u_n + \frac{\Delta{t}}{2}k_1,
    t_n+\frac{\Delta{t}}{2}\right) \\
k_3 &= f\left(u_n + \frac{\Delta{t}}{2}k_2,
    t_n+\frac{\Delta{t}}{2}\right) \\
k_4 &= f\left(u_n + \Delta{t}k_3,
    t_n+\Delta{t}\right) \\
\end{align}
$$

but an other representation is also use, by defining intermediate
solutions :

$$
\begin{align}
u_1 &= u_n \\
u_2 &= u_n + \frac{\Delta{t}}{2}f(u_1, t_1) \\
u_3 &= u_n + \frac{\Delta{t}}{2}f(u_2, t_2) \\
u_4 &= u_n + \Delta{t}f(u_3, t_3)
\end{align}
$$

and a final update

$$
u_{n+1} = u_n + \frac{\Delta{t}}{6}[
    f(u_1, t_1) + 2f(u_2, t_2) + 2f(u_3, t_3) + f(u_4, t_4)]
$$

with $t_1 = t_n$, $t_2 = t_3 = t_n + \Delta{t}/2$,
and $t_4 = t_n + \Delta{t}$.

This is actually closer to the Butcher table representation of RK4 :

$$
\begin{array}{c|cccc}
0 & \\
1/2 & 1/2 \\
1/2 & 0 & 1/2 \\
1 & 0 & 0 & 1 \\\hline
  & 1/6 & 1/3 & 1/3 & 1/6
\end{array}
$$

where each line in the upper par correspond to a "node" solution $u_m$, and the lower par corresponds to a "prolongation"

## Q-matrix form

From the Butcher table of RK4, we define the collocation matrix

$$
Q = \begin{pmatrix}
    0 & 0 & 0 & 0 \\
    1/2 & 0 & 0 & 0 \\
    0 & 1/2 & 0 & 0 \\
    0 & 0 & 1 & 0
\end{pmatrix},
$$

and the associated nodes
${\bf \tau}=[0, 1/2, 1/2, 1]^T$ and weights
${\bf w}=[1/6, 1/3, 1/3, 1/6]^T$.

Then we build the collocation problem

$$
(I-\Delta{t}Q\otimes f){\bf u} = {\bf u}_0,
$$

which is written in condensed form :

$$
\left[\begin{matrix}
1 & 0 & 0 & 0 \\
- \frac{\Delta{t} f}{2} & 1 & 0 & 0 \\
0 & - \frac{\Delta{t} f}{2} & 1 & 0 \\
0 & 0 & - \Delta{t} f & 1
\end{matrix}\right]
\begin{pmatrix}
u_1 \\ u_2 \\ u_3 \\ u_4
\end{pmatrix}
=
\begin{pmatrix}
u_0 \\ u_0 \\ u_0 \\ u_0
\end{pmatrix}
$$

