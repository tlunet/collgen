# From Nodes-to-Nodes (N2N) and Zero-to-Nodes (Z2N) formulation

## Generic Nodes-to-Nodes formulation

Let's consider one time-step of a collocation method on the following nodes $\tau_1, \tau_2, \dots, \tau_M$ to solve

$$
\frac{du}{dt} = \lambda u, \quad t \in [t_0, t_1], \quad u(t_0) = u_0,
$$

We note ${\bf u}$ the nodes solution vector as usuall,
and $\Delta{\tau}_m = \tau_m-\tau_{m-1}, \quad \tau_0 = 0$.

If we apply a given Runge-Kutta method with stability function
$R$ to compute the solution from nodes to nodes, then ${\bf u}$
is the solution of the system :

$$
\begin{pmatrix}
g_1^{-1} & & & \\
-1 & g_2^{-1} & \\
& \ddots & \ddots
\end{pmatrix}
\begin{pmatrix}
u_1 \\ u_2 \\ \vdots
\end{pmatrix}
=
\begin{pmatrix}
u_0 \\ 0 \\ \vdots
\end{pmatrix}
= u_0{\bf e}_1,
$$

with $g_m :=  R(\lambda\Delta{\tau}_m)$ the amplification factor associated to the time-integration between $\tau_{m-1}$ and $\tau_m$.

Here the right-hand side vector is sparse and the matrix translate the computation from one node to an other, hence we denote this the **Nodes-to-Nodes formulation (N2N)**.

This differs from the classical way of writing collocation method (_e.g_ in pySDC) using a dense right-hand side vector with $u_0$
for each element, _i.e_ $[u_0, u_0, \dots, u_0]=u_0{\bf e}$.
We denote that case by **Zero-to-Nodes formulation (Z2N)**.

## From Nodes-to-Nodes to Zero-to-Nodes formulation

### Backward Euler

Replacing the stability function, we have :

$$
\begin{pmatrix}
1-\lambda\Delta{\tau}_1 & & & \\
-1 & 1-\lambda\Delta{\tau}_2 & \\
& \ddots & \ddots
\end{pmatrix}
\begin{pmatrix}
u_1 \\ u_2 \\ \vdots
\end{pmatrix}
=
\begin{pmatrix}
u_0 \\ 0 \\ \vdots
\end{pmatrix}
$$

To get back the Z2N formulation, we add the first line of the linear system to the second, then the newly second line the the third, and so on ... which is equivalent to multiplying both side on the left by

$$
{\bf T} = \begin{pmatrix}
1 & & &\\
1 & 1 & &\\
1 & 1& 1 &\\
\vdots & \vdots & & \ddots
\end{pmatrix}.
$$

This yields the Z2N form of the problem

$$
\begin{pmatrix}
1-\lambda\Delta{\tau}_1 & & & \\
-\lambda\Delta{\tau}_1 & 1-\lambda\Delta{\tau}_2 & \\
\vdots & \ddots & \ddots
\end{pmatrix}
\begin{pmatrix}
u_1 \\ u_2 \\ \vdots
\end{pmatrix}
=
\begin{pmatrix}
u_0 \\ u_0 \\ \vdots
\end{pmatrix},
$$

that we can write equivalently

$$
(I-\lambda\Delta{t}Q_{\Delta,BE}) {\bf u} = u_0{\bf e}
$$

### Forward Euler

Replacing the stability function, we have :

$$
\begin{pmatrix}
(1+\lambda\Delta{\tau}_1)^{-1} & & & \\
-1 & (1+\lambda\Delta{\tau}_2)^{-1} & \\
& \ddots & \ddots
\end{pmatrix}
\begin{pmatrix}
u_1 \\ u_2 \\ \vdots
\end{pmatrix}
=
\begin{pmatrix}
u_0 \\ 0 \\ \vdots
\end{pmatrix},
$$

then multiplying each side by

$$
\begin{pmatrix}
1+\lambda\Delta{\tau}_1 & & & \\
& 1+\lambda\Delta{\tau}_2 & \\
& & \ddots
\end{pmatrix}
$$

we get

$$
\begin{pmatrix}
1 & & & \\
-(1+\lambda\Delta{\tau}_1) & 1 & \\
& \ddots & \ddots
\end{pmatrix}
\begin{pmatrix}
u_1 \\ u_2 \\ \vdots
\end{pmatrix}
= (1+\lambda\Delta{\tau}_1)
\begin{pmatrix}
u_0 \\ 0 \\ \vdots
\end{pmatrix},
$$

and finally switching to Z2N formulation (_i.e_ multiplying by $T$) yields

$$
\begin{pmatrix}
1 & & & \\
-\lambda\Delta{\tau}_1 & 1 & \\
\vdots & \ddots & \ddots
\end{pmatrix}
\begin{pmatrix}
u_1 \\ u_2 \\ \vdots
\end{pmatrix}
= (1+\lambda\Delta{\tau}_1)
\begin{pmatrix}
u_0 \\ u_0 \\ \vdots
\end{pmatrix}.
$$

Again, this can be written equivalently

$$
(I-\lambda\Delta{t} Q_{\Delta,FE}){\bf u} = (1+\lambda\Delta{\tau}_1) u_0{\bf e}.
$$

This is equivalent as the default Z2N formulation for the collocation methods, except this additional factor on the right-hand side

> :scroll: Note that for Lobatto or Radau-Left nodes, this term can be ignored since it's cancelled.
> But in general, any Runge-Kutta method that have a non identity numerator in its stability function will have to deal with this term
> (_e.g_ trapezoidal rule has a $\Delta{\tau}_1/2$, etc ...).

## From Zero-to-Nodes to Nodes-to-Nodes formulation

Since the $T$ matrix defined above is invertible, we can simply transform a collocation problem written in Z2N form into N2N form
by multiplying on the left with

$$
T^{-1} =
\begin{pmatrix}
1 & & &\\
-1 & 1 & &\\
 & -1 & 1 &\\
 & & \ddots & \ddots
\end{pmatrix},
$$

that is :

$$
T^{-1} (I - \lambda \Delta{t}Q){\bf u} = u_0 T^{-1}{\bf e}
= u_0 {\bf e}_1.
$$

which gives us this Z2N $Q$-matrix form for the collocation problem :

$$
(I - \lambda \Delta{t}T^{-1}Q - I_L){\bf u}
= u_0 {\bf e}_1,
$$

where $I_L$ is a matrix containing only ones on its first lower diagonal.