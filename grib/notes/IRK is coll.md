# The generic implicit Runge-Kutta formulation
Let 

$$
y_t = f(t, y), y_0 = y(0),
$$

where $f$ is continuously differentiable. Let $t_n = n\Delta t$ be the discretization of the time domain and $y_n \approx y(t_n)$.
A generic $s$-stage Runge-Kutta method is then defined as 

$$
k_i = f(t_n + c_i \Delta t, y_n + \Delta t \sum_{j = 1}^s a_{ij} k_j)
$$

and a prolongation equation

$$
y_{n+1} = y_n + \Delta t \sum_{i = 1}^s b_i k_i.
$$

The vectors values $b_i, c_i$ and $a_{ij}$ are formally stored in the Butcher tableau as:

$$
\begin{array}{c|cccc}
c_1 & a_{11} & a_{12} & \dots & a_{1s} \\
c_2 & a_{21} & a_{22} & \dots & a_{2s} \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
c_s & a_{s1} & a_{s2} & \dots & a_{ss} \\
\hline
& b_1 & b_2 & \dots & b_s
\end{array},
$$

or, written compactly as

$$
\begin{array}{c|c}
{\bf c} & {\bf A} \\
\hline
 & {\bf b}^T
\end{array}.
$$

## Change of variables

Now, let us define 

$$
u_i = y_n + \Delta t \sum_{j = 1}^s a_{ij} k_j,
$$

and the definitions of $k_i$ now read

$$
k_i = f(t_n + c_i \Delta t, u_i).
$$

All together, written line by line in a matrix formulation, we have

$$
 \begin{bmatrix}
 u_1 \\
 u_2 \\
 \vdots \\
 u_s 
 \end{bmatrix}=
 \begin{bmatrix}
 y_n \\
 y_n \\
 \vdots \\
 y_n
 \end{bmatrix} + \Delta t {\bf A}
 \begin{bmatrix}
     f(\tau_1, u_1) \\
     f(\tau_2, u_2) \\
     \vdots \\
     f(\tau_s, u_s)
 \end{bmatrix},
$$

where we define $\tau_i = t_n + c_i \Delta t$. If we look more closely, that is exactly the formulation of the collocation problem, 
where $y_n$ denotes the approximation of the solution in point $t_n$ and we are solving on $t_n + \tau_i$ and later on, for $t_{n+1}$ with 
the prolongation equation.
The prolongation equation now looks like 

$$
\begin{align*}
    y_{n+1} &= y_n + \Delta t \sum_{i = 1}^s b_i k_i\\
            &= y_n + \Delta t \sum_{i = 1}^s b_i f(\tau_i, u_i).
\end{align*}
$$

TODO: we have to prove that we can do this substitution $k \rightarrow u$. The way I saw it somewhere is that they 
prove that the collocation problem $F({\bf u}) = {\bf u}$ is a contraction (we have to find on which compact set) and then by the Banach 
fixed point theorem, we get that a solution exists. I guess that is enough? Or maybe with the inverse function theorem? 
(That's why I wrote $f$ has to be continuously differentiable). This part is a bit in the air... source: 
https://www.epfl.ch/labs/anchp/wp-content/uploads/2018/05/part2-1.pdf

# The collocation matrix

Now, let ${\bf Q}$ be a collocation matrix, meaning 

$$
[Q]_{ij} := \int_{t_n}^{\tau_i} c_j(x)dx,
$$ 

where $c_j$ are the Lagrange interpolation polynomials defined in points $\tau_1 < \dots < \tau_s$. Let us now examine what does 
$\bf Q$ do to a vector $(\tau_1^m, \dots, \tau_s^m)$, for $0 \leq m \leq s-1$. We have

$$
{\bf Q}
\begin{bmatrix}
    \tau_1^m \\
    \tau_2^m \\
    \vdots \\
    \tau_s^m
\end{bmatrix}= 
\begin{bmatrix}
    \int_{t_n}^{\tau_1} \sum_j \tau_j^m c_j(x)dx \\
    \int_{t_n}^{\tau_2} \sum_j \tau_j^m c_j(x)dx \\
    \vdots \\
    \int_{t_n}^{\tau_s} \sum_j \tau_j^m c_j(x)dx
\end{bmatrix}=
\begin{bmatrix}
    \int_{t_n}^{\tau_1} p^{(m)}(x)dx \\
    \int_{t_n}^{\tau_2} p^{(m)}(x)dx \\
    \vdots \\
    \int_{t_n}^{\tau_s} p^{(m)}(x)dx
\end{bmatrix},
$$

where $p^{(m)}$ is a polynomial of degree $s-1$ (because the Lagrange polynomials are of degree $s-1$) going through 
points $(\tau_1, \tau_1^m), \dots, (\tau_s, \tau_s^m)$. We can now conclude that $p^{(m)}(x) = x^m$ is a unique solution. 
Now, for $t_n = 0$, the integrals can be computed as

$$
\int_{t_n}^{\tau_i} p^{(m)}(x)dx = \int_{0}^{\tau_i} x^m dx = \frac{\tau_i^{m+1}}{m+1}.
$$

If we want the matrix $\bf A$ to be a collocation matrix, it is evident that the upper condition must be satisfied.
Another condition is that $\tau_i$ have to be distinct, otherwise we cannot form the polynomial Lagrange basis.
This proves the theorem 1.1.1 from the 'Collocation Methods for Volterra Integral
and Related Functional Equations' which is expressing when can the IRK method be obtained by collocation.
