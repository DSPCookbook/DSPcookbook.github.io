---
layout: single
title: Derivation of the Wiener Filter
modified: 2020-3-22
excerpt: Something with sinusoids.
categories: [Optimal-Filtering]
tags: [Statistical Signal Processing, Optimal Filtering]
classes: wide2
author: Poul
---

## Derivation of the time domain Wiener Filter

#### Derivation of the Wiener filter using linear algebra

The Wiener filter is a linear minimum mean-square-error (LMMSE) estimator, thus we seek a linear combination of the sequence $x(n)$ that minimizes the MSE. Let us define the mean-square-error to be

$\text{MSE} \triangleq \mathbb{E} \left[ \left( s(n) - \hat{s}(n) \right)^2 \right]$

where $s(n) - \hat{s}(n)$ is the error signal between the desired signal $s(n)$ and the estimated one $\hat{s}(n)$. Particularly for the Wiener filter is that $\hat{s}(n)$ is a linear combination of the noisy signal $x(n)$ such that $\hat{s}(n) = \textbf{w}^T \textbf{x}(n)$ (hence the name LMMSE). Therefore, the MSE for the Wiener filter is

$\text{MSE} = \mathbb{E} \left[ \left( s(n) - \textbf{w}^T \textbf{x}(n) \right)^2 \right]$.

As we are interested in obtaining $\textbf{w}$, we minimize the MSE and obtain the argument at the minima hence

$ \textbf{w}^\ast = \underset{\textbf{w}}{\text{arg min}} \quad \mathbb{E} \left[ \left( s(n) - \textbf{w}^T \textbf{x}(n) \right)^2 \right] $.

Obtaining the optimum solution $\textbf{w}^\ast$ is fairly straightforward given that we have some fundamental knowledge on multivariate calculus. We essentially only need to take the derivative with respect $\textbf{w}$, set it equal to 0 and solve for $\textbf{w}$ - just like the way you usually would solve an optimization problem. Of course, for the multivariate case i.e. where our solution of $\textbf{w}$ is a vector, we compute the gradient instead taking the derivative with respect to a coefficient in $\textbf{w}$. The gradient of the MSE is

$\nabla_w \text{MSE} = \nabla_w \mathbb{E} \left[ \left( s(n) - \textbf{w}^T \textbf{x}(n) \right)^2 \right] =  \mathbb{E} \left[ \nabla_w \left( s(n) - \textbf{w}^T \textbf{x}(n) \right)^2 \right]$

It is important to note that the expectation operator $\mathbb{E}[\cdot]$ is linear which means that we may move the gradient operator (which is also linear) inside the expectation operator. Taking the gradient of $\left( s(n) - \textbf{w}^T \textbf{x}(n) \right)^2$ is straightforward as we can simply use the chain rule. Using the chain rule, the gradient of the MSE is

$\nabla_w \text{MSE} = \mathbb{E} \left[ -2 \left( s(n) - \textbf{w}^T \textbf{x}(n) \right) \textbf{x}(n)  \right] = 0.$

We multiply into the parentesis and expand to

$\nabla_w \text{MSE} = \mathbb{E} \left[  -2 s(n)\textbf{x}(n) + 2 \textbf{x}(n) \textbf{x}^T(n) \textbf{w}  \right] = 0.$

and note that $\textbf{x}(n) \textbf{x}^T(n) \textbf{w} = (\textbf{w}^T \textbf{x}(n))\textbf{x}(n)$

$\nabla_w \text{MSE} = -\mathbb{E}\left[ s(n)\textbf{x}(n)\right] + \mathbb{E}\left[ \textbf{x}(n) \textbf{x}^T(n) \right] \textbf{w} = 0.$

and here it is also important to note that $\textbf{w}$ is not a random vector which means that we can move $\textbf{w}$ outside the expectation operator. We now have $\mathbb{E}\left[ s(n)\textbf{x}(n)\right]$ which is the cross-correlation _vector_ between $s(n)$ and $\textbf{x}(n)$ and $\mathbb{E}\left[ \textbf{x}(n) \textbf{x}^T(n) \right]$ which is the autocorrelation _matrix_. Let $R_{xx} = \mathbb{E}\left[ \textbf{x}(n) \textbf{x}^T(n) \right]$ and $r_{xs} = \mathbb{E}\left[ s(n)\textbf{x}(n)\right]$ then solving for $\textbf{w}$ becomes

$\text{w} = R_{xx}^{-1}r_{xs}$

Therefore

$\underset{\textbf{w}}{\text{arg min}} \quad \mathbb{E} \left[ \left( s(n) - \textbf{w}^T \textbf{x}(n) \right)^2 \right] = R_{xx}^{-1}r_{xs}$



#### Derivation of the Wiener filter without linear algebra

We recommend you to get comfortable with multivariate calculus and linear algebra as it is by far more convenient that deriving the solution to the Wiener filter without linear algebra. And as you will see, you most likely will end up using linear algebra to get to the final solution anyway.

Nonetheless, the optimization problem we are trying to solve in this case is

$ w_0^\ast,...,w_{I-1}^\ast = \underset{w_0,...,w_{I-1}}{\text{arg min}} \quad \mathbb{E} \left[ \left( s(n) - \sum\limits_{i=0}^{M-1} w_i x(n-i) \right)^2 \right] $

To solve this optimization problem, we take the partial derivative with respect to the Wiener filter coefficients. That is

$\frac{\partial}{\partial w_j} \text{MSE} = \mathbb{E} \left[ \frac{\partial}{\partial w_j} \left( s(n) - \sum\limits_{i=0}^{M-1} w_i x(n-i) \right)^2 \right]$

Again we use the chain rule such that

$\frac{\partial}{\partial w_j} \text{MSE} = \mathbb{E} \left[ -2 \left( s(n) - \sum\limits_{i=0}^{M-1} w_i x(n-i) \right)x(n-j)  \right] = 0$

$ = \mathbb{E} \left[  -s(n) x(n-j) + \sum\limits_{i=0}^{M-1} w_i x(n-i) x(n-j)  \right] $

$ = -\mathbb{E} \left[  s(n) x(n-j) \right] + \sum\limits_{i=0}^{M-1} w_i \mathbb{E} \left[ x(n-i) x(n-j)  \right] $

where $r_{sx_j} = \mathbb{E} \left[  s(n) x(n-j) \right]$ is the cross-correlation between $s(n)$ and $x(n-j)$, and $r_{x_i x_j} = \mathbb{E} \left[ x(n-i) x(n-j)  \right]$ is the autocorrelation between $x(n-i)$ and $x(n-j)$. We now have

$ \sum\limits_{i=0}^{M-1} w_i r_{x_i x_j} = r_{sx_j}$

which is

$ w_0 r_{x_0 x_j} + ... + w_j r_{x_j x_j} + ... + w_{M-1} r_{x_{M-1} x_j} = r_{sx_j}$

and as we can see, the solution of $w_j$ depends on all other coefficients of the Wiener filter also (no surprise here). However, if we take the partial derivative for all $M$ coefficient we end up with $M$ linear equations with $M$ unknowns i.e.

$\frac{\partial}{\partial w_0} \text{MSE} = 0 \Rightarrow w_0 r_{x_0 x_0} + ... + w_{M-1} r_{x_{M-1} x_0} = r_{sx_0}$

$\vdots$

$\frac{\partial}{\partial w_{M-1}} \text{MSE} = 0 \Rightarrow w_0 r_{x_0 x_{M-1}} + ... + w_{M-1} r_{x_{M-1} x_{M-1}} = r_{sx_{M-1}}$.

On this form linear algebra comes in very handy as we may put the equations into matrix-vector form making it easy to solve.
