---
layout: recipe
title: The Wiener filter
modified: 2020-3-17
excerpt: Something with sinusoids.
categories: [Optimal-Filtering]
tags: [Statistical Signal Processing, Optimal Filtering]

---

## Time Domain Wiener Filter

Generally for the Wiener filter, we are interested in retrieving a desired signal in noise. We consider the following signal model

$x(n) = s(n) + v(n),$

where:
- $x(n)$ is the noisy observed signal.
- $v(n)$ is the noise.
- $s(n)$ is the desired signal.
- $n$ is the time index.

**Overall Objective**: Retrieve $s(n)$ given $x(n)$

The Wiener filter finds the optimal filter coefficients that solves following optimization problem

$ w_0^\ast,...,w_{I-1}^\ast = \underset{w_0,...,w_{I-1}}{\text{arg min}} \quad \mathbb{E} \left[ \left( s(n) - \sum\limits_{i=0}^{I-1} w_i x(n-i) \right)^2 \right] $

or equivalently

$ \textbf{w}^\ast = \underset{\textbf{w}}{\text{arg min}} \quad \mathbb{E} \left[ \left( s(n) - \textbf{w}^T \textbf{x}(n) \right)^2 \right] $

where $\textbf{w} = [w_0,...,w_{I-1}]^T$ and $\textbf{x}(n) = [x(n),x(n-1)...,x(n-I+1)]^T$


### Case 1: AR(1) process in white Gaussian noise

#### Assumptions:
1. $s(n)$ and $v(n)$ are uncorrelated random processes.
2. $v(n)$ is white Gaussian noise $v(n) \sim \mathcal{N}(0,\sigma_v^2)$.
3. $s(n)$ is an AR(1) random process such that
	$s(n) = \alpha s(n-1) + u(n)$
4. $u(n)$ is white Gaussian process $u(n) \sim \mathcal{N}(0,\sigma_u^2)$.
5. _Known_: The noise variance $\sigma_v^2$
6. _Known_: The variance of $u(n)$ i.e. $\sigma_u^2$
7. _Known_: The time constant $\alpha$

#### Solution:

$\textbf{w} = R_{XX}^{-1} r_{xs}$

#### Step-by-step:

1. Set filter length e.g. $M=20$.
2. Set autocorrelation of the noise to

	$r_{vv}(k) = \sigma_v^2 \delta (k), \quad \text{for }k = 0,...,M-1$

	where $\delta (k)$ is the Kronecker delta function and $k$ is the lag-index.
3. Compute the autocorrelation of $s(n)$ for $M$-lags

	$r_{ss}(k) = \frac{\alpha^{k}}{1-\alpha^2} \sigma_s^2, \quad \text{for }k = 0,...,M-1 $

4. Compute the autocorrelation of $x(n)$

	$r_{xx}(k) = r_{ss}(k) + r_{vv}(k)$

5. Form the autocorrelation matrix $R_{xx}$ of $x(n)$
6. Form the cross-correlation vector between $x(n)$ and $s(n)$

	$r_{xs} = r_{ss}$

7. Compute the Wiener filter coefficients

	$\text{w} = R_{xx}^{-1}r_{xs}$


#### MATLAB code

<details>
<summary>
Click here to see the code
</summary>

```matlab
%   Copyright 2020: DSPCookbook
clc, clear, close all

% Number of samples
N       = 1000;

% Statistics of signals
varu    = 1;
varv    = 1000;
a       = 0.999;

% Generate the desired signal
u       = sqrt(varu)*randn(N,1);
s       = 0; % Initialize s to be the mean
for n = 2:N
    s(n,1) = a*s(n-1,1) + u(n);
end

% Generate the noise signal
v       = sqrt(varv)*randn(N,1);

% Generate the noisy observed signal
x       = s + v;


%% Obtaining the Wiener filter

% Desired Wiener filter size
M       = 100;
tau     = (0:(M-1))';

% Autocorrelation vector of noise v(n)
rvv     = zeros(M,1);
rvv(1)  = varv;             % Autocorrelation is a Kronecker Delta function

% Autocorrelation vector of desired signal s(n)
rss     = (a.^abs(tau))/(1-a^2)*varu;

% Autocorrelation of noisy observation x(n)
rxx     = rss + rvv;        % Autocorrelation vector
Rxx     = toeplitz(rxx);    % Autocorrelation matrix

% Cross-correlation vector between x(n) and s(n)
rxs     = rss;

% Compute Wiener filter coefficients
w       = inv(Rxx)*rxs;


%% Perform filtering
for n = M:N
    y(n,1) = w'*flipud(x((n-M+1):n));
end

%% Plot result
plot(x)
grid, xlabel('time index n'), ylabel('Amplitude'), hold on
plot(y)
plot(s)
legend('x(n)','y(n)','s(n)')
title('Wiener filtering')

```

</details>


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

#### What to do on embedded devices?

You have decided that you would like to implement the Wiener filter into an embedded device (let's say a small microcontroller). You are very exciting to get it up an running, but as soon as you start coding and implementing the Wiener filter, you realize that MATLAB is not running the microcontroller, and you are scared of the thought of inverting a matrix that is bigger than 2x2. Of course your matrix is bigger than 2x2, as you had big plans of making a 50-order Wiener filter which makes the matrix you are supposed to invert 50x50. As your motivation and hopes lay in ruins, you are about to throw the Wiener filter into the trashbin and start on something else - but wait! Fortunately for you, you can avoid the matrix inversion by implementing an numerical solver that solves the optimization problem iteratively. It will not get you the exact optimal solution, but given that optimization problem for the LMMSE estimator is convex and quadratic, a numerical solver can be made very efficient and converge to the optimum solution at a high rate.



... More to come!


#### Practical C code implementation of the Wiener filter

This is a sippet of the code. You may download the full code here...

```c
void main(void)
{
  FILE     *fpv_out, *fpy_out,
  *fpu_out, *fps_out, *fpx_out;               		
  char     vname[] 	= "v.txt";  
  char     uname[] 	= "u.txt";
  char     sname[] 	= "s.txt";  
  char     xname[] 	= "x.txt";  
  char     yname[] 	= "y.txt";  
  fpv_out 			= fopen(vname, "w");
  fpu_out 			= fopen(uname, "w");
  fps_out 			= fopen(sname, "w");
  fpx_out 			= fopen(xname, "w");
  fpy_out 			= fopen(yname, "w");

  // Set random seed
  rand_val(1);

  // Set parameters for simulation  
  int      N 		= 2000;         
  double   stdv 	= 10;     
  double   stdu 	= 1;
  double   alpha 	= 0.99;
  int      M 		= 50;
  double   w[M];
  // Allocate array for correlation vector
  double   rss[M];           
  double   rvv[M];
  double   rxx[M];
  double   rxs[M];
  double   Rxx[M][M];

  // Allocate array for signal
  double   v[N];           
  double   u[N];
  double   s[N];
  double   x[N];
  double   y[N];

  // Generate random numbers
  for (int i=0; i<N; i++)
  {
    v[i] = norm(0, stdv);
	u[i] = norm(0, stdu);
    fprintf(fpv_out, "%f \n", v[i]);
    fprintf(fpu_out, "%f \n", u[i]);
  }


  // Generate the desired signal s(n) to be an AR(1) process
  s[0] = 0;
  for (int n=1; n<N; n++)
  {
  	y[n] = 0;
  	s[n] = alpha*s[n-1] + u[n];
  	x[n] = s[n] + v[n];
    fprintf(fps_out, "%f \n", s[n]);
    fprintf(fpx_out, "%f \n", x[n]);
  }

  // Compute autocorrelation of v(n)
  for (int k=0; k<M; k++)
  {
  	w[k] 	= 0;
  	rvv[k] 	= (k == 0) ? pow(stdv,2) : 0;
  }


  for (int k=0; k<M; k++)
  {
  	// Compute autocorrelation of s(n)
  	rss[k] = pow(alpha,k)/(1-pow(alpha,2))*pow(stdu,2);

  	// Cross-correlation of x(n) and s(n)
  	rxs[k] = rss[k];

  	// Compute autocorrelation of s(n)
  	rxx[k] = rss[k] + rvv[k];
  }

  // Autocorrelation matrix of x(n)
  for (int i=0; i<M; i++)
  {
  	for (int j=0; j<M; j++)
  	{
  		Rxx[i][j] = rxx[abs(i-j)];
  		// printf("%f\t",Rxx[i][j] );
  	}
  	// printf("\n");
  }


  int N_iter 	= 20;
  double beta 	= 0.001;
  double temp 	= 0;
  for (int k=0; k<N_iter; k++)
  {
  	for (int i=0; i<M; i++)
	{
	  // Compute the Rxx*w
	  for (int j=0; j<M; j++)
	  {
	  	temp = temp + Rxx[i][j]*w[j];
	  }
	  w[i] = w[i] - beta*(-rxs[i] + temp);
	  temp = 0;
	  //printf("%0.4f\t",w[i]);
	}
	//printf("\n");
  }


  for (int n=M-1; n<N; n++)
  {
  	y[n] = 0;
  	for (int i=0; i<M; i++)
	{
		y[n] = y[n] + w[i]*x[n-i];
	}
	fprintf(fpy_out, "%f \n", y[n]);
  }

  fclose(fpv_out);
  fclose(fpu_out);
  fclose(fps_out);
  fclose(fpx_out);
  fclose(fpy_out);
}

```
