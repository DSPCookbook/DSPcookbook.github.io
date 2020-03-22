---
layout: recipe
title: Time domain Wiener filter
modified: 2020-3-22
excerpt: Wiener filter
categories: [Optimal-Filtering]
tags: [Statistical Signal Processing, Optimal Filtering]

---

[//]: # "Comment"

## 1) AR(1) process in white Gaussian noise

The objective for the Wiener filter is to retrieve a desired signal $s(n)$ in noise $v(n)$ given the observed signal $x(n)$ i.e.:

$x(n) = s(n) + v(n).$

The time domain Wiener filter is a linear minimum mean square error estimator (LMMSEE) with optimal filter coefficients $\textbf{w}^\ast$ that satisfy following optimization problem:

$ \textbf{w}^\ast = \underset{\textbf{w}}{\text{arg min}} \quad \mathbb{E} \left[ \left( s(n) - \textbf{w}^T \textbf{x}(n) \right)^2 \right] $,

where $\textbf{w} = [w_0,...,w_{M-1}]^T \in \mathbb{R}^{M}$ and $\textbf{x}(n) = [x(n),x(n-1)...,x(n-M+1)]^T \in \mathbb{R}^{M}$ such that the LMMSE estimate of $s(n)$ given $\textbf{x}(n)$ is:

$\hat{s}(n) = \textbf{w}^{\ast T}\textbf{x}(n)$.



## 2) Assumptions


| Assumption   					 			|| Mathematical description   			| 
|-------------------------------------------||:------------------------------------:|
| $s(n)$ and $v(n)$ are uncorrelated. 		|| $\mathbb{E}[s(n)v(n)] = 0$ 			|
| $v(n)$ is a white Gaussian noise. 		|| $v(n) \sim \mathcal{N}(0,\sigma_v^2)$|
| $s(n)$ is an AR(1) process. 		 		|| $s(n) = \alpha s(n-1) + u(n)$ 		|
| $u(n)$ is a white Gaussian process. 		|| $u(n) \sim \mathcal{N}(0,\sigma_u^2)$|
| The noise variance $\sigma_v^2$. 			|| _Known_ 								|
| The variance of $u(n)$ i.e. $\sigma_u^2$. || _Known_ 								|
| The time constant $\alpha$. 				|| _Known_ 								|




## 3) Solution

The Wiener filter is given as

$\textbf{w}^\ast = R_{XX}^{-1} r_{xs}$,

where $R_{XX} \in \mathbb{R}^{M\times M}$ is the autocorrelation matrix of $x(n)$ and $r_{xs} \in \mathbb{R}^M$ is the cross-correlation vector between $s(n)$ and $x(n)$.


## 4) Step-by-step guide

1. Set filter length e.g. $M=100$.
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

	$\text{w}^{\ast} = R_{xx}^{-1}r_{xs}$

8. Perform filtering

	$\hat{s}(n) = \text{w}^{\ast T} \textbf{x}(n)$


## 5) MATLAB code

A MATLAB example is provided where the signals are artificially generated.

```matlab
%   Copyright 2020: DSPCookbook
clc, clear, close all

% Constants in simuation
N       = 1000;
varu    = 1;
varv    = 1000;
a       = 0.999;

% Generate the noisy observed signal
u       = sqrt(varu)*randn(N,1);
s       = 0; % Initialize s to be the mean
for n = 2:N
    s(n,1) = a*s(n-1,1) + u(n);
end
v       = sqrt(varv)*randn(N,1);
x       = s + v;

% Step 1
M       = 100;
tau     = (0:(M-1))';

% Step 2
rvv     = zeros(M,1);
rvv(1)  = varv;             % Autocorrelation is a Kronecker Delta function

% Step 3
rss     = (a.^abs(tau))/(1-a^2)*varu;

% Step 4
rxx     = rss + rvv;        % Autocorrelation vector

% Step 5
Rxx     = toeplitz(rxx);    % Autocorrelation matrix

% Step 6
rxs     = rss;

% Step 7
w       = inv(Rxx)*rxs;

% Step 8
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


------

## Derivation of the Wiener filter

For the derivation of the time domain Wiener filter, check out the [extra material]({{ site.baseurl }}{% link _posts/Optimal-Filtering/2020-03-20-wiener-filter-derivation-extra.md %}).



## C code implementation of the Wiener filter

A C-code implementation of the simulation can be found [here]({{ site.baseurl }}{% link _posts/Optimal-Filtering/2020-03-20-wiener-filter-c-code-extra.md %})





