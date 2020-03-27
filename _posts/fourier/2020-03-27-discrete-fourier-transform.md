---
layout: recipe
title: Discrete Fourier transform (DFT)
modified: 2020-03-27
excerpt: Something with sinusoids.
categories: [Fourier]
tags: [Fourier, Spectral Analysis, DFT]
---


## 1) Introduction

The discrete Fourier transform is often used for spectral analysis of sequences $x(n)$ of finite length e.g. discrete time series under the assumption that $x(n)$ is $N$-periodic ($N$-periodic refers to that $x(n)$ repeats itself for every $N$ sample). The discrete Fourier transform is defined as

$$ X(k) = \sum_{n=0}^{N-1} x(n)\exp\left( - \frac{i 2\pi}{N} kn \right),  $$

where $k$ is the frequency bin index, $n$ is the time index, $N$ is the length of the sequence, and $\sqrt{-1} = i$.

## 2) Assumptions


| Assumption   					 				|| Mathematical description   			| 
|-----------------------------------------------||:------------------------------------:|
| 1. $x(n)$ is $N$-periodic.    				|| $x(n+mN) = x(n)$, $m \in \mathbb{Z}$ |
| 2. $X(k)$ is $N$-periodic. 					|| $X(k+mN) = X(k)$, $m \in \mathbb{Z}$	|


## 3) Step-by-step guide: Plotting the spectrum of $x(n)$

1. Select frequency bins to $k = 0,1,...,N-1$.

2. Compute the discrete Fourier transform for all selected frequency bins.

	$$X(k) = \sum\limits_{n=0}^{N-1} x(n)\exp\left( - \frac{i 2\pi}{N} kn \right), \quad k = 0,1,...,N-1 $$

3. Compute the power spectrum.

	$$P_x (k)= |X(k)|^2, \quad k = 0,1,...,N-1 $$

## 4) MATLAB code

A MATLAB example is provided with example signals.

```matlab
%   Copyright 2020: DSPCookbook
clc, clear, close all

% Generate example signal
N       = 2000;
n       = (0:N-1)';
f       = [10 40 80 150];
fs      = 1e3;
x       = sum(cos(2*pi*f.*n/fs),2);

% Step 1
K = 0:(N-1);

% Step 2
for k = K
    X(k+1,1) = sum(x.*exp(-1j*2*pi*k*(0:N-1)'/N));
end

% Step 3
Px = abs(X).^2;

% Plot Spectrum
plot((0:(N-1))*fs/N,Px)
grid, xlabel('Frequency [Hz]'),ylabel('Amplitude')

```


## 5) C code

A C code implementation of the MATLAB simulation is provided.

```c
#include <stdio.h>              
#include <stdlib.h>             
#include <math.h>               
#include <complex.h>         

//----- Defines -------------------------------------------------------------
#define PI         3.14159265   // The value of pi

//===== Main program ========================================================
void main(void)
{
  FILE     *fpx_out, *fpX_out;
  char     xname[] 	= "signal.txt";  
  char     Xname[] 	= "spectrum.txt";  
  fpx_out = fopen(xname, "w");
  fpX_out = fopen(Xname, "w");

  // Set parameters for simulation  
  int      N 		  = 2000;         
  double   f[] 	  = {10,40,80,150};     
  double   fs 	  = 1000;

  // Allocate array for signal
  double   x[N]; 
  double   complex X[N];
  double   Px[N];


  // Generate the signal x(n) to be sum of sinusoids
  for (int n=0; n<N; n++)
  {
  	x[n] = cos(2*PI*f[0]*n/fs)+cos(2*PI*f[1]*n/fs)+cos(2*PI*f[2]*n/fs)+cos(2*PI*f[3]*n/fs);
    fprintf(fpx_out, "%f \n", x[n]);
  }

  // Step 1 + 2 + 3
  for (int k=0; k<N; k++)
  {
    	double complex temp = 0;

      // Step 2
    	for (int n=0; n<N; n++)
    	{
    		temp = temp + x[n]*cexp(-1*I * 2*PI/N*k*n);
    	}
      X[k] = temp;

      // Step 3
      Px[k] = pow(abs(X[k]),2);
	    fprintf(fpX_out, "%f \n", Px[k]);
  }

  fclose(fpx_out);
  fclose(fpX_out);
}
```