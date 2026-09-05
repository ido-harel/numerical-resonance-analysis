# Noisy Resonance Analysis

Numerical simulation and analysis of a **driven damped harmonic oscillator under Gaussian noise**, implemented in Python.

The project combines numerical integration, frequency-domain analysis, parameter estimation, and Monte Carlo uncertainty analysis.

## Model

The simulated oscillator is described by:

```math
\frac{d^2x}{dt^2}
+ 2\gamma \frac{dx}{dt}
+ \omega_0^2 x(t)
=
\frac{F_0}{m}\cos(\Omega t)
+ \sigma \eta(t)
```

## Numerical Integration

The second-order ODE is converted into a first-order system and solved using:

- Forward Euler
- Adaptive RK45

The two methods are compared in terms of numerical results, function evaluations, and runtime.

## Frequency-Domain Analysis

The simulated displacement signal is analyzed using the **Fast Fourier Transform (FFT)**.

A one-sided **Power Spectral Density (PSD)** is calculated and used to estimate the resonance frequency.

## Parameter Estimation

The PSD is fitted to a Lorentzian-style model:

```math
M(f)
=
\frac{A}
{\left(f^2-f_0^2\right)^2
+4\gamma_{\mathrm{fit}}^2 f^2}
+B
```

The parameters are estimated by minimizing the Sum of Squared Errors:

```math
\mathrm{SSE}
=
\sum_i
\left[
M(f_i)-S_x(f_i)
\right]^2
```

Optimization is performed using the **Nelder-Mead algorithm**.

The quality factor is then calculated as:

```math
Q
=
\frac{\omega_0}{2\gamma_{\mathrm{fit}}}
=
\frac{\pi f_0}{\gamma_{\mathrm{fit}}}
```

## Monte Carlo Analysis

The full simulation and fitting process is repeated for **50 noise realizations** to estimate the uncertainty in the recovered resonance frequency and quality factor.

Non-physical fitting results are filtered before calculating the final statistics.

## Technologies

- Python
- NumPy
- SciPy
- Matplotlib

## Methods Used

- Numerical ODE integration
- Euler and RK45 methods
- FFT and Power Spectral Density
- Numerical optimization
- Lorentzian fitting
- Monte Carlo simulation
- Statistical uncertainty estimation

## Background

Developed as the final project for the **Numerical Methods for Physics** course at Tel Aviv University.
