# Noisy Resonance Analysis

Numerical simulation and signal analysis of a **driven damped harmonic oscillator under Gaussian noise**, implemented in Python.

The project simulates a noisy resonance-spectroscopy experiment, solves the oscillator dynamics using different numerical integration methods, analyzes the resulting signal in the frequency domain, estimates physical parameters through numerical fitting, and evaluates their uncertainty using Monte Carlo simulations.

## Overview

The simulated system is a driven damped harmonic oscillator described by

$$
\frac{d^2x}{dt^2}
+ 2\gamma \frac{dx}{dt}
+ \omega_0^2 x(t)
=
\frac{F_0}{m}\cos(\Omega t)
+ \sigma \eta(t)
$$

where:

- $\omega_0$ is the natural angular frequency
- $\gamma$ is the damping rate
- $\Omega$ is the driving angular frequency
- $\frac{F_0}{m}$ is the driving-force amplitude per unit mass
- $\sigma$ controls the noise strength
- $\eta(t)$ is Gaussian white noise

The goal of the project is to simulate the oscillator's noisy displacement and recover physical information about the system from the resulting time-series data.

## Analysis Pipeline

### 1. Gaussian Noise Generation

A reproducible Gaussian white-noise signal is generated using NumPy.

The noise is statistically examined using:

- Sample mean
- Standard deviation
- Mean-to-standard-deviation ratio
- Time-domain visualization
- Histogram of noise samples

This verifies that the generated signal behaves approximately as zero-mean Gaussian noise.

### 2. Numerical ODE Integration

The second-order differential equation is rewritten as a first-order system:

$$
\frac{dx}{dt} = v
$$

$$
\frac{dv}{dt}
=
-2\gamma v
-\omega_0^2x
+\frac{F_0}{m}\cos(\Omega t)
+\sigma\eta(t)
$$

The system is solved using two numerical methods:

- **Forward Euler**
- **Adaptive RK45**

The resulting displacement $x(t)$ is compared between the two methods.

The project also compares the methods in terms of:

- Number of function evaluations
- Execution time
- Resulting time-domain trajectory

This provides a practical comparison between a simple fixed-step integration method and an adaptive higher-order solver.

## 3. Frequency-Domain Analysis

The displacement obtained from the RK45 solution is transformed into the frequency domain using the **Fast Fourier Transform (FFT)**.

A one-sided **Power Spectral Density (PSD)** estimate is calculated from the signal.

The PSD is used to identify the oscillator's resonance peak and estimate its resonance frequency $f_0$.

The analysis focuses on the physically relevant frequency range in order to distinguish the resonance from other spectral components, including the external driving frequency.

## 4. Resonance Parameter Estimation

The calculated PSD is fitted to a Lorentzian-style resonance model:

$$
M(f)
=
\frac{A}
{\left(f^2-f_0^2\right)^2
+4\gamma_{\mathrm{fit}}^2f^2}
+B
$$

where:

- $A$ controls the amplitude of the resonance
- $f_0$ is the fitted resonance frequency
- $\gamma_{\mathrm{fit}}$ is the fitted damping parameter
- $B$ represents the spectral background

The model parameters are determined by minimizing the **Sum of Squared Errors (SSE)**:

$$
\mathrm{SSE}
=
\sum_i
\left[
M(f_i)-S_x(f_i)
\right]^2
$$

Optimization is performed using the **Nelder-Mead algorithm** through `scipy.optimize.minimize`.

From the fitted resonance frequency and damping rate, the oscillator's quality factor is estimated as

$$
Q
=
\frac{\omega_0}{2\gamma_{\mathrm{fit}}}
=
\frac{\pi f_0}{\gamma_{\mathrm{fit}}}
$$

This stage demonstrates how physical parameters can be extracted from noisy frequency-domain measurements.

## 5. Monte Carlo Uncertainty Analysis

To study the effect of random noise on the estimated resonance parameters, the simulation and fitting procedure is repeated for **50 noise realizations**.

For each accepted realization, the analysis estimates:

- Resonance frequency $f_0$
- Quality factor $Q$

The mean and standard deviation of the resulting estimates are then calculated.

To avoid including non-physical fitting results, the code performs validity checks on parameters such as:

- Resonance frequency
- Damping rate
- Fit amplitude
- Background level
- Quality factor

The resulting distribution of $Q$ values is visualized using a histogram.

## Project Structure

```text
noisy-resonance-analysis/
├── resonance_analysis.py
├── outputs/
│   ├── A_noise_time.png
│   ├── A_noise_hist.png
│   ├── B_time_series_compare.png
│   ├── C_psd_peak.png
│   ├── D_fit.png
│   └── E_Q_hist.png
└── README.md
