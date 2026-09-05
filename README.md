# Noisy Resonance Analysis

Numerical simulation and signal analysis of a **driven damped harmonic oscillator under Gaussian noise**, implemented in Python.

The project simulates a noisy resonance-spectroscopy experiment, solves the oscillator dynamics using multiple numerical integration methods, analyzes the resulting signal in the frequency domain, estimates physical parameters through model fitting, and evaluates their uncertainty using Monte Carlo simulations.

## Overview

The simulated system is a driven damped harmonic oscillator described by

\[
x''(t) + 2\gamma x'(t) + \omega_0^2 x(t)
=
\frac{F_0}{m}\cos(\Omega t) + \sigma \eta(t),
\]

where:

- \(\omega_0\) is the natural angular frequency
- \(\gamma\) is the damping rate
- \(\Omega\) is the driving frequency
- \(\sigma\) controls the noise strength
- \(\eta(t)\) represents Gaussian white noise

The goal is to recover physical information about the resonator from the simulated noisy time-domain signal.

## Analysis Pipeline

The project is divided into five stages.

### 1. Gaussian Noise Generation

A reproducible Gaussian white-noise signal is generated and statistically validated using:

- Sample mean
- Standard deviation
- Mean-to-standard-deviation ratio
- Noise histogram
- Time-domain visualization

### 2. Numerical ODE Integration

The second-order oscillator equation is rewritten as a first-order system and solved numerically using two methods:

- **Forward Euler**
- **Adaptive RK45**

The implementations are compared in terms of:

- Resulting displacement \(x(t)\)
- Number of function evaluations
- Execution time

This provides a direct comparison between a simple fixed-step numerical method and an adaptive higher-order solver.

### 3. Frequency-Domain Analysis

The RK45 displacement signal is transformed into the frequency domain using the **Fast Fourier Transform (FFT)**.

A one-sided **Power Spectral Density (PSD)** is calculated from the signal in order to identify its dominant frequency components.

The resonance frequency is estimated by locating the main spectral peak in the relevant frequency range.

### 4. Resonance Parameter Estimation

The PSD is fitted to a Lorentzian-style resonance model:

\[
M(f) =
\frac{A}
{(f^2-f_0^2)^2 + 4\gamma^2 f^2}
+ B
\]

The model parameters are estimated by minimizing the **Sum of Squared Errors (SSE)** between the measured PSD and the model.

Optimization is performed using the **Nelder-Mead algorithm**.

From the fitted parameters, the oscillator's quality factor is estimated as

\[
Q = \frac{\pi f_0}{\gamma}.
\]

This stage demonstrates how physical parameters can be recovered from noisy frequency-domain data.

### 5. Monte Carlo Uncertainty Analysis

To quantify the effect of noise on the estimated parameters, the full simulation and fitting process is repeated for **50 independent noise realizations**.

For each accepted realization, the analysis estimates:

- Resonance frequency \(f_0\)
- Quality factor \(Q\)

The resulting distributions are used to compute the mean and standard deviation of the estimates.

Non-physical fitting results are filtered using validity checks on parameters such as the damping rate, resonance frequency, and quality factor.

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
