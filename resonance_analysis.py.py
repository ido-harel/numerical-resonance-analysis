# name:Ido Harel
# ID: 314937400

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

## Part A 

# set the parameters according to my ID "314937400"
w0 = 2 * np.pi * 2.8        
gamma = 1/4               
F0_m = 1.4                  
Omega = 2 * np.pi * 0.8     
sigma = 0.35                
noise_seed = 7400           

# A.1 - choosing dt and defind time vector
T = 2 * np.pi / w0
dt = 0.001  
T_samp = 20
t = np.arange(0.0, T_samp + dt, dt)
# next line only used to get the first 0.1s of the white noise for the report
# t_x = t[:100]

# A.2 - Generate noise and save plots as path
N = len(t)
# next line only used to get the first 0.1s of the white noise for the report
# N = len(t_x)
np.random.seed(noise_seed)
eta = np.random.normal(0.0, 1.0, N)

# plots as path
outdir = Path("outputs")
outdir.mkdir(exist_ok = True)

# plot eta(t) 
plt.figure()
# next line only used to get the first 0.1s of the white noise for the report
# plt.plot(t_x, eta, linewidth = 0.8)
plt.plot(t, eta, linewidth = 0.8)
plt.xlabel("t [s]")
plt.ylabel(r"$\eta(t)$")
plt.title("Gaussian white noise")
plt.tight_layout()
plt.savefig(outdir / "A_noise_time", dpi = 300)
plt.close()

# A.3 - calculate and plot histogram
print("\nPart A:")
mu = float(eta.mean())
sd = float(eta.std(ddof=0))
ratio = abs(mu)/sd 
print(f"mean = {mu:.6f}")
print(f"std  = {sd:.6f}")
print(f"mean/std = {ratio:.6f}")

# histogram
plt.figure()
plt.hist(eta, bins = 30)
plt.xlabel(r"$\eta$(t)")
plt.ylabel("probability density")
plt.title("Histogram of noise samples")
plt.tight_layout()
plt.savefig(outdir / "A_noise_hist", dpi = 300)
plt.close()

print("outputs/A_noise_time.png")
print("outputs/A_noise_hist.png")

## Part B
from scipy.integrate import solve_ivp
import time

# we want the values for every eta(t) so the RK45 method will be accuarate,
# thus we need the interpulation helper methd
# we are using it in euler even tho we dont need the interpulation there

# helper 
def eta_t(t_q):
    return float(np.interp(t_q, t, eta))

# defenition of the problem - this is the mathematical modle for the integration methods
def f(t_n, y):
    x, v = y
    dxdt = v
    dvdt = -2.0 * gamma * v - (w0**2) * x + F0_m * np.cos(Omega * t_n) + sigma * eta_t(t_n)
    return np.array([dxdt, dvdt], dtype = float)

# initial conditions
x0 = 0.0
v0 = 0.0

# B.2 - Integration methods
# Euler 
y_e = np.zeros((N, 2), dtype = float)
y_e[0] = [x0, v0]

# we want to check if there are time differences therefore we use "time.perf_counter" method
start = time.perf_counter()
for n in range(N - 1):
    y_e[n + 1] = y_e[n] + dt * f(t[n], y_e[n])
e_runtime = time.perf_counter() - start
e_cnt = N - 1
x_e = y_e[:, 0]

# RK45 
start = time.perf_counter()
# acccording to page 17 in rec 11 we saw the function below 
sol = solve_ivp(
    fun = f,
    t_span = (t[0], t[-1]),
    y0 = [x0, v0],
    method = "RK45",
    t_eval = t,      
    rtol = 1e-6,
    atol = 1e-8
)
rk_runtime = time.perf_counter() - start

x_rk = sol.y[0]
rk_cnt = sol.nfev

# B.3 - plot x(t) for both methods 
plt.figure()
plt.plot(t, x_e, linewidth=1, label="Euler")
plt.plot(t, x_rk, linewidth=1, label="RK45")
plt.xlabel("t [s]")
plt.ylabel("x(t)")
plt.title("x(t) comparison: Euler vs RK45")
plt.legend()
plt.tight_layout()
plt.savefig(outdir / "B_time_series_compare", dpi = 300)
plt.close()

# B.4 - compare function evaluations and runtime
print("\nPart B:")
print(f"Euler: nfev = {e_cnt}, runtime = {e_runtime:.4f} s")
print(f"RK45 : nfev  = {rk_cnt}, runtime = {rk_runtime:.4f} s")
print("outputs/B_time_series_compare.png")

## Part C  
x_fft = x_rk 
t_fft = t
N_fft = len(x_fft)

# C.1 - compute FFT - based PSD
# Sampling frequency
fs = 1.0 / dt

# FFT 
X = np.fft.rfft(x_fft)
freq = np.fft.rfftfreq(N_fft, d = dt)

# one sided PSD estimate 
S_x = (1.0 / (fs * N_fft)) * (np.abs(X) ** 2)

# the rfft method takes only the positive frequencies thus we have half of the total power 
# so we need to nultiply by 2 the bins that are not doubled
if N_fft % 2 == 0:
    S_x[1:-1] *= 2.0
else:
    S_x[1:] *= 2.0

# C.2 - identify resonance peak f0 
# as I explained in the report I had to change the search range after my first attempt.
# the code for the first attempt
# idx_peak = np.argmax(S_x[1:]) + 1
# f0 = freq[idx_peak]
# print(f"\nPart C: resonance peak, estimate f0 = {f0:.4f} Hz")

# the fixed code
freq_range = freq > 1
idx_peak = np.argmax(S_x[freq_range])
f0 = freq[freq_range][idx_peak]
print("\nPart C:")
print(f"resonance peak, estimate f0 ≈ {f0:.4f} Hz")

# Plot PSD and mark peaks
plt.figure()
plt.plot(freq, S_x, linewidth=1)
plt.axvline(f0, linestyle = "--", linewidth = 1, label = f"f0 = {f0:.3f} Hz")
plt.xlim(0, 5)  # I cut the x axil to get a better look on the peak
plt.xlabel("f [Hz]")
plt.ylabel(r"$S_x(f)$  [x$^2$/Hz]")
plt.title("Power Spectral Density")
plt.legend()
plt.tight_layout()
plt.savefig(outdir / "C_psd_peak", dpi = 300)
plt.close()

print("outputs/C_psd_peak.png")

## Part D
from scipy.optimize import minimize

# D.1 - modle functions and initial guesses
# lorenzian func
def lorentzian(f, A, f0, gamma_fit, B):
    return A / ((f**2 - f0**2)**2 + (4 * gamma_fit**2 * f**2)) + B

# SSE func
def sse(params, f_data, s_data):
    A, f0, gamma_fit, B = params
    model_vals = lorentzian(f_data, A, f0, gamma_fit, B)
    return np.sum((model_vals - s_data)**2)

# again we select a range to ignore the drive frequency
mask = (freq > 1) & (freq < 5)
f_fit = freq[mask]
s_fit = S_x[mask]

# initial guesses
# I'll explain in the report how did I choose them
f0 = 2.750
gamma_fit = 1/4
B = np.mean(S_x[-50:]) 
peak = np.max(s_fit)
guess_A = (peak - B) * (2 * gamma_fit * f0)**2
initial_params = [guess_A, f0, gamma_fit, B]

# D.2 - minimization
# we use Nelder-Mead method (explained in the report why to choose that nethod)
res = minimize(sse, initial_params, args = (f_fit, s_fit), method='Nelder-Mead')
A_final, f0_final, gamma_final, B_final = res.x

# D.3 - quality factor Q
Q_factor = (np.pi * f0_final) / gamma_final

print(f"\nPart D")
print(f"Estimated f0: {f0_final:.4f} Hz")
print(f"Estimated gamma: {gamma_final:.4f}")
print(f"Quality Factor Q: {Q_factor:.2f}")
print(f"Minimization Success: {res.success}")

# D.4 - plot fit 
plt.figure(figsize = (9, 6))
plt.scatter(freq, S_x, s=2, color='red', alpha=0.5, label='Data')
plt.plot(f_fit, lorentzian(f_fit, *res.x), linewidth=1.5, label="Lorentzian Fit")
plt.axvline(f0_final, linestyle="--", linewidth=1, label=f"f0={f0_final:.2f} Hz")
plt.xlim(0, 5)
plt.ylim(0, peak * 1.4)
plt.xlabel("f [Hz]")
plt.ylabel(r"$S_x(f)$")
plt.title(f"Resonance Fit: f0={f0_final:.2f} Hz, Q={Q_factor:.2f}")
plt.legend()
plt.tight_layout()
plt.savefig(outdir / "D_fit.png", dpi = 300)
plt.close()

print("outputs/D_fit.png")

## Part E 
n = 50
# easy to add to a list and then convert to arrays
f0_est = []
Q_est = []
np.random.seed(noise_seed)

# at first when i run the code for part E, 
# I got a problem in the Q values they where negative and some were huge
# that is not fisicall because it means i got gamma < 0
# so i bounded and added terms to avoid this problem
# also I bounded the frequencies as I did in the previouse parts
# and count the accepted and rejected values to be in control and not reject all exept from few
gamma_min = 1e-6          
f0_min, f0_max = 1.0, 5.0 
A_min = 0.0
B_min = 0.0
Q_max = 800              

accepted = 0
rejected = 0

for i in range(n):
    eta_mc = np.random.normal(0.0, 1.0, N)

    def f_mc(t_n, y):
        x, v = y
        curr_eta = np.interp(t_n, t, eta_mc)
        dvdt = -2.0 * gamma * v - (w0**2) * x + F0_m * np.cos(Omega * t_n) + sigma * curr_eta
        return [v, dvdt]

    sol_mc = solve_ivp(f_mc, (t[0], t[-1]), [x0, v0], t_eval = t, method="RK45")
    x_ss = sol_mc.y[0]

    X_mc = np.fft.rfft(x_ss)
    freq_mc = np.fft.rfftfreq(N, d = dt)
    S_mc = (1.0 / (fs * N)) * (np.abs(X_mc) ** 2)

    # correct one-sided doublin, as I did in part C
    if N % 2 == 0:
        S_mc[1:-1] *= 2.0
    else:
        S_mc[1:] *= 2.0

    S_fit_mc = S_mc[mask]

    # Nelder–Mead fit (unconstrained), but we enforce "bounds" by rejecting bad outputs
    res_mc = minimize(
        sse,
        x0=[A_final, f0_final, gamma_final, B_final],
        args=(f_fit, S_fit_mc),
        method='Nelder-Mead'
    )

    if not res_mc.success:
        rejected += 1
        continue

    A_mc, f0_mc, gamma_mc, B_mc = res_mc.x

    # bound checks for gamma
    if (A_mc <= A_min) or (B_mc < B_min) or (f0_mc <= 0) or (gamma_mc <= 0):
        rejected += 1
        continue

    if (gamma_mc < gamma_min) or (f0_mc < f0_min) or (f0_mc > f0_max):
        rejected += 1
        continue

    Q_mc = (np.pi * f0_mc) / gamma_mc
    
    # bound check for Q
    if (not np.isfinite(Q_mc)) or (Q_mc <= 0) or (Q_mc > Q_max):
        rejected += 1
        continue

    # accept
    f0_est.append(f0_mc)
    Q_est.append(Q_mc)
    accepted += 1

f0_est = np.array(f0_est)
Q_est = np.array(Q_est)

print("\nPart E:")
print(f"Accepted runs: {accepted}/{n}")


print("Monte Carlo Results:")
print(f"f0: {np.mean(f0_est):.4f} ± {np.std(f0_est, ddof=1):.4f} Hz")
print(f"Q:  {np.mean(Q_est):.2f} ± {np.std(Q_est, ddof=1):.2f}")

# histogram
plt.figure()
plt.hist(Q_est, bins=15, color='lightgreen', edgecolor='black')
plt.axvline(np.mean(Q_est), color='red', linestyle='dashed', label='Mean Q')
plt.xlabel("Quality Factor Q")
plt.ylabel("Frequency")
plt.title(f"Histogram of Q Estimates ({accepted} / {n} accepted runs)")
plt.legend()
plt.tight_layout()
plt.savefig(outdir / "E_Q_hist.png", dpi = 300)
plt.close()

print("outputs/E_Q_hist.png")
