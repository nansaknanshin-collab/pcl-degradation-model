import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
import os

# =========================================================
# DIMENSIONAL MODEL
#   x   in mm
#   t   in hours
#   De0 in mm/h
#   k*  in 1/h
#
# FIXES vs original:
#   FIX 1: L_THICKNESS_MM = 0.25  (half of 0.5 mm film, not 0.2)
#   FIX 2: Porosity uses only polymer solid (C+A), not (C+A+EC+EA)
#           - physically: EC/EA are surface-bound enzyme, not polymer
#           - gravimetric: film is washed before weighing -> bound
#             enzyme is removed
#           - consistent with Wang et al. (2008) whose porosity
#             contains only polymer species
#   FIX 3: WL observable uses only (C+A) as polymer solid
#   FIX 4: Xc observable uses C/(C+A) — unbound polymer only
#   FIX 5: ALPHA_E = 4.5 from Wang et al. (2008) homogenisation
#           (was 1.0, which underestimates porosity-diffusion coupling)
# =========================================================

# ---------------------------
# USER SETTINGS
# ---------------------------
MIN_ERR_PCT = 0.1
USE_WEIGHTS = True

N_GRID   = 25
MAX_NFEV = 800
SEED     = 1

# FIX 1: half-thickness of 0.5 mm film = 0.25 mm
L_THICKNESS_MM = 0.25

# Fixed ICs (dimensionless volume fractions, sum = 1)
C0_POLY = 0.302
A0_POLY = 0.698
E_INIT  = 0.0
P_INIT  = 0.0

# Dirichlet bulk enzyme at x=0 (normalised by bath concentration)
E_BULK  = 1

# FIX 5: porosity-diffusion coupling from Wang et al. (2008)
# homogenisation gives alpha = 4.5 for randomly distributed pores
# below the percolation threshold
ALPHA_E = 1

np.random.seed(SEED)
download_dir = os.getcwd()

# =========================================================
# DATA (time in hours)
# =========================================================
P_raw = np.array([
    [72, 87.56, 0.1],
    [64, 86.12, 3.0],
    [56, 86.06, 1.5],
    [48, 86.03, 1.4],
    [40, 84.80, 1.4],
    [32, 83.75, 1.5],
    [24, 81.47, 1.4],
    [20, 80.47, 3.0],
    [16, 75.67, 1.5],
    [12, 66.66, 1.2],
    [8,  54.85, 1.5],
    [4,  38.12, 1.2],
    [0,   0.00, 0.0]
], dtype=float)

P_raw      = P_raw[np.argsort(P_raw[:, 0])]
t_P_h      = P_raw[:, 0]
WL_obs_pct = P_raw[:, 1]
err_WL_pct = P_raw[:, 2]
sigma_WL   = np.maximum(err_WL_pct, MIN_ERR_PCT)

shi_lipase = np.array([
    [0,  30.2, 0.9],
    [4,  21.6, 0.8],
    [8,  17.1, 0.6],
    [16,  7.2, 0.9],
    [24,  6.9, 0.6],
    [40,  0.0, 0.0],
    [64,  0.0, 0.0]
], dtype=float)

shi_lipase = shi_lipase[np.argsort(shi_lipase[:, 0])]
t_C_h      = shi_lipase[:, 0]
Xc_obs_pct = shi_lipase[:, 1]
err_Xc_pct = shi_lipase[:, 2]
sigma_Xc   = np.maximum(err_Xc_pct, MIN_ERR_PCT)

# Union time grid for one simulation
t_all_h = np.unique(np.concatenate([t_P_h, t_C_h]))

# Spatial grid (mm)
x  = np.linspace(0.0, L_THICKNESS_MM, N_GRID)
dx = x[1] - x[0]

# =========================================================
# VARIABLE-COEFFICIENT DIFFUSION OPERATOR
# d/dx( D(x) du/dx )
# BCs: u(0) = u_left (Dirichlet), du/dx(L) = 0 (Neumann)
# =========================================================
def diffusion_dirichlet_neumann_varD(u, D, u_left):
    u = u.copy()
    u[0] = u_left
    D = np.asarray(D, dtype=float)

    flux = np.zeros(N_GRID + 1, dtype=float)
    D_iface = 0.5 * (D[:-1] + D[1:])
    flux[1:N_GRID] = D_iface * (u[1:] - u[:-1]) / dx
    flux[N_GRID] = 0.0   # Neumann BC

    div = np.zeros_like(u)
    div[1:-1] = (flux[2:N_GRID] - flux[1:N_GRID-1]) / dx
    div[-1] = 2.0 * (flux[N_GRID] - flux[N_GRID-1]) / dx
    div[0] = 0.0   # Dirichlet node
    return div

# =========================================================
# PDE RHS (METHOD OF LINES)
# =========================================================
def rhs_pde(t, y, pars):
    k1, km1, k3, km3, kconv, kdegC, kdegA, De0 = pars

    Y = y.reshape(6, N_GRID)
    E, C, A, EC, EA, P = Y

    E = E.copy()
    E[0] = E_BULK   # Dirichlet BC

    # FIX 2: porosity = fraction of domain NOT occupied by polymer solid
    denom0 = C0_POLY + A0_POLY
    phi = 1.0 - (C + A) / max(denom0, 1e-12)
    phi = np.clip(phi, -0.5, 2.0)

    # FIX 5: ALPHA_E from Wang et al. homogenisation
    De = np.maximum(De0 * (1.0 + ALPHA_E * phi), 1e-30)
    diffE = diffusion_dirichlet_neumann_varD(E, De, E_BULK)

    dE = (diffE
          - k1   * E * C  + km1  * EC
          - k3   * E * A  + km3  * EA
          + kconv  * EC
          + kdegC * EC    + kdegA * EA)

    dC  = -k1  * E * C  + km1  * EC
    dA  = -k3  * E * A  + km3  * EA  + kconv * EC
    dEC =  k1  * E * C  - km1  * EC  - kconv * EC - kdegC * EC
    dEA =  k3  * E * A  - km3  * EA  - kdegA * EA
    dP  =  kdegC * EC   + kdegA * EA   # no diffusion

    dE[0] = 0.0   # hold Dirichlet node fixed

    return np.vstack([dE, dC, dA, dEC, dEA, dP]).reshape(-1)

# =========================================================
# SIMULATION -> OBSERVABLES
# =========================================================
def simulate_observables(times_h, pars):
    E0  = np.full(N_GRID, E_INIT)
    C0v = np.full(N_GRID, C0_POLY)
    A0v = np.full(N_GRID, A0_POLY)
    EC0 = np.zeros(N_GRID)
    EA0 = np.zeros(N_GRID)
    P0  = np.full(N_GRID, P_INIT)

    y0 = np.vstack([E0, C0v, A0v, EC0, EA0, P0]).reshape(6, N_GRID)
    y0[0, 0] = E_BULK
    y0 = y0.reshape(-1)

    times_h = np.asarray(times_h, dtype=float)

    sol = solve_ivp(
        fun=lambda t, y: rhs_pde(t, y, pars),
        t_span=(float(times_h[0]), float(times_h[-1])),
        y0=y0,
        t_eval=times_h,
        method="BDF",
        rtol=1e-6,
        atol=1e-9,
        max_step=0.5
    )
    if not sol.success:
        raise RuntimeError("PDE solver failed: " + str(sol.message))

    Yt = sol.y.T.reshape(len(sol.t), 6, N_GRID)
    C  = Yt[:, 1, :]
    A  = Yt[:, 2, :]

    # FIX 3: WL from polymer solid only
    M0 = np.trapz(C[0, :] + A[0, :], x) / L_THICKNESS_MM
    M  = np.trapz(C + A, x, axis=1) / L_THICKNESS_MM
    WL_pct = 100.0 * (1.0 - M / max(M0, 1e-12))

    # FIX 4: Xc from polymer solid only
    denom = np.trapz(C + A, x, axis=1) / L_THICKNESS_MM
    Cbar  = np.trapz(C, x, axis=1) / L_THICKNESS_MM
    Xc_pct = 100.0 * Cbar / np.maximum(denom, 1e-12)

    return WL_pct, Xc_pct

# =========================================================
# FITTING
# =========================================================
param_names = ["k1", "km1", "k3", "km3", "kconv", "kdegC", "kdegA", "De0"]

lb = np.array([1e-6]*7 + [1e-12], dtype=float)
ub = np.array([50.0]*7 + [50.0],  dtype=float)

theta_lb = np.log(lb)
theta_ub = np.log(ub)

p0     = np.array([0.01, 0.1, 0.1, 1.0, 0.03, 0.001, 0.10, 1e-3], dtype=float)
theta0 = np.log(np.clip(p0, lb, ub))

def unpack_pars(theta):
    return np.exp(theta)

def residuals_theta(theta):
    pars = unpack_pars(theta)
    try:
        WL_all, Xc_all = simulate_observables(t_all_h, pars)
    except Exception:
        return np.full(t_P_h.size + t_C_h.size, 1e6, dtype=float)

    WL_pred = np.interp(t_P_h, t_all_h, WL_all)
    Xc_pred = np.interp(t_C_h, t_all_h, Xc_all)

    rW = (WL_pred - WL_obs_pct) / sigma_WL if USE_WEIGHTS else WL_pred - WL_obs_pct
    rX = (Xc_pred - Xc_obs_pct) / sigma_Xc if USE_WEIGHTS else Xc_pred - Xc_obs_pct

    return np.concatenate([rW, rX])

# =========================================================
# RUN
# =========================================================
print("\n=== DIMENSIONAL fit (x in mm, t in h, De0 in mm^2/h) ===")
print(f"L_THICKNESS_MM = {L_THICKNESS_MM} mm  (half of 0.5 mm film)")
print(f"N_GRID = {N_GRID} | E_BULK = {E_BULK} (fixed)")
print(f"ALPHA_E = {ALPHA_E} (Wang et al. 2008 homogenisation, fixed)")
print(f"Porosity: phi = 1 - (C+A)/(C0+A0)  [polymer solid only]")
print(f"WL: 1 - mean(C+A)/mean(C0+A0)       [polymer solid only]")
print(f"Xc: mean(C)/mean(C+A)                [polymer solid only]")
print(f"Weighted fit: {USE_WEIGHTS} | error floor: {MIN_ERR_PCT}%")
print("Fitted parameters:", ", ".join(param_names))

res = least_squares(
    fun=residuals_theta,
    x0=theta0,
    bounds=(theta_lb, theta_ub),
    method="trf",
    max_nfev=MAX_NFEV,
    ftol=1e-10,
    xtol=1e-10,
    gtol=1e-10
)

theta_hat = res.x
pars_hat  = unpack_pars(theta_hat)

print("\nLeast-squares status:", res.status, res.message)
print("\nEstimated parameters:")
for n, v in zip(param_names, pars_hat):
    unit = "mm^2/h" if n == "De0" else "1/h"
    print(f"  {n:6s} = {v:.6g}  [{unit}]")

# =========================================================
# DIAGNOSTICS
# =========================================================
def r2_score(y_obs, y_hat):
    ss_res = np.sum((y_obs - y_hat) ** 2)
    ss_tot = np.sum((y_obs - np.mean(y_obs)) ** 2)
    return (1.0 - ss_res / ss_tot) if ss_tot > 0 else np.nan

WL_all_hat, Xc_all_hat = simulate_observables(t_all_h, pars_hat)
WL_fit = np.interp(t_P_h, t_all_h, WL_all_hat)
Xc_fit = np.interp(t_C_h, t_all_h, Xc_all_hat)

r2_WL = r2_score(WL_obs_pct, WL_fit)
r2_Xc = r2_score(Xc_obs_pct, Xc_fit)
print(f"\nR^2 (Weight loss) = {r2_WL:.6f}")
print(f"R^2 (Xc)          = {r2_Xc:.6f}")

if USE_WEIGHTS:
    r = residuals_theta(theta_hat)
    chi2 = np.sum(r ** 2)
    dof = r.size - theta_hat.size
    if dof > 0:
        print(f"chi^2 = {chi2:.6g},  dof = {dof},  reduced chi^2 = {chi2/dof:.6g}")

# Dense grid for smooth plots
t_plot_h = np.linspace(0.0, float(t_all_h[-1]), 600)
WL_plot, Xc_plot = simulate_observables(t_plot_h, pars_hat)

# Residuals at experimental degradation times
WL_res = WL_fit - WL_obs_pct
Xc_res = Xc_fit - Xc_obs_pct

# =========================================================
# EXPORT RESIDUALS TO EXCEL
# =========================================================
df_wl = pd.DataFrame({
    "Time_h": t_P_h,
    "WL_Observed_pct": WL_obs_pct,
    "WL_Predicted_pct": WL_fit,
    "WL_Residual_pct": WL_res,
    "WL_ErrorBar_pct": err_WL_pct,
    "WL_Sigma_used_pct": sigma_WL
})

df_xc = pd.DataFrame({
    "Time_h": t_C_h,
    "Xc_Observed_pct": Xc_obs_pct,
    "Xc_Predicted_pct": Xc_fit,
    "Xc_Residual_pct": Xc_res,
    "Xc_ErrorBar_pct": err_Xc_pct,
    "Xc_Sigma_used_pct": sigma_Xc
})

t_union = np.unique(np.concatenate([t_P_h, t_C_h]))
df_combined = pd.DataFrame({"Time_h": t_union})

df_combined = df_combined.merge(
    df_wl[["Time_h", "WL_Observed_pct", "WL_Predicted_pct", "WL_Residual_pct"]],
    on="Time_h", how="left"
)

df_combined = df_combined.merge(
    df_xc[["Time_h", "Xc_Observed_pct", "Xc_Predicted_pct", "Xc_Residual_pct"]],
    on="Time_h", how="left"
)

excel_out = os.path.join(download_dir, "residuals_at_degradation_times.xlsx")
with pd.ExcelWriter(excel_out, engine="openpyxl") as writer:
    df_wl.to_excel(writer, sheet_name="WeightLoss_Residuals", index=False)
    df_xc.to_excel(writer, sheet_name="Crystallinity_Residuals", index=False)
    df_combined.to_excel(writer, sheet_name="Combined_Residuals", index=False)

print(f"Residual tables saved to Excel: {excel_out}")

# =========================================================
# PLOTS
# =========================================================

# ---- Plot 1: Weight loss ----
plt.figure(figsize=(7, 5))
plt.errorbar(
    t_P_h, WL_obs_pct, yerr=err_WL_pct, fmt="o",
    color="black", capsize=4, label="PCL-Lipase data"
)
plt.plot(t_plot_h, WL_plot, "r", linewidth=3, label="Model fit")
plt.xlabel("Time (hours)", fontsize=14)
plt.ylabel("Weight loss (%)", fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=11)
plt.grid(True, alpha=0.3)
plt.tight_layout()
out1 = os.path.join(download_dir, "fit_weight_loss.png")
plt.savefig(out1, dpi=600)
plt.show()
print(f"Saved: {out1}")

# ---- Plot 2: Crystallinity ----
plt.figure(figsize=(7, 5))
plt.errorbar(
    t_C_h, Xc_obs_pct, yerr=err_Xc_pct, fmt="s",
    color="black", capsize=4, label=r"$\chi_c$ data"
)
plt.plot(t_plot_h, Xc_plot, "r", linewidth=3, label="Model fit")
plt.xlabel("Time (hours)", fontsize=14)
plt.ylabel(r"$\chi_c$ (%)", fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=11)
plt.grid(True, alpha=0.3)
plt.tight_layout()
out2 = os.path.join(download_dir, "fit_crystallinity.png")
plt.savefig(out2, dpi=600)
plt.show()
print(f"Saved: {out2}")

# ---- Plot 3: WL residuals ----
plt.figure(figsize=(7, 4))
plt.axhline(0.0, linewidth=1, color="black")
plt.plot(t_P_h, WL_res, "o", color="green", label="Residuals")
plt.ylim([-10, 10])
plt.xlabel("Time (hours)", fontsize=14)
plt.ylabel("Weight loss residual (%)", fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid(True, alpha=0.3)
plt.legend(fontsize=10, loc="upper left", frameon=False)
plt.tight_layout()
out3 = os.path.join(download_dir, "residuals_weight_loss.png")
plt.savefig(out3, dpi=600)
plt.show()
print(f"Saved: {out3}")

# ---- Plot 4: Xc residuals ----
plt.figure(figsize=(7, 4))
plt.axhline(0.0, linewidth=1, color="black")
plt.plot(t_C_h, Xc_res, "o", color="blue", label="Residuals")
plt.ylim([-10, 10])
plt.xlabel("Time (hours)", fontsize=14)
plt.ylabel(r"$\chi_c$ residual (%)", fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid(True, alpha=0.3)
plt.legend(fontsize=10, loc="upper left", frameon=False)
plt.tight_layout()
out4 = os.path.join(download_dir, "residuals_crystallinity.png")
plt.savefig(out4, dpi=600)
plt.show()
print(f"Saved: {out4}")
