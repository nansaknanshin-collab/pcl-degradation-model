import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import pandas as pd
import os

# =========================================================
# THICKNESS VARIATION STUDY
# Fixed parameters from least-squares calibration to Shi et al.
# =========================================================

PARS_HAT = np.array([
    0.651327,   # k1    [1/h]
    0.474921,   # km1   [1/h]
    0.300349,   # k3    [1/h]
    0.0282859,  # km3   [1/h]
    0.575375,   # kconv [1/h]
    0.00116811, # kdegC [1/h]
    0.00492959, # kdegA [1/h]
    1.65766,    # De0   [mm^2/h]
])

N_GRID   = 25
C0_POLY  = 0.302
A0_POLY  = 0.698
E_INIT   = 0.0
P_INIT   = 0.0
E_BULK   = 0.5
ALPHA_E  = 1.0

# These are MODEL HALF-THICKNESSES L [mm], because the PDE is solved on x in [0, L]
MODEL_HALF_THICKNESSES = [0.1, 0.25, 0.5, 3.0, 4.0]

# Corresponding FULL PHYSICAL FILM THICKNESSES [mm]
FULL_THICKNESSES = [2 * L for L in MODEL_HALF_THICKNESSES]

LABELS = [
    rf"Film thickness = {full:.2f} mm" + (" (case study)" if i == 1 else "")
    for i, full in enumerate(FULL_THICKNESSES)
]

COLORS = ["tab:green", "tab:red", "tab:blue", "tab:purple", "tab:orange"]

# Index of the calibration case in the lists above
# L = 0.25 mm half-thickness <=> 0.50 mm full thickness
CAL_IDX = 1

t_plot_h = np.linspace(0.0, 72.0, 600)
download_dir = os.getcwd()

# =========================================================
# DIFFUSION OPERATOR
# =========================================================
def diffusion_op(u, D, u_left, dx, n):
    u = u.copy()
    u[0] = u_left

    flux = np.zeros(n + 1)
    D_iface = 0.5 * (D[:-1] + D[1:])

    flux[1:n] = D_iface * (u[1:] - u[:-1]) / dx
    flux[n] = 0.0

    div = np.zeros(n)
    div[1:-1] = (flux[2:n] - flux[1:n-1]) / dx
    div[-1] = 2.0 * (flux[n] - flux[n-1]) / dx
    div[0] = 0.0

    return div

# =========================================================
# PDE RHS
# =========================================================
def make_rhs(L, n):
    x = np.linspace(0.0, L, n)
    dx = x[1] - x[0]

    def rhs(t, y, pars):
        k1, km1, k3, km3, kconv, kdegC, kdegA, De0 = pars

        Y = y.reshape(6, n)
        E, C, A, EC, EA, P = Y

        E = E.copy()
        E[0] = E_BULK

        denom0 = C0_POLY + A0_POLY
        phi = 1.0 - (C + A) / max(denom0, 1e-12)
        phi = np.clip(phi, -0.5, 2.0)

        De = np.maximum(De0 * (1.0 + ALPHA_E * phi), 1e-30)
        diffE = diffusion_op(E, De, E_BULK, dx, n)

        dE = (
            diffE
            - k1 * E * C + km1 * EC
            - k3 * E * A + km3 * EA
            + kconv * EC
            + kdegC * EC + kdegA * EA
        )
        dC = -k1 * E * C + km1 * EC
        dA = -k3 * E * A + km3 * EA + kconv * EC
        dEC = k1 * E * C - km1 * EC - kconv * EC - kdegC * EC
        dEA = k3 * E * A - km3 * EA - kdegA * EA
        dP = kdegC * EC + kdegA * EA

        dE[0] = 0.0

        return np.vstack([dE, dC, dA, dEC, dEA, dP]).reshape(-1)

    return rhs, x, dx

# =========================================================
# SIMULATE
# =========================================================
def simulate(L, pars, times):
    n = N_GRID
    rhs, x, dx = make_rhs(L, n)

    E0 = np.full(n, E_INIT)
    C0v = np.full(n, C0_POLY)
    A0v = np.full(n, A0_POLY)

    y0 = np.vstack([
        E0,
        C0v,
        A0v,
        np.zeros(n),
        np.zeros(n),
        np.full(n, P_INIT)
    ]).reshape(-1)

    y0[0] = E_BULK

    sol = solve_ivp(
        fun=lambda t, y: rhs(t, y, pars),
        t_span=(float(times[0]), float(times[-1])),
        y0=y0,
        t_eval=times,
        method="BDF",
        rtol=1e-6,
        atol=1e-9,
        max_step=0.5
    )

    if not sol.success:
        raise RuntimeError(f"Solver failed (L={L}): {sol.message}")

    Yt = sol.y.T.reshape(len(sol.t), 6, n)
    C = Yt[:, 1, :]
    A = Yt[:, 2, :]

    M0 = np.trapz(C[0] + A[0], x) / L
    M = np.trapz(C + A, x, axis=1) / L
    WL_pct = 100.0 * (1.0 - M / max(M0, 1e-12))

    denom = np.trapz(C + A, x, axis=1) / L
    Cbar = np.trapz(C, x, axis=1) / L
    Xc_pct = 100.0 * Cbar / np.maximum(denom, 1e-12)

    return WL_pct, Xc_pct

# =========================================================
# RUN SWEEP
# =========================================================
print("Running thickness variation study...")
results = []

for L, full in zip(MODEL_HALF_THICKNESSES, FULL_THICKNESSES):
    print(f"  Simulating half-thickness L = {L:.3f} mm (full thickness = {full:.2f} mm)...", end=" ")
    WL, Xc = simulate(L, PARS_HAT, t_plot_h)
    results.append((WL, Xc))
    tau = L**2 / PARS_HAT[-1]
    print(f"done  (tau_diff = {tau:.4f} h)")

# =========================================================
# PLOT 1: ALL THICKNESSES — WEIGHT LOSS
# =========================================================
fig, ax = plt.subplots(figsize=(9, 9))

for i, ((WL, _), label, color) in enumerate(zip(results, LABELS, COLORS)):
    ls = "--" if i == CAL_IDX else "-"
    ax.plot(t_plot_h, WL, color=color, linewidth=3, linestyle=ls, label=label)

ax.set_xlabel("Time (hours)", fontsize=16)
ax.set_ylabel("Weight loss (%)", fontsize=16)
ax.legend(fontsize=10, framealpha=0.9)
ax.set_xlim([0, 72])
ax.set_ylim([0, 100])
ax.grid(True, alpha=0.3)

plt.tight_layout()
out1 = os.path.join(download_dir, "thickness_WL.png")
plt.savefig(out1, dpi=600)
plt.show()
print(f"Saved: {out1}")

# =========================================================
# PLOT 2: ALL THICKNESSES — CRYSTALLINITY
# =========================================================
fig, ax = plt.subplots(figsize=(9, 9))

for i, ((_, Xc), label, color) in enumerate(zip(results, LABELS, COLORS)):
    ls = "--" if i == CAL_IDX else "-"
    ax.plot(t_plot_h, Xc, color=color, linewidth=3, linestyle=ls, label=label)

ax.set_xlabel("Time (hours)", fontsize=16)
ax.set_ylabel(r"$\chi_c$ (%)", fontsize=16)
ax.legend(fontsize=10, framealpha=0.9)
ax.set_xlim([0, 72])
ax.set_ylim([0, 35])
ax.grid(True, alpha=0.3)

plt.tight_layout()
out2 = os.path.join(download_dir, "thickness_Xc.png")
plt.savefig(out2, dpi=600)
plt.show()
print(f"Saved: {out2}")

# =========================================================
# PLOT 3: EACH THICKNESS vs CALIBRATION — WEIGHT LOSS
# =========================================================
non_cal_indices = [i for i in range(len(MODEL_HALF_THICKNESSES)) if i != CAL_IDX]

fig, axes = plt.subplots(2, 2, figsize=(9, 9), sharey=True)
axes = axes.flatten()

WL_cal, Xc_cal = results[CAL_IDX]

for ax, idx in zip(axes, non_cal_indices):
    WL_other, _ = results[idx]

    ax.plot(
        t_plot_h, WL_cal,
        color="tab:red", linewidth=3.0, linestyle="--",
        label=LABELS[CAL_IDX]
    )
    ax.plot(
        t_plot_h, WL_other,
        color=COLORS[idx], linewidth=3.0, linestyle="-",
        label=LABELS[idx]
    )

    ax.set_xlabel("Time (hours)", fontsize=16)
    ax.set_ylabel("Weight loss (%)", fontsize=16)
    ax.set_xlim([0, 72])
    ax.set_ylim([0, 100])
    ax.legend(fontsize=9, framealpha=0.9)
    ax.grid(True, alpha=0.3)

plt.tight_layout()
out3 = os.path.join(download_dir, "comparison_WL.png")
plt.savefig(out3, dpi=600, bbox_inches="tight")
plt.show()
print(f"Saved: {out3}")

# =========================================================
# PLOT 4: EACH THICKNESS vs CALIBRATION — CRYSTALLINITY
# =========================================================
fig, axes = plt.subplots(2, 2, figsize=(9, 9), sharey=True)
axes = axes.flatten()

for ax, idx in zip(axes, non_cal_indices):
    _, Xc_other = results[idx]

    ax.plot(
        t_plot_h, Xc_cal,
        color="tab:red", linewidth=3.0, linestyle="--",
        label=LABELS[CAL_IDX]
    )
    ax.plot(
        t_plot_h, Xc_other,
        color=COLORS[idx], linewidth=3.0, linestyle="-",
        label=LABELS[idx]
    )

    ax.set_xlabel("Time (hours)", fontsize=16)
    ax.set_ylabel(r"$\chi_c$ (%)", fontsize=16)
    ax.set_xlim([0, 72])
    ax.set_ylim([0, 35])
    ax.legend(fontsize=9, framealpha=0.9)
    ax.grid(True, alpha=0.3)

plt.tight_layout()
out4 = os.path.join(download_dir, "comparison_Xc.png")
plt.savefig(out4, dpi=600, bbox_inches="tight")
plt.show()
print(f"Saved: {out4}")

# =========================================================
# SUMMARY TABLE
# =========================================================
print("\n=== Summary: Final values at t = 72 h ===")
print(f"{'Full thickness':>18} | {'WL(72h) %':>10} | {'Xc(72h) %':>10} | {'tau_diff (h)':>12}")
print("-" * 65)

for L, full, (WL, Xc) in zip(MODEL_HALF_THICKNESSES, FULL_THICKNESSES, results):
    tau = L**2 / PARS_HAT[-1]
    print(f"{full:>16.2f} mm | {WL[-1]:>10.2f} | {Xc[-1]:>10.2f} | {tau:>12.4f}")

# =========================================================
# EXPORT TO EXCEL
# =========================================================
excel_path = os.path.join(download_dir, "thickness_predictions.xlsx")

with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
    for L, full, (WL, Xc) in zip(MODEL_HALF_THICKNESSES, FULL_THICKNESSES, results):
        sheet_name = f"full_{full:.2f}mm".replace(".", "_")
        df = pd.DataFrame({
            "Time (h)": t_plot_h,
            "Weight Loss (%)": WL,
            "Xc (%)": Xc,
            "Model half-thickness L (mm)": [L] * len(t_plot_h),
            "Film thickness (mm)": [full] * len(t_plot_h),
        })
        df.to_excel(writer, sheet_name=sheet_name, index=False)

    summary = pd.DataFrame({"Time (h)": t_plot_h})
    for full, (WL, Xc) in zip(FULL_THICKNESSES, results):
        summary[f"WL - {full:.2f} mm (%)"] = WL
    for full, (WL, Xc) in zip(FULL_THICKNESSES, results):
        summary[f"Xc - {full:.2f} mm (%)"] = Xc

    summary.to_excel(writer, sheet_name="Summary", index=False)

print(f"Saved Excel: {excel_path}")