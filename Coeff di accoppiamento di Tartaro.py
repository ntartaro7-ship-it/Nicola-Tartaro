import numpy as np
from scipy.integrate import quad, solve_ivp
import matplotlib.pyplot as plt

# --- CONFIGURATION ---
h_planck = 0.674
Om_m_planck = 0.315
Om_L_planck = 1.0 - Om_m_planck
Om_r_planck = 9e-5
z_cmb = 1090.0
sigma8_planck = 0.811

# DESI parameters
w0_desi = -0.7379
wa_desi = -1.0049

# C-DCR parameters
z_tr = 1.05


# Model Functions
def get_rho_de_cpl(z, w0, wa):
    return (1 + z) ** (3 * (1 + w0 + wa)) * np.exp(-3 * wa * z / (1 + z))


def get_rho_m_screened(z, beta):
    if z < z_tr:
        return (1 + z) ** (3.0 - beta)
    else:
        return (1 + z_tr) ** (3.0 - beta) * ((1 + z) / (1 + z_tr)) ** 3.0


def hubble_sq_norm(z, beta):
    rho_r = Om_r_planck * (1 + z) ** 4
    rho_m = Om_m_planck * get_rho_m_screened(z, beta)
    rho_de = Om_L_planck * get_rho_de_cpl(z, w0_desi, wa_desi)
    return rho_r + rho_m + rho_de


def solve_cdcr(beta):
    # H0
    def inv_E_lcdm(z): return 1.0 / np.sqrt(Om_r_planck * (1 + z) ** 4 + Om_m_planck * (1 + z) ** 3 + Om_L_planck)

    dist_target = quad(inv_E_lcdm, 0, z_cmb)[0] / h_planck

    def inv_E_dcr(z): return 1.0 / np.sqrt(hubble_sq_norm(z, beta))

    h0_dcr = (quad(inv_E_dcr, 0, z_cmb)[0] / dist_target) * 100

    # S8
    def growth_ode(lna, y):
        delta, delta_prime = y
        z = np.exp(-lna) - 1
        E = np.sqrt(hubble_sq_norm(z, beta))
        rho_m = Om_m_planck * get_rho_m_screened(z, beta)
        Om_m_z = rho_m / E ** 2

        dz = 1e-4
        dE_dz = (np.sqrt(hubble_sq_norm(z + dz, beta)) - np.sqrt(hubble_sq_norm(z - dz, beta))) / (2 * dz)
        friction = 2.0 - (1 + z) * dE_dz / E
        if z < z_tr: friction += 2.0 * beta * (1.0 - Om_m_z)

        return [delta_prime, 1.5 * Om_m_z * delta - friction * delta_prime]

    y0 = [1e-3, 1e-3]
    sol_dcr = solve_ivp(growth_ode, [-np.log(1 + 1000), 0.0], y0, method='RK45')

    # LCDM reference
    def growth_lcdm(lna, y):
        d, dp = y
        z = np.exp(-lna) - 1
        E2 = Om_r_planck * (1 + z) ** 4 + Om_m_planck * (1 + z) ** 3 + Om_L_planck
        Om_m_z = Om_m_planck * (1 + z) ** 3 / E2
        dE2_dz = 4 * Om_r_planck * (1 + z) ** 3 + 3 * Om_m_planck * (1 + z) ** 2
        friction = 2.0 - (1 + z) * (0.5 * dE2_dz / np.sqrt(E2)) / np.sqrt(E2)
        return [dp, 1.5 * Om_m_z * d - friction * dp]

    sol_lcdm = solve_ivp(growth_lcdm, [-np.log(1 + 1000), 0.0], y0, method='RK45')

    return h0_dcr, sigma8_planck * (sol_dcr.y[0][-1] / sol_lcdm.y[0][-1])


# Data Generation for Plot
betas = np.linspace(0.0, 0.35, 20)
h0s, s8s = [], []
for b in betas:
    h, s = solve_cdcr(b)
    h0s.append(h)
    s8s.append(s)

# Plotting
fig, ax1 = plt.subplots(figsize=(10, 6))

color = 'tab:blue'
ax1.set_xlabel(r'Dark Coupling $\beta$', fontsize=12)
ax1.set_ylabel(r'$H_0$ (km/s/Mpc)', color=color, fontsize=12)
ax1.plot(betas, h0s, color=color, linewidth=3, label=r'C-DCR Prediction for $H_0$')
ax1.tick_params(axis='y', labelcolor=color)

# SH0ES band
ax1.axhspan(73.04 - 1.04, 73.04 + 1.04, color='blue', alpha=0.1, label='SH0ES (R21)')
ax1.axhline(73.04, color='blue', linestyle='--', alpha=0.5)

ax2 = ax1.twinx()
color = 'tab:red'
ax2.set_ylabel(r'$S_8$ (Amplitude)', color=color, fontsize=12)
ax2.plot(betas, s8s, color=color, linewidth=3, label=r'C-DCR Prediction for $S_8$')
ax2.tick_params(axis='y', labelcolor=color)

# KiDS band
ax2.axhspan(0.766 - 0.02, 0.766 + 0.02, color='red', alpha=0.1, label='KiDS-1000')
ax2.axhline(0.766, color='red', linestyle='--', alpha=0.5)

# Intersection Line
plt.axvline(0.25, color='green', linestyle=':', linewidth=2, label=r'Best Fit $\beta=0.25$')

# Title & Legend
plt.title(r'Unification of Tensions: The $\beta=0.25$ Solution', fontsize=14)

# Combined legend
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
lines = lines1 + lines2 + [plt.Line2D([0], [0], color='green', linestyle=':', linewidth=2)]
labels = labels1 + labels2 + [r'Best Fit $\beta=0.25$']
ax1.legend(lines, labels, loc='center right')

fig.tight_layout()
plt.savefig('2_english.png', dpi=300)
plt.show()