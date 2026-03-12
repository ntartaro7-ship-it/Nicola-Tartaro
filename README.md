import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.integrate import quad, solve_ivp
import os

# ==========================================
# COSTANTI E SETUP
# ==========================================
C_LIGHT = 299792.458
H0_GUESS = 73.13
FILENAME = r"C:\Users\nicol\PyCharmMiscProject\NGC5055_rotmod.dat"

OM_M = 0.3153
OM_L = 1.0 - OM_M
S8_PLANCK = 0.811
Z_TRANS = 0.85
W0_PHANTOM = -1.15
WA_PHANTOM = -0.15


# ==========================================
# 1. SCALA GALATTICA: ESTRAZIONE SIGMA_BETA
# ==========================================
def load_data():
    if not os.path.exists(FILENAME): return None
    cols = ['Rad', 'Vobs', 'errV', 'Vgas', 'Vdisk', 'Vbul', 'SBdisk', 'SBbul']
    df = pd.read_csv(FILENAME, sep=r'\s+', engine='python', comment='#', names=cols)
    return df[df['Rad'] > 0.1]


df = load_data()


def modello_unificato_v_total(rad, beta):
    # Legge di scala unificata D-DCR
    a_crit = (beta ** 2 / 2.0) * C_LIGHT * H0_GUESS

    # Massa barionica (beta funge da M/L ratio del disco)
    v_bar_sq = (df['Vgas'] ** 2) + (beta * df['Vdisk'] ** 2) + (0.6 * df['Vbul'] ** 2)
    g_bar = np.maximum(v_bar_sq / rad, 1e-12)

    # Profilo RAR
    x = np.sqrt(g_bar / a_crit)
    g_obs = g_bar / (1 - np.exp(-x))
    return np.sqrt(g_obs * rad)


print("--- 1. SCALA GALATTICA ---")
if df is not None:
    # FIT CON LIMITI FISICI: beta deve essere tra 0.1 e 1.0 (masse positive e sensate)
    popt, pcov = curve_fit(modello_unificato_v_total, df['Rad'], df['Vobs'],
                           p0=[0.4404], bounds=(0.1, 1.0),
                           sigma=df['errV'], absolute_sigma=True)
    beta_real = popt[0]
    sigma_beta = np.sqrt(pcov[0][0])
    print(f"Beta ottimizzato: {beta_real:.4f}")
    print(f"Sigma_beta REALE: {sigma_beta:.4f}\n")
else:
    beta_real = 0.4404
    sigma_beta = 0.0
    print("ERRORE: File NGC5055_rotmod.dat non trovato.\n")


# ==========================================
# 2. SCALA COSMOLOGICA: ESTRAZIONE SIGMA_D
# ==========================================
def get_coupling(z, beta_val):
    val = beta_val * (1.0 - (z / Z_TRANS) ** 2)
    smooth_step = 0.5 * (1.0 - np.tanh(25.0 * (z - Z_TRANS)))
    return max(0.0, val) * smooth_step


def E_DCR(z, Om_m_dcr):
    exp_term = 3 * (1 + W0_PHANTOM + WA_PHANTOM) * np.log(1 + z) - 3 * WA_PHANTOM * (z / (1 + z))
    rho_de = (1.0 - Om_m_dcr) * np.exp(exp_term)
    return np.sqrt(Om_m_dcr * (1 + z) ** 3 + rho_de)


def dlnE_dlnA_DCR(z, Om_m_dcr):
    dz = 1e-5
    E_plus = E_DCR(z + dz, Om_m_dcr)
    E_minus = E_DCR(z - dz, Om_m_dcr)
    E_current = E_DCR(z, Om_m_dcr)
    dE_dz = (E_plus - E_minus) / (2 * dz)
    return -(1 + z) * dE_dz / E_current


# ==========================================
# 2. SCALA COSMOLOGICA CORRETTA: ESTRAZIONE SIGMA_D
# ==========================================
def calcola_s8_dinamico(D_drag_val, beta_val):
    om_m_dcr = 0.1430 / ((H0_GUESS / 100) ** 2)

    # ODE per D-DCR
    def growth_ode_DCR(N, y):
        delta, d_delta_dN = y
        a = np.exp(N)
        z = (1.0 / a) - 1.0
        Ez = E_DCR(z, om_m_dcr)
        Omega_m_z = (om_m_dcr * (1 + z) ** 3) / (Ez ** 2)
        beta_eff = get_coupling(z, beta_val)

        phi_prime = beta_eff * Omega_m_z
        # Frizione Lineare Effettiva
        disformal_friction = D_drag_val * (beta_val ** 2) * abs(phi_prime)

        dlnE = dlnE_dlnA_DCR(z, om_m_dcr)
        friction = 2.0 + dlnE + disformal_friction
        G_eff_factor = 1.0 + 2.0 * (beta_eff ** 2)
        source = 1.5 * Omega_m_z * G_eff_factor
        return [d_delta_dN, source * delta - friction * d_delta_dN]

    # ODE per LCDM (Baseline)
    def growth_ode_LCDM(N, y):
        delta, d_delta_dN = y
        a = np.exp(N)
        z = (1.0 / a) - 1.0
        Ez = np.sqrt(OM_M * (1 + z) ** 3 + OM_L)
        Omega_m_z = (OM_M * (1 + z) ** 3) / Ez ** 2
        dlnE = -(1 + z) * (1.5 * OM_M * (1 + z) ** 2) / (Ez ** 2)
        friction = 2.0 + dlnE
        source = 1.5 * Omega_m_z
        return [d_delta_dN, source * delta - friction * d_delta_dN]

    a_ini = 1.0 / 1001.0
    sol_dcr = solve_ivp(growth_ode_DCR, [np.log(a_ini), 0], [a_ini, a_ini], method='LSODA', rtol=1e-8, atol=1e-10)
    sol_lcdm = solve_ivp(growth_ode_LCDM, [np.log(a_ini), 0], [a_ini, a_ini], method='LSODA', rtol=1e-8, atol=1e-10)

    # Rapporto di crescita esatto come nel tuo script originale
    ratio = sol_dcr.y[0][-1] / sol_lcdm.y[0][-1]
    sigma8_dcr = (S8_PLANCK / np.sqrt(OM_M / 0.3)) * ratio
    return sigma8_dcr * np.sqrt(om_m_dcr / 0.3)


print("--- 2. SCALA COSMOLOGICA ---")
D_nom = 13.0543
beta_fisso = 0.4404  # Usiamo il prior teorico globale
eps = 0.5

s8_nominale = calcola_s8_dinamico(D_nom, beta_fisso)
s8_plus = calcola_s8_dinamico(D_nom + eps, beta_fisso)
s8_minus = calcola_s8_dinamico(D_nom - eps, beta_fisso)

derivata_S8_su_D = (s8_plus - s8_minus) / (2 * eps)
sigma_S8_obs = 0.017  # Errore survey DES/KiDS

print(f"S8 Calcolato con D={D_nom}: {s8_nominale:.4f}")
print(f"Derivata d(S8)/d(D): {derivata_S8_su_D:.6f}")

if abs(derivata_S8_su_D) > 1e-6:
    sigma_D = abs(1.0 / derivata_S8_su_D) * sigma_S8_obs
    print(f"Sigma_D REALE: {sigma_D:.4f}\n")
else:
    print("Derivata ancora troppo debole.")
