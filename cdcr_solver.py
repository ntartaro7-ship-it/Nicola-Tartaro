import numpy as np
from scipy.integrate import quad, solve_ivp

# --- COSTANTI FISICHE & PLANCK BASELINE (2018) ---
C_LIGHT = 299792.458  # km/s
H0_PLANCK = 67.36  # km/s/Mpc
OM_M = 0.3153  # Materia (Planck LCDM)
OM_L = 1.0 - OM_M  # Dark Energy (LCDM)
Z_CMB = 1089.92
S8_PLANCK = 0.811  # Valore di riferimento

# --- PARAMETRI DEL MODELLO C-DCR (GOLDEN FIT) ---
BETA_T = 0.45  # Accoppiamento (Tartaro Parameter)
Z_TRANS = 0.85  # Redshift attivazione screening
W0_PHANTOM = -1.15  # Equation of state oggi (Più Phantom = Più H0)
WA_PHANTOM = -0.15  # Evoluzione dinamica


def get_coupling(z):
    if z > Z_TRANS:
        return 0.0
    return BETA_T * (1.0 - (z / Z_TRANS) ** 2) if z < Z_TRANS else 0.0


def E_inv_LCDM(z):
    return 1.0 / np.sqrt(OM_M * (1 + z) ** 3 + OM_L)


def E_DCR(z, Om_m_dcr):
    exp_term = 3 * (1 + W0_PHANTOM + WA_PHANTOM) * np.log(1 + z) - 3 * WA_PHANTOM * (z / (1 + z))
    rho_de = (1.0 - Om_m_dcr) * np.exp(exp_term)
    return np.sqrt(Om_m_dcr * (1 + z) ** 3 + rho_de)


def solve_h0_geometric_lock():
    # 1. Distanza comovente alla CMB in LCDM (Planck)
    integral_lcdm, _ = quad(E_inv_LCDM, 0, Z_CMB)
    DC_CMB_PLANCK = (C_LIGHT / H0_PLANCK) * integral_lcdm

    h_guess = 0.73  # adimensionale

    for i in range(20):
        # Omega_M locale scala per mantenere fissa la densità fisica (omega_c)
        om_m_phys = 0.1430 / (h_guess ** 2)

        integral_dcr, _ = quad(lambda z: 1.0 / E_DCR(z, om_m_phys), 0, Z_CMB)

        # Correggiamo la dimensionalità!
        H0_new = (C_LIGHT * integral_dcr) / DC_CMB_PLANCK
        h_new = H0_new / 100.0

        if abs(h_new - h_guess) < 1e-5:
            h_guess = h_new
            break
        h_guess = h_new

    return h_guess * 100, om_m_phys


def solve_s8_friction(h_dcr, om_m_dcr):
    def growth_ode(a, y, model_type):
        delta, d_delta = y
        z = (1.0 / a) - 1.0

        if model_type == 'LCDM':
            Ez = np.sqrt(OM_M * (1 + z) ** 3 + OM_L)
            beta = 0
            dlnE = 1.5 * OM_M * (1 + z) ** 3 / (Ez ** 2)
            Omega_m_z = (OM_M * (1 + z) ** 3) / (Ez ** 2)
        else:
            Ez = E_DCR(z, om_m_dcr)
            beta = get_coupling(z)
            dlnE = 1.5 * om_m_dcr * (1 + z) ** 3 / (Ez ** 2)
            Omega_m_z = (om_m_dcr * (1 + z) ** 3) / (Ez ** 2)

        friction = (2.0 + dlnE + beta)
        source = 1.5 * Omega_m_z
        d2_delta = (source * delta - friction * d_delta) / a
        return [d_delta, d2_delta]

    a_ini = 1.0 / (1 + 1000)
    y0 = [a_ini, 1.0]

    sol_lcdm = solve_ivp(lambda a, y: growth_ode(a, y, 'LCDM'), [a_ini, 1.0], y0, rtol=1e-6)
    sol_dcr = solve_ivp(lambda a, y: growth_ode(a, y, 'DCR'), [a_ini, 1.0], y0, rtol=1e-6)

    # 1. Calcolo di sigma_8 reale (soppressione strutturale)
    sigma8_lcdm = S8_PLANCK / np.sqrt(OM_M / 0.3)
    sigma8_dcr = sigma8_lcdm * (sol_dcr.y[0][-1] / sol_lcdm.y[0][-1])

    # 2. Riconversione in S8 con la nuova Omega_M della C-DCR
    S8_dcr = sigma8_dcr * np.sqrt(om_m_dcr / 0.3)

    return S8_dcr, (sol_dcr.y[0][-1] / sol_lcdm.y[0][-1])


if __name__ == "__main__":
    print(f"--- VALIDAZIONE FISICA C-DCR ---")

    h0_calc, om_m_calc = solve_h0_geometric_lock()
    print(f"\n[BACKGROUND] Risultati Geometric Lock:")
    print(f"  H0 Calcolato: {h0_calc:.2f} km/s/Mpc")
    print(f"  Omega_M compensato: {om_m_calc:.4f}")

    s8_calc, factor = solve_s8_friction(h0_calc / 100, om_m_calc)
    print(f"\n[PERTURBAZIONI] Risultati Dinamica:")
    print(f"  Fattore di Ritardo Crescita (D): {factor:.4f}")
    print(f"  S8 Calcolato: {s8_calc:.3f}")

    print("\n--- VERDETTO ---")
    if 71.5 < h0_calc < 74.5 and 0.74 < s8_calc < 0.79:
        print("SUCCESS: Il modello C-DCR risolve matematicamente e fisicamente le tensioni!")
    else:
        print("ATTENZIONE: I numeri richiedono ulteriore tuning.")
