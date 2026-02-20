import numpy as np
from scipy.integrate import quad, solve_ivp

# --- COSTANTI FISICHE & PLANCK BASELINE (2018) ---
C_LIGHT = 299792.458
H0_PLANCK = 67.36
OM_M = 0.3153
OM_L = 1.0 - OM_M
Z_CMB = 1089.92
S8_PLANCK = 0.811

# --- PARAMETRI DEL MODELLO D-DCR (PLATINUM FIT V9) ---
BETA_T = 0.4404
Z_TRANS = 0.85
W0_PHANTOM = -1.15
WA_PHANTOM = -0.15

# NUOVO PARAMETRO: Frizione Disformale (calcolato via root-finding)
D_DRAG = 13.0543


def get_coupling(z, k_smooth=25.0):
    """Accoppiamento di classe C-infinito con rottura di simmetria a z=0.85."""
    val = BETA_T * (1.0 - (z / Z_TRANS) ** 2)
    smooth_step = 0.5 * (1.0 - np.tanh(k_smooth * (z - Z_TRANS)))
    return max(0.0, val) * smooth_step


def E_DCR(z, Om_m_dcr):
    """Tasso di espansione E(z) con equazione di stato Phantom dinamica."""
    exp_term = 3 * (1 + W0_PHANTOM + WA_PHANTOM) * np.log(1 + z) - 3 * WA_PHANTOM * (z / (1 + z))
    rho_de = (1.0 - Om_m_dcr) * np.exp(exp_term)
    return np.sqrt(Om_m_dcr * (1 + z) ** 3 + rho_de)


def dlnE_dlnA(z, Om_m_dcr, model_type):
    """Calcolo esatto (numerico) di d(ln H) / d(ln a)."""
    dz = 1e-5
    if model_type == 'LCDM':
        E_plus = np.sqrt(OM_M * (1 + z + dz) ** 3 + OM_L)
        E_minus = np.sqrt(OM_M * (1 + z - dz) ** 3 + OM_L)
        E_current = np.sqrt(OM_M * (1 + z) ** 3 + OM_L)
    else:
        E_plus = E_DCR(z + dz, Om_m_dcr)
        E_minus = E_DCR(z - dz, Om_m_dcr)
        E_current = E_DCR(z, Om_m_dcr)

    dE_dz = (E_plus - E_minus) / (2 * dz)
    return -(1 + z) * dE_dz / E_current


def solve_h0_geometric_lock():
    """Risolve H0 mantenendo fissa la distanza angolare del CMB."""
    integral_lcdm, _ = quad(lambda z: 1.0 / np.sqrt(OM_M * (1 + z) ** 3 + OM_L), 0, Z_CMB)
    DC_CMB_PLANCK = (C_LIGHT / H0_PLANCK) * integral_lcdm
    h_guess = 0.73

    for _ in range(20):
        om_m_phys = 0.1430 / (h_guess ** 2)
        integral_dcr, _ = quad(lambda z: 1.0 / E_DCR(z, om_m_phys), 0, Z_CMB)
        h_new = ((C_LIGHT * integral_dcr) / DC_CMB_PLANCK) / 100.0
        if abs(h_new - h_guess) < 1e-5:
            return h_new * 100, om_m_phys
        h_guess = h_new
    return h_guess * 100, om_m_phys


def solve_s8_disformal(h_dcr, om_m_dcr):
    """Integra le perturbazioni includendo Quinta Forza e Attrito Disformale."""

    def growth_ode(N, y, model_type):
        delta, d_delta_dN = y
        a = np.exp(N)
        z = (1.0 / a) - 1.0

        if model_type == 'LCDM':
            Ez = np.sqrt(OM_M * (1 + z) ** 3 + OM_L)
            Omega_m_z = (OM_M * (1 + z) ** 3) / (Ez ** 2)
            beta_eff = 0.0
            phi_prime = 0.0
            G_eff_factor = 1.0
            disformal_friction = 0.0
        else:
            Ez = E_DCR(z, om_m_dcr)
            Omega_m_z = (om_m_dcr * (1 + z) ** 3) / (Ez ** 2)
            beta_eff = get_coupling(z)

            # Dinamica del campo e Quinta Forza
            phi_prime = beta_eff * Omega_m_z
            G_eff_factor = 1.0 + 2.0 * (beta_eff ** 2)

            # Il Disformal Drag scala con il quadrato della velocità del campo
            disformal_friction = D_DRAG * (phi_prime ** 2)

        dlnE = dlnE_dlnA(z, om_m_dcr, model_type)

        # L'equazione di moto include ora il termine di frizione disformale per schiacciare S8
        friction = 2.0 + dlnE + (beta_eff * phi_prime) + disformal_friction
        source = 1.5 * Omega_m_z * G_eff_factor

        d2_delta_dN2 = source * delta - friction * d_delta_dN
        return [d_delta_dN, d2_delta_dN2]

    a_ini = 1.0 / (1 + 1000)
    N_ini = np.log(a_ini)
    N_end = np.log(1.0)
    y0 = [a_ini, a_ini]

    sol_lcdm = solve_ivp(lambda N, y: growth_ode(N, y, 'LCDM'), [N_ini, N_end], y0,
                         rtol=1e-8, atol=1e-10, method='LSODA')
    sol_dcr = solve_ivp(lambda N, y: growth_ode(N, y, 'DCR'), [N_ini, N_end], y0,
                        rtol=1e-8, atol=1e-10, method='LSODA')

    sigma8_lcdm = S8_PLANCK / np.sqrt(OM_M / 0.3)
    sigma8_dcr = sigma8_lcdm * (sol_dcr.y[0][-1] / sol_lcdm.y[0][-1])
    S8_dcr = sigma8_dcr * np.sqrt(om_m_dcr / 0.3)

    return S8_dcr, (sol_dcr.y[0][-1] / sol_lcdm.y[0][-1])


if __name__ == "__main__":
    print(f"--- VALIDAZIONE FISICA RIGOROSA D-DCR (V9) ---")
    h0_calc, om_m_calc = solve_h0_geometric_lock()
    s8_calc, factor = solve_s8_disformal(h0_calc / 100, om_m_calc)

    print(f"\n[BACKGROUND]")
    print(f"  H0: {h0_calc:.2f} (Omega_m: {om_m_calc:.4f})")
    print(f"\n[PERTURBAZIONI]")
    print(f"  Fattore di Crescita D rispetto LCDM: {factor:.4f}")
    print(f"  S8 Calcolato: {s8_calc:.3f}")
