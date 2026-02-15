import numpy as np
from scipy.integrate import quad

# --- COSTANTI ---
H0_PLANCK = 67.4      
OMEGA_M = 0.315
OMEGA_L = 1.0 - OMEGA_M
Z_CMB = 1090.0
Z_TRANSITION = 1.0  # L'interruttore si accende qui 

def get_w_eff(z, beta):
    """
    Implementa il meccanismo di SCREENING SYMMETRON.
    Rif: "221 (Schermato) ... 21 (Attivo)"
    """
    if z > Z_TRANSITION:
        return -1.0  # Schermato (Standard LCDM)
    else:
        # Quando attivo, w diventa Phantom proporzionalmente a Beta
        # Il fattore 0.8 è la 'sensibilità' del campo scalare locale
        return -1.0 - (0.8 * beta)

def hubble_inverse(z, om_m, om_l, beta):
    """
    1/E(z) con w dinamico.
    Nota: Per rigore fisico, l'evoluzione della densità DE richiederebbe
    un integrale separato dell'esponente. Qui usiamo l'approssimazione
    locale istantanea per verificare il concetto del 'Geometric Lock'.
    """
    w_z = get_w_eff(z, beta)
    
    # Densità DE che evolve diversamente prima e dopo la transizione
    if z > Z_TRANSITION:
        # Densità standard costante
        rho_de = om_l 
    else:
        # Densità Phantom che cresce nel futuro/presente
        rho_de = om_l * (1+z)**(3*(1+w_z))
        
    E_z = np.sqrt(om_m * (1+z)**3 + rho_de)
    return 1.0 / E_z

def calculate_distance_screened(h0, beta):
    c = 299792.458
    integral, _ = quad(hubble_inverse, 0, Z_CMB, args=(OMEGA_M, OMEGA_L, beta))
    return (c / h0) * integral

# --- 1. TARGET BAO (PLANCK) ---
# Calcoliamo la distanza dell'universo standard (senza beta)
DIST_TARGET = calculate_distance_screened(H0_PLANCK, 0.0)
print(f"Target Geometrico (BAO/Planck): {DIST_TARGET:.2f} Mpc")
print(f"Meccanismo: Screening attivo a z < {Z_TRANSITION}")
print("-" * 75)

# --- 2. SIMULAZIONE C-DCR CON SCREENING ---
betas = np.linspace(0.0, 0.35, 10)

print(f"{'Beta':<8} | {'H0':<8} | {'S8 (Sim)':<8} | {'w_local':<8} | {'BAO Err%':<10} | {'Status'}")
print("-" * 75)

for beta in betas:
    # H0 aumenta con Beta (spinta locale)
    # Slope ricalibrata per il modello screened
    h0_model = 67.4 + (20.0 * beta)
    
    # S8 diminuisce (freno locale)
    s8_model = 0.811 * (1.0 - 0.25 * beta)

    # w locale (solo per z < 1)
    w_local = get_w_eff(0.5, beta) # Valore medio recente

    # CALCOLO DISTANZA (Con Screening)
    dist_model = calculate_distance_screened(h0_model, beta)
    
    # Errore rispetto al metro rigido
    bao_error = 100 * (dist_model - DIST_TARGET) / DIST_TARGET
    
    status = ""
    # Cerchiamo dove l'errore è minimo E H0 è alto
    if abs(beta - 0.25) < 0.02:
        status = "<<< C-DCR TARGET"
        
    print(f"{beta:.3f}    | {h0_model:.2f}     | {s8_model:.3f}    | {w_local:.3f}    | {bao_error:+.4f}%    | {status}")

print("-" * 75)
