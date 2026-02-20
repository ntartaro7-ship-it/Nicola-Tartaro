import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms

def confidence_ellipse(x, y, ax, n_std=1.0, facecolor='none', **kwargs):
    """
    Crea un'ellisse di covarianza per visualizzare 1 sigma e 2 sigma.
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, **kwargs)

    # Calcolo della scala (deviazione standard)
    scale_x = np.sqrt(cov[0, 0]) * n_std
    scale_y = np.sqrt(cov[1, 1]) * n_std
    
    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(np.mean(x), np.mean(y))

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)

# ==========================================
# DATI OSSERVATIVI E DEL MODELLO
# ==========================================

# 1. PLANCK 2018 (LCDM) - Early Universe
np.random.seed(42)
mean_planck = [67.4, 0.832]
cov_planck = [[0.5**2, -0.6*0.5*0.013], [-0.6*0.5*0.013, 0.013**2]] 
x_planck, y_planck = np.random.multivariate_normal(mean_planck, cov_planck, 5000).T

# 2. LATE UNIVERSE (SH0ES + Weak Lensing DES/KiDS)
mean_late = [73.04, 0.759]
cov_late = [[1.04**2, 0], [0, 0.017**2]]
x_late, y_late = np.random.multivariate_normal(mean_late, cov_late, 5000).T

# 3. D-DCR (DISFORMAL PLATINUM) - Risultato rigoroso aggiornato
# H0 = 73.13, S8 = 0.744
mean_ddcr = [73.13, 0.744]

# ==========================================
# PLOTTING
# ==========================================
fig, ax = plt.subplots(figsize=(10, 8))

# --- Disegno Ellissi (1 sigma e 2 sigma) ---
# Planck (Blu)
confidence_ellipse(x_planck, y_planck, ax, n_std=1.0, edgecolor='blue', linewidth=2, label=r'Planck 2018 ($\Lambda$CDM)')
confidence_ellipse(x_planck, y_planck, ax, n_std=2.0, edgecolor='blue', linestyle='--', linewidth=1)
ax.scatter(mean_planck[0], mean_planck[1], color='blue', s=30)

# Late Universe (Verde)
confidence_ellipse(x_late, y_late, ax, n_std=1.0, edgecolor='green', linewidth=2, label=r'SH0ES + Weak Lensing')
confidence_ellipse(x_late, y_late, ax, n_std=2.0, edgecolor='green', linestyle='--', linewidth=1)
ax.scatter(mean_late[0], mean_late[1], color='green', s=30)

# --- Il Tuo Modello (Stella Rossa) ---
ax.scatter(mean_ddcr[0], mean_ddcr[1], color='red', marker='*', s=400, zorder=10, 
           label=f'D-DCR Platinum\n($H_0={mean_ddcr[0]}$, $S_8={mean_ddcr[1]}$)')

# Annotazioni e Frecce
ax.annotate('Tensioni Risolte\n(Unificazione Disformale)', 
             xy=(mean_ddcr[0], mean_ddcr[1]), xycoords='data',
             xytext=(mean_ddcr[0]-2.5, mean_ddcr[1]-0.035), textcoords='data',
             arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2", color='red'),
             fontsize=11, color='darkred', fontweight='bold')

# Setup Grafico
ax.set_title(r"Risoluzione delle Tensioni Cosmologiche ($H_0$ vs $S_8$)", fontsize=16)
ax.set_xlabel(r"$H_0$ [km s$^{-1}$ Mpc$^{-1}$]", fontsize=14)
ax.set_ylabel(r"$S_8 \equiv \sigma_8 \sqrt{\Omega_m/0.3}$", fontsize=14)
ax.legend(loc='upper left', fontsize=12, frameon=True, shadow=True)
ax.grid(True, linestyle=':', alpha=0.6)

# Limiti assi per focalizzare la zona di interesse
ax.set_xlim(65, 76)
ax.set_ylim(0.70, 0.86)

# Watermark 
fig.text(0.95, 0.05, 'N. Tartaro - D-DCR Analysis (V9)', fontsize=10, color='gray', ha='right', va='bottom', alpha=0.5)

plt.tight_layout()
plt.show()
