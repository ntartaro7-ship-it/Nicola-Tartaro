import numpy as np
import matplotlib.pyplot as plt
from getdist import plots, MCSamples

# ==========================================
# 1. GENERAZIONE MOCK MCMC
# ==========================================
# In una pubblicazione reale, la funzione np.loadtxt()
# sostituirà questa sezione per importare le catene originali.
np.random.seed(42)

# Planck 2018 (LCDM)
mean_planck = [67.4, 0.832]
cov_planck = [[0.5**2, -0.6*0.5*0.013], [-0.6*0.5*0.013, 0.013**2]]
samples_planck = np.random.multivariate_normal(mean_planck, cov_planck, 20000)

# Late Universe (SH0ES + DES)
mean_late = [73.04, 0.759]
cov_late = [[1.04**2, 0], [0, 0.017**2]]
samples_late = np.random.multivariate_normal(mean_late, cov_late, 20000)

# ==========================================
# 2. INIZIALIZZAZIONE PIPELINE GETDIST
# ==========================================
# Nomenclatura compatibile con il rendering LaTeX
names = ['H0', 'S8']
labels = ['H_0', 'S_8']

mc_planck = MCSamples(samples=samples_planck, names=names, labels=labels, label='Planck 2018 ($\\Lambda$CDM)')
mc_late = MCSamples(samples=samples_late, names=names, labels=labels, label='SH0ES + Weak Lensing')

# ==========================================
# 3. PLOTTING PROFESSIONALE
# ==========================================
# getdist gestisce automaticamente i contorni 1-sigma e 2-sigma
g = plots.get_single_plotter(width_inch=8)
g.plot_2d([mc_planck, mc_late], 'H0', 'S8', filled=True, colors=['#1f77b4', '#2ca02c'], alphas=[0.7, 0.7])

# Otteniamo l'asse matplotlib sottostante per le customizzazioni
ax = g.subplots[0, 0]

# Punto teorico D-DCR (Platinum Fit V9)
H0_ddcr = 73.13
S8_ddcr = 0.744

# Tracciamento incrociato e stella rossa
ax.axvline(H0_ddcr, color='darkred', linestyle='--', alpha=0.4, zorder=1)
ax.axhline(S8_ddcr, color='darkred', linestyle='--', alpha=0.4, zorder=1)
ax.scatter(H0_ddcr, S8_ddcr, color='red', marker='*', s=400, zorder=10, edgecolor='black',
           label=f'D-DCR Platinum\n($H_0={H0_ddcr}$, $S_8={S8_ddcr}$)')

# Annotazione della risoluzione
ax.annotate('Tensioni Risolte\n(Unificazione Disformale)',
            xy=(H0_ddcr, S8_ddcr), xycoords='data',
            xytext=(H0_ddcr - 3.5, S8_ddcr + 0.025), textcoords='data',
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2", color='darkred', lw=1.5),
            fontsize=11, color='darkred', fontweight='bold')

# Estetica finale del grafico
ax.set_title("Risoluzione delle Tensioni Cosmologiche", fontsize=16, pad=20)
ax.legend(loc='upper right', fontsize=12, frameon=True, shadow=True)
ax.grid(True, linestyle=':', alpha=0.6)

# Limiti degli assi per inquadrare la rottura di simmetria rispetto a LCDM
ax.set_xlim(65, 76)
ax.set_ylim(0.70, 0.86)

# Watermark
plt.figtext(0.95, 0.05, 'N. Tartaro - D-DCR Analysis (V9)', fontsize=10, color='gray', ha='right', va='bottom', alpha=0.5)

plt.show()
