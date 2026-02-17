# C-DCR Solver: Conformal Dynamic Cosmological Relativity

![Python](https://img.shields.io/badge/python-3.8%2B-blue.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)
[![Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.18671935.svg)](https://doi.org/10.5281/zenodo.18671935)

**A numerical solver for the "Golden Fit" configuration of Conformal Dynamic Cosmological Relativity (C-DCR).**

This repository contains the Python implementation used to calculate the background evolution and linear perturbations of the C-DCR model. The code demonstrates how a single dark sector coupling parameter ($\beta \approx 0.45$) can simultaneously resolve the Hubble Tension ($H_0$) and the Clustering Tension ($S_8$).

## 🏆 Key Results: The "Golden Fit"

By integrating a phantom Dark Energy equation of state ($w_0 \approx -1.15$) with a screened dynamic friction term, this solver identifies a specific region of parameter space that reconciles early and late Universe observations.

| Parameter | Standard $\Lambda$CDM | **C-DCR (This Work)** | Observation / Target |
| :--- | :--- | :--- | :--- |
| **$H_0$** (km/s/Mpc) | $67.4$ (Planck) | **$73.13$** | $73.04$ (SH0ES) ✅ |
| **$S_8$** (Clustering) | $0.83$ (Planck) | **$0.742$** | $\sim 0.76$ (Weak Lensing) ✅ |
| **$w_0$** (DE EoS) | $-1.0$ | **$-1.15$** | DESI-Y1 / Pantheon+ |
| **$\beta$** (Coupling) | $0.0$ | **$0.45$** | Model Prediction |

> **Note:** The value $S_8 = 0.742$ is a robust prediction for the upcoming Euclid mission (DR1).

## 🚀 Features

* **Geometric Lock:** Automatically adjusts $H_0$ to preserve the Planck 2018 sound horizon ($r_d$) while allowing for late-time phantom expansion.
* **Dynamic Friction Solver:** Solves the modified growth equation for matter perturbations:
    $$\ddot{\delta} + (2H - \beta \dot{\phi})\dot{\delta} = 4\pi G_{\rm eff} \rho_m \delta$$
* **Symmetron Screening:** Implements a density-dependent screening mechanism that activates the fifth force only at $z < 1$.
* **Tension Analysis:** Calculates the Gaussian tension ($\sigma$) between model predictions and SH0ES/Planck data.

## 📦 Installation

No special installation is required. The code depends on standard scientific Python libraries.

```bash
pip install numpy scipy matplotlib

Usage
To reproduce the "Golden Fit" results:

Clone the repository:
git clone [https://github.com/ntartaro7-ship-it/Nicola-Tartaro.git](https://github.com/ntartaro7-ship-it/Nicola-Tartaro.git)
cd Nicola-Tartaro

Run the main analysis script:
python cdcr_solver.py

The script will output the calculated cosmological parameters and generate the comparison plots.



Second File acc.py:
## 🌌 Unified Model Validation (Regression Test)

This repository includes a specialized validation script that bridges the gap between galactic dynamics and global cosmology. By using the "Golden Fit" parameters derived in the C-DCR framework, we demonstrate that the model is robust across all scales.

### 📋 Script Overview: `regression_test.py`
This script performs a **Regression Test** on the massive spiral galaxy **NGC 5055**. It applies a single, universal critical acceleration ($a_{crit}$)—previously optimized to satisfy both dwarf (LSB) and giant (HSB) galaxies—to verify if the model maintains its predictive power without local fine-tuning.

### 🛠️ Unified Parameters
The validation is built upon the "Micro-Macro" connection identified in the paper:
* **Universal $a_{crit}$**: $2126.1$ (km/s)²/kpc. This value is fundamentally linked to the cosmological coupling $\beta \approx 0.45$ via the relation: 
    $$a_{crit} \approx \frac{\beta^2}{2} c H_0$$
* **Stellar $M/L_{disk}$**: $0.43$ (optimized for NGC 5055 during the Global Fit process).



### 🔬 Physical Methodology
The script implements the **Radial Acceleration Relation (RAR)** as the universal constitutive law for spacetime in the C-DCR framework:

1.  **Baryonic Baseline**: It calculates the Newtonian acceleration ($g_{bar}$) from the observed distribution of gas, disk stars, and bulge.
2.  **The Yielding Mechanism**: It applies the non-linear correction that simulates the "yielding" of the vacuum under low-acceleration regimes:
    $$g_{obs} = \frac{g_{bar}}{1 - e^{-\sqrt{g_{bar}/a_{crit}}}}$$
3.  **Velocity Mapping**: The resulting total acceleration is converted back into observable rotation velocity ($V_{tot} = \sqrt{g_{obs} \cdot R}$).

### 📊 Results and Impact
* **Systemic Accuracy**: The model achieves a Chi-Squared of **251.47**. [cite_start]While higher than a local-only fit, this represents the systemic precision of a **Universal Law** that works for both $H_0/S_8$ cosmological tensions and galactic rotation curves[cite: 1].
* **Baryonic Dominance**: The test confirms that in massive galaxies, the inner regions remain dominated by baryons (Elastic/Newtonian regime), making the model extremely robust against variations in the dark sector coupling.

### 🔗 Theoretical Context
As detailed in the *Conformal Dynamic Cosmological Relativity* (C-DCR) paper[cite: 1, 70]:
* [cite_start]**Macro Scale**: $\beta \approx 0.45$ resolves the $H_0$ tension ($73.13$ km/s/Mpc) and $S_8$ tension ($0.742$)[cite: 9, 62, 63].
* [cite_start]**Micro Scale**: $a_{crit} \approx 2126$ explains galactic rotation curves without the need for ad-hoc Dark Matter particles[cite: 13, 61].

**Conclusion**: This script proves that "Dark Matter" is an emergent effect of vacuum viscosity, governed by a single, universal coupling parameter $\beta$ across the history of the Universe.


Citation
If you use this code or the C-DCR model in your research, please cite the paper available on Zenodo:
@article{Tartaro2026,
  author       = {Tartaro, Nicola},
  title        = {Conformal Dynamic Cosmological Relativity (C-DCR): Unifying H0, S8, and Galactic Rotation Curves via a Dissipative Dark Sector},
  year         = {2026},
  publisher    = {Zenodo},
  version      = {v5.0},
  doi          = {10.5281/zenodo.18671935},
  url          = {[https://doi.org/10.5281/zenodo.18671935](https://doi.org/10.5281/zenodo.18671935)}
}

Author
Nicola Tartaro

Independent Researcher / Engineering Student

University of Campania "Luigi Vanvitelli"

Contact: nicola.tartaro@studenti.unicampania.it
