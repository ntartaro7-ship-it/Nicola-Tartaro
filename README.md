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
