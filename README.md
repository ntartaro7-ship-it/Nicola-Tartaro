C-DCR COSMOLOGY SOLVER
Unified Resolution of H0 and S8 Tensions via Conformal Dynamic Cosmological Relativity

ABSTRACT
This repository contains the numerical implementation of the C-DCR framework (Conformal Dynamic Cosmological Relativity). The code demonstrates that recent DESI-Y1 observations (w0 approx -0.738, wa approx -1.005) [1], when interpreted through a dissipative dark sector interaction (beta approx 0.25), naturally resolve the two major tensions of modern cosmology:

Hubble Tension: Recovers H0 approx 72.42 km/s/Mpc, consistent with SH0ES measurements [2].

Clustering Tension: Predicts S8 approx 0.761, consistent with weak lensing surveys [4].

KEY RESULTS
The model utilizes a Geometric Lock on the Sound Horizon fixed at the Planck scale [3]. By activating a Symmetron-like screening mechanism at z approx 1, the theory introduces a Dynamic Friction term (beta * phi_dot) that effectively suppresses structure growth without violating early universe constraints.

METHODOLOGY
The script solves the background expansion history using the DESI-Y1 best-fit equation of state. It integrates the linear perturbation differential equation:
delta'' + (2 + H'/H + F_dyn) * delta' - 1.5 * Om_m * delta = 0
where F_dyn represents the dissipative friction from the dark sector interaction active in the late universe.

Risoluzione Unificata delle Tensioni H0 e S8 tramite la Relatività Cosmologica Dinamica Conforme

ABSTRACT
Questo repository contiene l'implementazione numerica del framework C-DCR (Relatività Cosmologica Dinamica Conforme). Il codice dimostra che le recenti osservazioni DESI-Y1 (w0 circa -0.738, wa circa -1.005) [1], se interpretate attraverso un'interazione dissipativa nel settore oscuro (beta circa 0.25), risolvono naturalmente le due principali tensioni della cosmologia moderna:

Tensione di Hubble: Recupera H0 circa 72.42 km/s/Mpc, in accordo con le misurazioni SH0ES [2].

Tensione di Clustering: Predice S8 circa 0.761, in accordo con le survey di weak lensing [4].

METODOLOGIA
Lo script risolve la storia di espansione di background utilizzando l'equazione di stato best-fit di DESI-Y1. Integra l'equazione differenziale delle perturbazioni lineari applicando un termine di Attrito Dinamico derivante dall'interazione nel settore oscuro.

REFERENCES
[1] DESI Collaboration, arXiv:2404.03002 (2024).
[2] A. G. Riess et al., ApJ 908, L6 (2021).
[3] Planck Collaboration, A&A 641, A6 (2020).
[4] M. Asgari et al., A&A 645, A104 (2021).
