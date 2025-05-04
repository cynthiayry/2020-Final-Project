# Causal Inference with BCF, BART, and Horseshoe GLM

This repository contains the final project for **DATA2020** at Brown University. The goal of this project is to estimate heterogeneous treatment effects (HTEs) under confounding using:

- **Bayesian Causal Forest (BCF)**
- **Bayesian Additive Regression Trees (BART)**
- **Generalized Linear Model with Horseshoe Prior (GLM-HS)**

We benchmark these methods on simulated datasets that vary in sample size, treatment effect heterogeneity, and prognostic function nonlinearity.

---

## Key Findings
- **BCF** consistently achieves the lowest RMSE and most reliable coverage across settings.
- **Horseshoe GLM** underperforms on ATE estimates due to model misfit but provides interpretable coefficient-level inference.
- **X₅** is selected as a significant variable in every setting, while **X₁** is rarely selected.
- Larger sample sizes improve coverage, interval length, and variable detection accuracy.

---
