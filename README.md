<div align="center">

# Extending the Idea of Bayesian Grouped Regression under Logistic Setup

### *Bayesian grouped shrinkage methods for logistic regression in R*

<p>
  <img src="https://img.shields.io/badge/Language-R-276DC3?style=for-the-badge&logo=r&logoColor=white" />
  <img src="https://img.shields.io/badge/Model-Bayesian%20Logistic%20Regression-7B2CBF?style=for-the-badge" />
  <img src="https://img.shields.io/badge/Methods-BGHS%20%7C%20GRASP%20%7C%20GIGG-C1121F?style=for-the-badge" />
  <img src="https://img.shields.io/badge/MCMC-Gibbs%20%2B%20MHwG-2A9D8F?style=for-the-badge" />
</p>

### Author

<h3 align="center" style="color:#7B2CBF; font-size:32px; margin-bottom:6px;">
  Tathastu Ghosh
</h3>

<p align="center">
  <em>M.Sc. in Statistics</em><br>
  Department of Statistics, University of Calcutta
</p>

</div>

---

## Overview

This repository contains the **R implementation** developed for my Master’s dissertation on extending **Bayesian grouped shrinkage regression** to the **logistic regression** setting.

The project studies grouped Bayesian shrinkage methods for **binary-response regression** using **Pólya–Gamma data augmentation** and model-specific posterior computation strategies.

### Methods implemented

* **BGHS** — Bayesian Grouped Horseshoe
* **GRASP** — Grouped Regression with Adaptive Shrinkage Priors
* **GIGG** — Group Inverse-Gamma Gamma Shrinkage

---

## Highlights

* Bayesian grouped shrinkage under **logistic regression**
* **Pólya–Gamma augmentation** for binary outcomes
* Posterior simulation using **Gibbs sampling**
* **Metropolis–Hastings within Gibbs (MHwG)** updates where required
* Code for **proposal tuning**, **trace plots**, and **ACF diagnostics**
* Developed as the computational component of a **Master’s dissertation in Statistics**

---

## Repository Structure

```text
Extending-the-idea-of-Bayesian-Grouped-Regression-under-Logistic-Setup/
│
├── Code files/
│   ├── BGHS implementation
│   ├── GRASP implementation
│   ├── GIGG implementation
│
├── Trace plots and ACF plots/
│   └── MCMC diagnostic scripts / outputs
│
└── README.md
```

---

## Repository Contents

### 📁 `Code files/`

Contains the main **R implementations** of the posterior sampling algorithms for:

* **BGHS**
* **GRASP**
* **GIGG**

These scripts form the core computational component of the dissertation.

### ⚙️ `Tuning of proposal distribution under MHwG/`

Contains code related to **proposal tuning** for the **Metropolis–Hastings within Gibbs** steps used in models where certain full conditional distributions are not available in closed form.

### 📈 `Trace plots and ACF plots/`

Contains scripts and/or outputs for **MCMC diagnostics**, including:

* trace plots,
* autocorrelation function (ACF) plots,
* visual convergence assessment.

---

## Methods Covered

| Method    | Description                                                      |
| --------- | ---------------------------------------------------------------- |
| **BGHS**  | Bayesian grouped horseshoe model under logistic regression       |
| **GRASP** | Grouped adaptive shrinkage model with MHwG-based updates         |
| **GIGG**  | Group inverse-gamma gamma shrinkage model for grouped predictors |

---

## Computational Framework

The posterior computation in this repository is based on **Markov chain Monte Carlo (MCMC)** methods, combining:

* **Pólya–Gamma latent variable augmentation** for the logistic likelihood,
* **Gibbs sampling** for tractable full conditional distributions,
* **Metropolis–Hastings within Gibbs (MHwG)** updates for parameters requiring non-conjugate sampling steps.

---

## Typical Workflow

1. Open the required script from **`Code files/`**
2. Load the necessary R packages
3. Prepare the dataset in the format expected by the script
4. Specify initial values, tuning values, number of iterations, and burn-in
5. Run the Gibbs / MHwG sampler
6. Use the diagnostic scripts to inspect trace plots and ACF plots

---

## Software and Packages

The code is written in **R**. Depending on the script, the following packages may be required:

* `BayesLogit`
* `MASS`
* `MCMCpack`
* `coda`
* `ggplot2`

---

## Project Context

This repository serves as the **computational companion** to my dissertation and documents the implementation of grouped Bayesian shrinkage methods under a logistic regression framework.

It is intended to provide:

* the main sampler implementations for **BGHS**, **GRASP**, and **GIGG**,
* supporting code for **proposal tuning**,
* and **diagnostic tools** for posterior simulation assessment.

---

## Notes

This repository is part of an academic dissertation project. The code has been developed for **research and educational purposes** and may continue to be refined as the dissertation evolves.
