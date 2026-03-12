# GRB-afterglow-stacking-analysis

[![arXiv](https://img.shields.io/badge/arXiv-2603.09179-b31b1b.svg)](https://arxiv.org/abs/2603.09179)

This repository contains the data and analysis code for the paper:

**"Detection of afterglow emission up to 100 GeV through a stacking analysis of gamma-ray bursts"**  
Authors: [Shi Chen], [Co-author Yuan Qiang Yi-Qing Guo, Ben-Zhong Dai, He Gao, Bing Zhang]  
Time: Tue, 10 Mar 2026 04:33:48 UTC 

## 📁 Repository Structure

- `Data/`: All data files used in the analysis
  - `processed/`: Final stacked SEDs and Light Curves results (machine-readable)
- `Catalogues/`: Gamma-ray Burst Sample Catalogues
- `Code/`: Analysis scripts
  - `stacking/`: GRB stacking pipeline
  - `fitting/`: Spectral fitting codes
  - `plotting/`: Figure generation scripts
- `Model/`: The results of forward shock model

## 🚀 Quick Start

1. Clone this repository:
   ```bash
   git clone https://github.com/shichen144/GRB-afterglow-stacking-analysis.git
   cd GRB-afterglow-stacking-analysis.git
## Software Requirements

This analysis relies on two main software components:

### 1. Fermi Science Tools
The raw Fermi data were processed using the official [Fermi Science Tools](https://fermi.gsfc.nasa.gov/ssc/data/analysis/) (version v11r5p3). 
These tools are required if you wish to reproduce the data reduction steps from the original LAT/LLE files.

### 2. Python Analysis Scripts
All subsequent stacking, fitting, and plotting were performed with Python 3.8+ using standard scientific libraries.

## Teh Forward Shock Model

The afterglow emission in this work is modeled using a **forward shock model** (also known as the external shock model) following the standard framework of GRB afterglow physics. This model calculates the broadband spectra and light curves of GRB afterglows based on the dynamics of a relativistic blast wave propagating into the ambient medium.

### Model Description

The forward shock model assumes that the GRB ejecta drives a relativistic shock into the circumburst medium. Electrons are accelerated at the shock front and cool via synchrotron and synchrotron self-Compton (SSC) radiation. The model includes the following key physical processes:
- Hydrodynamical evolution of the blast wave (adiabatic/radiative phases)
- Electron acceleration with a power-law energy distribution
- Synchrotron emission and self-absorption
- Synchrotron self-Compton (SSC) emission
- Jet break effects when the Lorentz factor drops below 1/θ

### Model Parameters

The model requires the following input parameters:

| Parameter | Symbol | Description |
|-----------|--------|-------------|
| **E₀** | `E0` | Total kinetic energy of the forward shock [erg] |
| **n₀** | `n0` | Ambient medium density at radius R₀ [cm⁻³] |
| **εₑ** | `epsilon_e` | Electron energy equipartition fraction |
| **ε_B** | `epsilon_B` | Magnetic field energy equipartition fraction |
| **p** | `p` | Electron injection spectral index (N(γ) ∝ γ⁻ᵖ) |
| **Γ₀** | `LF0` | Initial Lorentz factor of the ejecta |
| **θ** | `theta` | Jet half-opening angle [rad] |
| **z** | `redshift` | Burst redshift |

### Model Results Provided in This Repository

The forward shock model used in this analysis is **proprietary** and its source code is not publicly available. However, to ensure transparency and reproducibility, we provide the full model outputs used in our paper.
