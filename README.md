PCL Enzymatic Degradation — Reaction-Diffusion Model
Computational analysis code for the manuscript:  
"Reaction-diffusion modelling of enzymatic degradation of semicrystalline PCL films"
---
Overview
This repository contains code for fitting, uncertainty quantification, sensitivity analysis, and thickness variation studies of a 1D reaction-diffusion PDE model describing lipase-mediated degradation of poly(ε-caprolactone) (PCL) films. The model is calibrated against weight loss and DSC crystallinity data from Shi et al. (2020).
The enzyme transport and reaction scheme tracks six species spatially resolved across the film half-thickness: free enzyme (E), crystalline polymer (C), amorphous polymer (A), enzyme–crystalline complex (EC), enzyme–amorphous complex (EA), and degradation product (P). Porosity-dependent diffusion follows the homogenisation result of Wang et al. (2008).
---
Repository Structure
```
.
├── Model_estimation.py   # Least-squares parameter estimation (Python)
├── thickness.py          # Film thickness variation study (Python)
├── DimensionalPlot.m     # Bayesian inference via DE-MCz MCMC (MATLAB)
├── sobolshi.m            # Sobol global sensitivity analysis (MATLAB)
└── plottt.m              # Publication figure: De0 total Sobol index (MATLAB)
```
---
File Descriptions
`Model_estimation.py`
Fits the 8-parameter PDE model to the Shi et al. dataset using weighted nonlinear least squares (`scipy.optimize.least_squares`, trust-region reflective method). Parameters are estimated in log-space to enforce positivity. Outputs R² and reduced chi-squared diagnostics, residual tables as `.xlsx`, and publication-quality fit plots at 600 DPI.
Observables fitted:
Weight loss (WL): spatial mean of polymer solid (C + A), normalised by initial mass
Degree of crystallinity (Xc): spatial mean of C / (C + A)
`thickness.py`
Parametric sweep over five film half-thicknesses (L = 0.1, 0.25, 0.5, 3.0, 4.0 mm) using the calibrated parameter set from `Model_estimation.py`. Computes WL and Xc trajectories for each thickness, reports diffusion time constants (τ = L² / De0), and exports summary plots and an Excel workbook.
`DimensionalPlot.m`
Bayesian posterior inference using the DE-MCz (Differential Evolution Markov Chain with snooker update) algorithm. Runs 5 chains × 8 walkers for 500 iterations (thinning factor 5) in log-parameter space. Produces posterior CI bar charts for all 8 parameters (Figures 9–11) and posterior predictive uncertainty bands for WL and Xc. Includes a full pairwise posterior scatter plot.
`sobolshi.m`
Variance-based global sensitivity analysis using the Saltelli estimator for first-order (Si) and total-order (Ti) Sobol indices. Samples the 8-parameter space via a scrambled Sobol quasi-random sequence (M = 256 base samples; total PDE solves ≈ (p+2)·M = 2560). The model output is the spatially averaged product P(t) evaluated at 12 time points (0–72 h). Results are saved as CSV, Excel, and 600 DPI PNG figures.
`plottt.m`
Standalone publication figure plotting the total-order Sobol index of De0 across five film thicknesses, demonstrating the transition from reaction-limited to diffusion-limited regimes. Exports a square 1200 DPI PNG.
---
Model Parameters
Symbol	Description	Units
k1	Enzyme–crystalline binding rate	h⁻¹
k-1	Enzyme–crystalline unbinding rate	h⁻¹
k3	Enzyme–amorphous binding rate	h⁻¹
k-3	Enzyme–amorphous unbinding rate	h⁻¹
kconv	Crystalline-to-amorphous conversion rate	h⁻¹
kdegC	Crystalline complex degradation rate	h⁻¹
kdegA	Amorphous complex degradation rate	h⁻¹
De0	Reference enzyme diffusion coefficient	mm² h⁻¹
Fixed model constants:
Film half-thickness: L = 0.25 mm (0.5 mm full film, Shi et al.)
Initial crystalline fraction: C₀ = 0.302
Initial amorphous fraction: A₀ = 0.698
Boundary enzyme concentration: E_bulk = 1 (normalised)
Porosity–diffusion coupling: α_E = 1.0 (Wang et al. 2008)
Spatial grid: N = 25 nodes (method of lines)
---
Requirements
Python
Python ≥ 3.9
`numpy`, `scipy`, `matplotlib`, `pandas`, `openpyxl`
Install with:
```bash
pip install numpy scipy matplotlib pandas openpyxl
```
MATLAB
MATLAB R2021a or later
Statistics and Machine Learning Toolbox (for `sobolset`)
No additional toolboxes required for `DimensionalPlot.m` or `plottt.m`
---
Usage
Python scripts
Run from the project directory:
```bash
python Model_estimation.py
python thickness.py
```
All outputs (PNG figures, Excel files) are saved to the current working directory.
MATLAB scripts
From the MATLAB command window:
```matlab
DimensionalPlot      % Bayesian inference + posterior figures
sobolshi             % Sobol sensitivity analysis
plottt               % De0 thickness sensitivity figure
```
All outputs are saved to the system Downloads folder (auto-detected on Windows and macOS/Linux).
---
Outputs
Script	Output files
`Model_estimation.py`	`fit_weight_loss.png`, `fit_crystallinity.png`, `residuals_*.png`, `residuals_at_degradation_times.xlsx`
`thickness.py`	`thickness_WL.png`, `thickness_Xc.png`, `comparison_WL.png`, `comparison_Xc.png`, `thickness_predictions.xlsx`
`DimensionalPlot.m`	`Figure 9.png`, `Figure 10.png`, `Figure 11.png`, `Shi_WL_PredictionUncertainty.png`, `Shi_Xc_PredictionUncertainty.png`, `Shi_PairwisePosterior.png`
`sobolshi.m`	`sobol_first_order_P.png`, `sobol_total_order_P.png`, `sobol_Si_P.csv`, `sobol_Ti_P.csv`, `sobol_indices_P.xlsx`
`plottt.m`	`De0_square_1200dpi.png`
---
Reference Data
Experimental data (weight loss and DSC crystallinity) are from:
> Shi, K., Jing, J., Song, L., Su, T., and Wang, Z. (2020). Enzymatic hydrolysis of polyester: Degradation of poly(ε-caprolactone) by *Candida antarctica* lipase and *Fusarium solani* cutinase. *International Journal of Biological Macromolecules*, 144, 183–189. https://doi.org/10.1016/j.ijbiomac.2019.12.105
Data are embedded directly in the scripts as numeric arrays.
---
Authors
Nanshin Nansak, Leo Creedon, Denis O'Mahoney, Ramen Ghosh, Marion McAfee  
Centre for Mathematical Modelling and Intelligent Systems for Health and Environment (MISHE)  
Atlantic Technological University Sligo
