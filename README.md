# Uncertainty-Aware Survival Modeling of Pavement Deterioration with Connected Vehicle Data

**Author:** Cecil Swanzy, M.Sc. Statistics & Data Science

>  **This is my published Master's thesis.** Full PDF and reproducible R code are included.

---

## Overview

This thesis develops a **survival analysis framework** to predict when pavement sections will require maintenance, with explicit **uncertainty quantification** for decision-makers. Using high-frequency connected vehicle (CV) data from 3,707 Swedish road sections (3,357 deterioration events), I compare **Cox proportional hazards** and **Weibull** models with time-dependent covariates.

**Core contribution:** Moving beyond point estimates to provide maintenance planners with probabilistic predictions (e.g., "median time-to-maintenance is 17.7 weeks with a 95% CI of 16.9–18.7 weeks").

---

## Key Results

| Model | Median Time-to-Maintenance | 95% CI | Key Insight |
|-------|---------------------------|--------|-------------|
| Cox (Static only) | Not estimable | N/A | Relative risks only |
| Cox + CV Indicators | Not estimable | N/A | CV improves fit (LRT p<0.001) |
| **Weibull + CV Indicators** | **17.71 weeks** | **[16.87, 18.65]** | ✅ Decision-ready prediction |

### Key Findings

1. **CV indicators improve predictive performance:** LRT χ²(2) = 34.48, p < 0.001
2. **Model specification matters:** Cox = relative risks; Weibull = calendar-time predictions
3. **Uncertainty is measurable:** Section-level bootstrap (B=300), SD ≈ 0.48 weeks
4. **Hazard decreases over time:** Weibull shape γ = 0.694 < 1

---

## Methodology

### Data
- **3,707 road sections** from Nacka and Östersund, Sweden
- **144 weeks** of observation
- **3,357 deterioration events** (IRI ≥ 3.0 m/km)
- **CV indicators:** iri_mean, variance_mean, r_count_mean

### Models

| Component | Approach |
|-----------|----------|
| Models | Cox PH + Weibull |
| Time-dependent covariates | CV indicators lagged one week |
| Censoring | Right-censoring + delayed entry (counting-process) |
| Uncertainty | Section-level bootstrap (B=300) |
| Comparison | LRT (Cox) + AIC (Weibull) |

### Key Equations

**Cox Model:**
h_k(t) = h₀(t) exp(βᵀX_k^static + θᵀX_k^CV(t))


**Weibull Model:**
h_k(t) = λ_k γ t^(γ-1)
t₀.₅ = (log 2)^(1/γ) exp(βᵀX^static + θᵀX^CV)


---

## Repository Structure
pavement-deterioration-survival-modeling/
├── Executive summary.pdf
├── R code.R # Full R code
├── Uncertainty-Aware Survival Modeling of Pavement Deterioration with Connected Vehicle Data.pdf # Published thesis
└── README.md



## How to Reproduce

```bash
# Clone the repository
git clone https://github.com/CecilSwanzy/pavement-deterioration-survival-modeling.git

# In RStudio:
# 1. Set working directory to project root
# 2. Place data files (nacka.xlsx, ostersund.xlsx) in working directory
# 3. Run:
source("thesis_analysis.R")


Connect with Me
LinkedIn: www.linkedin.com/in/cecil-swanzy-191a60171

Email: mcskanty@gmail.com
