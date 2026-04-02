# PharmaTrialSim 🧬

A fully reproducible, cloud-deployed oncology clinical trial simulation
and survival analysis pipeline built in R, with AWS S3 data storage and
an interactive Shiny dashboard hosted on AWS EC2.

Built to demonstrate end-to-end pharmaceutical data science workflows:
causal inference under confounding, correlated failure time modeling,
and cloud deployment — all in a realistic multi-site trial setting.

---

## Live Dashboard

> http://<your-ec2-public-ip>:3838/pharmatrialsim

---

## Scientific Background

Real Phase III oncology trials present two core statistical challenges:

**1. Treatment Selection Bias**
In observational settings, patients are not randomly assigned to
treatment. Sicker patients (higher ECOG score, unfavorable biomarker
profile) may be systematically less likely to receive treatment,
confounding naive survival estimates. This project corrects for
confounding using **Inverse Probability Weighting (IPW)** via
propensity score estimation.

**2. Correlated Failure Times**
Patients enrolled at the same clinical site share unmeasured
environmental, demographic, and care-quality factors that induce
within-cluster correlation in their survival times. Ignoring this
leads to underestimated standard errors and anti-conservative inference.
This project addresses within-site correlation using a
**shared gamma frailty Cox model**.

---

## Data Simulation

- **n = 2,000** patients across **20 clinical sites** (100 per site)
- Site-level frailty drawn from a log-normal distribution (variance = 0.4)
- Covariates: age (Normal, mean=60), ECOG status (0/1/2),
  continuous biomarker expression
- Treatment assignment is **non-random**: depends on biomarker and
  ECOG status via a logistic propensity model — creating realistic
  confounding
- Survival times simulated from a **Weibull distribution** with:
  - True treatment HR = 0.65 (protective effect)
  - Site frailty multiplicatively scaling the baseline hazard
  - Max follow-up: 36 months

---

## Statistical Methods

### Propensity Score IPW (`WeightIt`)
```r
weightit(treatment ~ age + ecog + biomarker,
         method = "ps", estimand = "ATE")
```
ATE weights rebalance the sample to resemble a randomized trial,
removing confounding by age, ECOG status, and biomarker expression.

### IPW-Weighted Kaplan-Meier (`survival`)
```r
survfit(Surv(time, event) ~ treatment,
        data = df, weights = df$ipw)
```
Produces covariate-adjusted survival curves for treatment vs. control.

### Robust Weighted Cox Model (`survival`)
```r
coxph(Surv(time, event) ~ treatment + age + ecog + biomarker,
      weights = df$ipw, robust = TRUE)
```
Sandwich standard errors account for IPW weighting. Estimates adjusted
hazard ratios for all covariates.

### Shared Gamma Frailty Cox Model (`survival`)
```r
coxph(Surv(time, event) ~ treatment + age + ecog + biomarker +
        frailty(site, distribution = "gamma"))
```
Models within-site correlation via a site-level random effect.
Site-specific frailty estimates are extracted and visualized.

---

## AWS Architecture
```
Local R Scripts
      │
      ▼
AWS S3 (pharmatrialsim-bucket)
  ├── trial_data.csv       ← simulate_trial.R output
  └── results.rds          ← analyze.R output (KM, Cox, frailty, data)
      │
      ▼
AWS EC2 (t2.micro, Ubuntu 22.04)
  └── Shiny Server → app.R reads results.rds from S3 on startup
```

---

## Shiny Dashboard

Four interactive tabs:

| Tab | Content |
|---|---|
| **Kaplan-Meier** | IPW-adjusted KM curves with CIs and risk table; filterable by age range and ECOG status |
| **Cox Model** | HR estimates with 95% CI and robust p-values |
| **Covariate Balance** | IPW-weighted density plots for age, ECOG, and biomarker by treatment group |
| **Frailty Model** | Frailty-adjusted HR table + bar chart of site-level frailty estimates |

---

## Project Structure
```
PharmaTrialSim/
├── data_gen/
│   └── simulate_trial.R     # Simulate patients, upload to S3
├── analysis/
│   └── analyze.R            # IPW + KM + Cox + frailty, upload to S3
├── app/
│   └── app.R                # Shiny dashboard, reads from S3
└── README.md
```

---

## Requirements

### R Packages
```r
install.packages(c(
  "tidyverse", "simsurv", "survival", "survminer",
  "WeightIt", "aws.s3", "bslib", "broom", "shiny"
))
```

### AWS
- S3 bucket named `pharmatrialsim-bucket` (us-east-1)
- EC2 instance (t2.micro) with port 3838 open
- AWS credentials stored in `~/.Renviron`:
```
AWS_ACCESS_KEY_ID=your_key
AWS_SECRET_ACCESS_KEY=your_secret
AWS_DEFAULT_REGION=us-east-1
```

---

## Running Locally
```r
# Step 1 — Simulate data and upload to S3
source("data_gen/simulate_trial.R")

# Step 2 — Run analysis and upload results to S3
source("analysis/analyze.R")

# Step 3 — Launch dashboard locally
shiny::runApp("app/app.R")
```

---

## Deploying to EC2
```bash
# SSH into your instance
ssh -i your-key.pem ubuntu@<ec2-public-ip>

# Install R + Shiny Server (see deploy instructions)
# Copy app files
scp -i your-key.pem -r app/* ubuntu@<ec2-ip>:/srv/shiny-server/pharmatrialsim/

# Restart Shiny Server
sudo systemctl restart shiny-server
```

Full deployment walkthrough in [`deploy/DEPLOY.md`](deploy/DEPLOY.md).

---

## Author

**Ronit Gandhi**
Ph.D. Candidate in Biostatistics, University of Nebraska Medical Center
[github.com/ronit-gandhi](https://github.com/ronit-gandhi) ·
[linkedin.com/in/ronitg](https://linkedin.com/in/ronitg)
