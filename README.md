# TB_transmission_tPOC-RDT

Overview

This repository hosts a stochastic, individual-based compartmental model of tuberculosis (TB) transmission and progression, developed to evaluate the potential impact of a novel point-of-care rapid diagnostic test (tPOC-RDT) with enhanced rural diagnostic coverage. The model simulates TB dynamics under two scenarios: (1) baseline conditions reflecting current diagnostic access and delays, and (2) implementation of the tPOC-RDT.

The model utilises custom functions to transition individuals between disease states (compartments) based on calculated probabilities. It follows a stepwise method that mirrors the natural history of TB progression, with transitions repeated iteratively over the entire simulation period (in this case, 10 years). At each time step, transition probabilities are dynamically computed based on the group counts from the previous step, allowing the model to capture changes in the population structure over time.

This implementation predominantly uses base R, featuring original R code developed specifically for this project.

Model Description

The model is implemented in R and structured around seven compartments representing key states in TB infection and care:

S – Susceptible
E – Exposed (Latent infection)
P – Preclinical (non-infectious, early disease)
I – Infectious (active TB)
D – Diagnosed
T – Treated
R – Recovered
Individuals transition between these compartments according to biologically informed probabilities and rates, such as infection probability, progression rates, diagnostic rates, treatment uptake, and recovery. Stochastic transitions are modeled using Bernoulli processes applied at each daily time step.

Transition Flows Include:
Susceptible → Exposed via force of infection
Exposed → Preclinical via reactivation
Preclinical → Infectious via symptom progression
Infectious → Diagnosed via passive detection or diagnostic coverage
Diagnosed → Treated
Treated → Recovered
Recovered → Susceptible via loss of immunity or relapse
Exit/re-entry via mortality and birth modeled as a composite re-entry (recidivism) rate
Use Case

specifically below;
#' - inf.prob:         Probability of transmission per contact
#' - act.rate:         Number of contacts per person per time step
#' - piu.rate:         Rate of progression from latent to early infectious (E → P)
#' - delay.rate:       Rate of progression from preclinical to clinical TB (P → I)
#' - p_diag_coverage:  Proportion of individuals in the preclinical stage (P) eligible for diagnosis
#' - p_diag_rate:      Probability of diagnosis among eligible individuals in the preclinical stage (P → D)
#' - i_diag_coverage:  Proportion of individuals in the clinical stage (I) eligible for diagnosis
#' - i_diag_rate:      Probability of diagnosis among eligible individuals in the clinical stage (I → D)
#' - treat.rate:       Probability of treatment uptake once diagnosed (D → T)
#' - rec.rate:         Rate of recovery following treatment (T → R)
#' - recid.rate:       Rate of relapse or TB-related death, balanced by replacement into susceptible pool (R → S)

This model was created to compare baseline conditions to a theoretical implementation of a heat-stable, finger-prick tPOC-RDT with approximately 75–80% sensitivity and 95% specificity. Its major benefit is increased diagnostic coverage, particularly in rural populations with limited access to standard TB diagnostics.

Geographic Context
The model parameters and assumptions are loosely based on epidemiological patterns in the Amazonian region of Brazil, with rural/urban disparities reflected in future model versions.

Development Status

🚧 This model is under active development.
Key next steps include:
Validation and tuning of input parameters using empirical data
Disaggregation of rates by age, gender, and urban vs rural context
Separation of composite flows (e.g., disentangling mortality, relapse, and birth from the current recidivism rate)
Extension to include diagnostic delay stratified by setting and population

This repositary contains 3 files; two coding files that contain baseline and theoretical parameters along with a powerpoint walk through. 
