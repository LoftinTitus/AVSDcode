# AVSD Dynamical Systems Model (Codex Spec)

## Overview

This project implements a mechanistic, stochastic dynamical systems model for atrioventricular septal defect (AVSD) development.

Core ideas:
- AVSD is modeled as a threshold-crossing failure
- Development is governed by coupled stochastic differential equations (SDEs)
- Genetics influences parameters, not outcomes directly
- Outcome depends on whether the system reaches sufficient structural closure by time T

---

## Model Structure

### State Vector

u(t) = [C(t), D(t), G(t), E(t), A(t)]

| Variable | Meaning |
|----------|--------|
| C(t) | Cushion mesenchymal mass |
| D(t) | DMP contribution |
| G(t) | Septation gap (remaining) |
| E(t) | ECM maturation (0–1) |
| A(t) | AV canal competence (0–1) |

---

## Variable Transformations

To enforce constraints:

x_C = log(C)  
x_D = log(D)  
x_G = log(G)  

x_E = log(E / (1 - E))  
x_A = log(A / (1 - A))  

Recover:

C = exp(x_C)  
D = exp(x_D)  
G = exp(x_G)  

E = 1 / (1 + exp(-x_E))  
A = 1 / (1 + exp(-x_A))  

---

## Time Domain

t ∈ [0, T], where T ≈ 3.5 days

- t = 0 → E9.0 (mouse)
- Key windows:
  - EndMT: ~E9.5–E10.5
  - DMP: ~E9.5–E10.5
  - Closure: ~E10.5–E12.5

---

## Core Dynamics (SDE System)

General form:

du = f(u, t; θ) dt + Σ(u, t; θ) dW

---

### Cushion Growth

dC = [
  α_EMT * W_EMT(t)
  + r_C * C * φ_A(A) * φ_E(E)
  - δ_C * C
] dt + σ_C * C * dW_C

---

### DMP Contribution

dD = [
  α_DMP * W_DMP(t) * φ_A(A) * P_late(t)
  + r_D * D
  - δ_D * D
] dt + σ_D * D * dW_D

---

### Gap Closure

dG = [
  -k_G * W_close(t)
  * (ω_C * c + ω_D * d)
  * (ω_E * E + ω_A * A)
  + k_N * χ(N)
] dt + σ_G * dW_G

---

### ECM Maturation

dE = [
  k_E * (1 - E) * ψ_N(N) * ψ_C(c)
  - γ_E * E
] dt + σ_E * dW_E

---

### AV Canal Competence

dA = [
  k_A * W_AVC(t) * (1 - A)
  - γ_A * A
] dt + σ_A * dW_A

---

##  Supporting Functions

### Timing Window

W(t; t_on, t_off, Δ) =
  σ((t - t_on)/Δ) * σ((t_off - t)/Δ)

σ(x) = 1 / (1 + exp(-x))

---

### Late Penalty

P_late(t; t_dead, λ) =
  exp(-λ * max(t - t_dead, 0))

---

### Activation Functions

φ_A(A) = A^n / (A^n + K^n)  
φ_E(E) = E^n / (E^n + K^n)  

ψ_N(N) = N^n / (N^n + K^n)  
ψ_C(c) = c^n / (c^n + K^n)  

χ(N) = max(N - N_crit, 0)

---

##  Hemodynamics

Q(t) = Q0 * exp(k_Q * t)  

R_eff(t) = R0 + κ_R * g(t)  

τ(t) = (4 * μ * Q(t)) / (π * R_eff(t)^3)  

N(t) = τ(t) / τ*

---

##  Genotype → Parameter Mapping

log(θ_j) =
  log(θ_j^0)
  + β_j^(T21)
  + Σ B_jm * z_m
  + ε_j

Where:
- z_m = latent genetic axes
- β_j^(T21) = trisomy effect
- ε_j = stochastic variation

---

## Phenotype Definition

### Linear Model

S_lin = w_C * c + w_D * d + w_E * E + w_A * A

---

### Nonlinear Model

S_nonlin =
  c^(w_C) *
  d^(w_D) *
  E^(w_E) *
  A^(w_A)

---

### Threshold Rule

Y = 1 if S(T) ≥ τ else 0

Probabilistic version:

P(Y=1) = 1 / (1 + exp(-k * (S - τ)))

---

##  Simulation Pipeline

For each embryo:

1. Sample genetic modifiers z
2. Map to parameters θ
3. Solve SDE system
4. Recover physical variables
5. Compute:
   - G(T)
   - S_lin(T), S_nonlin(T)
   - phenotype Y

---

##  Population Simulation

- Sample N ≈ 1e3 – 1e5 embryos
- Run simulations independently (parallelizable)

Outputs:
- phenotype prevalence
- trajectory distributions
- sensitivity analysis

---

##  Numerical Methods

### Deterministic Phase

- Use stiff ODE solvers:
  - BDF
  - implicit Runge–Kutta

### Stochastic Phase

- Use SDE solvers:
  - Euler–Maruyama (baseline)
  - stochastic Runge–Kutta (preferred)

---

## Calibration Targets

Model should reproduce:

- AVSD prevalence:
  - Trisomy 21: ~15–33%
  - Euploid: low baseline
- Correct developmental timing
- Realistic trajectories
- Physiological shear stress ranges

---

##  Suggested Code Structure

This repository is implemented in Julia and follows the same conceptual split:

/src
  AVSDModel.jl
  /model
    state.jl
    dynamics.jl
    hemodynamics.jl
    activation.jl
    genotype.jl
  /simulation
    solver.jl
    sde_integrator.jl
    ensemble.jl
  /analysis
    phenotype.jl
    metrics.jl
    sensitivity.jl
  /config
    config.jl

/config
  parameters.yaml
  priors.yaml

/test
  runtests.jl

##  Julia Quick Start

1. `julia --project=.`
2. `using Pkg; Pkg.test()`
3. `include("examples/basic_run.jl")`

---

##  Design Principles

- Distributed failure, not single-cause
- Timing matters as much as magnitude
- Nonlinear interactions likely dominate
- Genetics shifts risk landscape, not deterministic outcomes
