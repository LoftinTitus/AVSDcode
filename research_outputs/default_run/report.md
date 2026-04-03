# AVSD Research Results

## Cohort Summary

```text
n_embryos            : 96.0
prevalence           : 0.0417
closure_fraction     : 0.0417
mean_closure_time    : 3.485
mean_gap             : 0.4838
std_gap              : 0.1834
mean_linear_score    : 0.6504
mean_nonlinear_score : 0.5865
mean_probability     : 0.2247
mean_terminal_shear  : 1.1805
mean_peak_shear      : 1.1832
```

# Calibration

## Objective

Score: 3.9553

## Penalties

- `closure_time`: 0.0085
- `euploid_prevalence`: 0.0
- `gap`: 0.3192
- `terminal_shear`: 0.0
- `trisomy21_prevalence`: 0.8659

## Euploid Summary

```text
n_embryos            : 96.0
prevalence           : 0.0312
closure_fraction     : 0.0312
mean_closure_time    : 3.4935
mean_gap             : 0.4992
std_gap              : 0.1593
mean_linear_score    : 0.6392
mean_nonlinear_score : 0.5764
mean_probability     : 0.2084
mean_terminal_shear  : 1.1574
mean_peak_shear      : 1.161
```

## Trisomy 21 Summary

```text
n_embryos            : 96.0
prevalence           : 0.0104
closure_fraction     : 0.0104
mean_closure_time    : 3.4978
mean_gap             : 0.5865
std_gap              : 0.154
mean_linear_score    : 0.6011
mean_nonlinear_score : 0.5418
mean_probability     : 0.1615
mean_terminal_shear  : 1.0252
mean_peak_shear      : 1.0298
```


# Posterior Summary

Acceptance rate: 0.2031
Mean score: Inf
MAP score: 3.9553
MAP log posterior: 3.8186
Retained samples: 68
Burn-in: 8
Max R-hat: NaN

## Posterior Means and Intervals

- `alpha_DMP`: 0.516 [0.3048, 0.63] (R-hat 1.6213)
- `alpha_EMT`: 0.6415 [0.366, 0.8] (R-hat 1.2647)
- `k_A`: 0.5414 [0.4291, 0.68] (R-hat 1.1934)
- `k_E`: 0.7784 [0.6144, 1.2474] (R-hat 1.7678)
- `k_G`: 0.9166 [0.7739, 1.209] (R-hat 1.3088)
- `threshold`: 0.8363 [0.7, 1.0526] (R-hat 1.0)

## MAP Parameters

- `alpha_DMP`: 0.3277
- `alpha_EMT`: 0.7152
- `k_A`: 0.49
- `k_E`: 0.6204
- `k_G`: 1.0096
- `threshold`: 0.8228


# Validation

Mean score: Inf

## Replicate Scores

- replicate 1: Inf
- replicate 2: 2.854032102831901e270
- replicate 3: 4.3742
- replicate 4: 3.268
- replicate 5: 2.6529845491393967e163

## Pass Rates

- `closure_time`: 0.0
- `euploid_prevalence`: 0.6
- `gap`: 0.0
- `terminal_shear`: 1.0
- `trisomy21_prevalence`: 0.0

## Mean Penalties

- `closure_time`: 0.0082
- `euploid_prevalence`: 0.0128
- `gap`: Inf
- `terminal_shear`: 0.0
- `trisomy21_prevalence`: 0.8468


## Sensitivity

Baseline score: Inf

### Local Effects

- `alpha_EMT`: Inf
- `k_A`: 8.46408424477275e191
- `k_E`: -20378.2595
- `alpha_DMP`: -16.6399
- `k_G`: -3.9323
