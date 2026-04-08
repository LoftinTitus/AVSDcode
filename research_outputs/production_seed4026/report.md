# AVSD Research Results

## Cohort Summary

```text
n_embryos            : 96.0
prevalence           : 0.125
closure_fraction     : 0.125
mean_closure_time    : 3.4696
mean_gap             : 0.509
std_gap              : 0.7678
mean_linear_score    : 0.8328
mean_nonlinear_score : 0.7151
mean_probability     : 0.2725
mean_terminal_shear  : 1.2566
mean_peak_shear      : 1.2822
```

# Calibration

## Objective

Score: 0.5079

## Penalties

- `closure_time`: 0.0074
- `euploid_prevalence`: 0.0
- `gap`: 0.3312
- `terminal_shear`: 0.0
- `trisomy21_prevalence`: 0.0

## Euploid Summary

```text
n_embryos            : 96.0
prevalence           : 0.0417
closure_fraction     : 0.0938
mean_closure_time    : 3.4643
mean_gap             : 0.5748
std_gap              : 1.0942
mean_linear_score    : 0.8472
mean_nonlinear_score : 0.7204
mean_probability     : 0.2403
mean_terminal_shear  : 1.2555
mean_peak_shear      : 1.307
```

## Trisomy 21 Summary

```text
n_embryos            : 96.0
prevalence           : 0.2604
closure_fraction     : 0.0417
mean_closure_time    : 3.4858
mean_gap             : 0.525
std_gap              : 0.1608
mean_linear_score    : 0.7599
mean_nonlinear_score : 0.6494
mean_probability     : 0.3701
mean_terminal_shear  : 1.0983
mean_peak_shear      : 1.1144
```


# Posterior Summary

Acceptance rate: 0.625
Mean score: 1.4041
MAP score: 0.5079
MAP log posterior: 10.4177
Retained samples: 1028
Burn-in: 128
Max R-hat: 1.3589

## Posterior Means and Intervals

- `alpha_DMP`: 0.618 [0.4235, 0.7809] (R-hat 1.0484)
- `alpha_EMT`: 0.8004 [0.5825, 1.0266] (R-hat 1.068)
- `k_A`: 0.6437 [0.5305, 0.8497] (R-hat 1.3589)
- `k_E`: 0.7484 [0.548, 0.9081] (R-hat 1.0503)
- `k_G`: 0.8217 [0.6868, 0.936] (R-hat 1.1144)
- `threshold`: 0.6525 [0.5791, 0.7223] (R-hat 1.095)

## MAP Parameters

- `alpha_DMP`: 0.6306
- `alpha_EMT`: 0.7618
- `k_A`: 0.6865
- `k_E`: 0.7169
- `k_G`: 0.805
- `threshold`: 0.6799


# Validation

Mean score: 3.1163

## Replicate Scores

- replicate 1: 1.5557
- replicate 2: 3.8196
- replicate 3: 0.8449
- replicate 4: 8.2195
- replicate 5: 1.1034
- replicate 6: 1.8599
- replicate 7: 7.0987
- replicate 8: 0.429

## Pass Rates

- `closure_time`: 0.0
- `euploid_prevalence`: 0.25
- `gap`: 0.0
- `terminal_shear`: 1.0
- `trisomy21_prevalence`: 0.875

## Mean Penalties

- `closure_time`: 0.0072
- `euploid_prevalence`: 0.6925
- `gap`: 0.6822
- `terminal_shear`: 0.0
- `trisomy21_prevalence`: 0.0012


## Sensitivity

Baseline score: 1.5557

### Local Effects

- `k_G`: 19.1933
- `k_A`: 8.9151
- `alpha_DMP`: -8.455
- `alpha_EMT`: 6.9731
- `k_E`: 2.4805
