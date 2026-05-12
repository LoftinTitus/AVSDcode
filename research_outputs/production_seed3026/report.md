# AVSD Research Results

## Cohort Summary

```text
n_embryos            : 96.0
prevalence           : 0.1146
closure_fraction     : 0.3125
mean_closure_time    : 3.3549
mean_gap             : 1.5368
median_gap           : 0.4209
std_gap              : 2.7581
mean_linear_score    : 0.842
mean_nonlinear_score : 0.7362
mean_probability     : 0.2759
mean_terminal_shear  : 1.1973
mean_peak_shear      : 1.5322
```

# Calibration

## Objective

Score: 0.007

## Penalties

- `closure_time`: 0.0
- `euploid_prevalence`: 0.0
- `gap`: 0.0047
- `terminal_shear`: 0.0
- `trisomy21_prevalence`: 0.0

## Euploid Summary

```text
n_embryos            : 96.0
prevalence           : 0.0677
closure_fraction     : 0.2865
mean_closure_time    : 3.3585
mean_gap             : 1.4536
median_gap           : 0.3909
std_gap              : 2.6865
mean_linear_score    : 0.8483
mean_nonlinear_score : 0.7389
mean_probability     : 0.2715
mean_terminal_shear  : 1.2075
mean_peak_shear      : 1.5169
```

## Trisomy 21 Summary

```text
n_embryos            : 96.0
prevalence           : 0.276
closure_fraction     : 0.1198
mean_closure_time    : 3.4556
mean_gap             : 0.7032
median_gap           : 0.4572
std_gap              : 1.3832
mean_linear_score    : 0.7749
mean_nonlinear_score : 0.6754
mean_probability     : 0.379
mean_terminal_shear  : 1.2224
mean_peak_shear      : 1.3028
```


# Posterior Summary

Acceptance rate: 0.5311
Mean score: 0.4493
MAP score: 0.007
MAP log posterior: 11.2466
Retained samples: 3804
Burn-in: 250
Max R-hat: 1.0479

## Posterior Means and Intervals

- `alpha_DMP`: 0.6252 [0.428, 0.836] (R-hat 1.0133)
- `alpha_EMT`: 0.8309 [0.5824, 1.0937] (R-hat 1.0038)
- `k_A`: 0.6775 [0.4847, 0.8803] (R-hat 1.0479)
- `k_E`: 0.7016 [0.498, 0.9103] (R-hat 1.0164)
- `k_G`: 0.9628 [0.7411, 1.2156] (R-hat 1.0289)
- `threshold`: 0.6798 [0.5831, 0.7752] (R-hat 1.0255)

## MAP Parameters

- `alpha_DMP`: 0.6
- `alpha_EMT`: 0.8
- `k_A`: 0.68
- `k_E`: 0.72
- `k_G`: 0.95
- `threshold`: 0.7


# Validation

Mean score: 0.6559

## Replicate Scores

- replicate 1: 0.1
- replicate 2: 0.0052
- replicate 3: 0.9492
- replicate 4: 0.0104
- replicate 5: 2.0352
- replicate 6: 0.1115
- replicate 7: 0.0038
- replicate 8: 2.0316

## Pass Rates

- `closure_time`: 1.0
- `euploid_prevalence`: 0.125
- `gap`: 0.375
- `terminal_shear`: 1.0
- `trisomy21_prevalence`: 1.0

## Mean Penalties

- `closure_time`: 0.0
- `euploid_prevalence`: 0.2167
- `gap`: 0.0039
- `terminal_shear`: 0.0
- `trisomy21_prevalence`: 0.0


## Sensitivity

Baseline score: 0.1

### Local Effects

- `alpha_DMP`: -45.249
- `alpha_EMT`: -33.9578
- `k_A`: -15.0708
- `k_G`: -1.6409
- `k_E`: -0.5806
