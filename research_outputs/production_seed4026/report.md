# AVSD Research Results

## Cohort Summary

```text
n_embryos            : 96.0
prevalence           : 0.1667
closure_fraction     : 0.25
mean_closure_time    : 3.3923
mean_gap             : 1.2615
median_gap           : 0.3733
std_gap              : 2.4284
mean_linear_score    : 0.8449
mean_nonlinear_score : 0.7378
mean_probability     : 0.2931
mean_terminal_shear  : 1.2105
mean_peak_shear      : 1.4813
```

# Calibration

## Objective

Score: 0.0004

## Penalties

- `closure_time`: 0.0
- `euploid_prevalence`: 0.0
- `gap`: 0.0002
- `terminal_shear`: 0.0
- `trisomy21_prevalence`: 0.0

## Euploid Summary

```text
n_embryos            : 96.0
prevalence           : 0.0677
closure_fraction     : 0.3125
mean_closure_time    : 3.3602
mean_gap             : 1.3538
median_gap           : 0.3557
std_gap              : 2.617
mean_linear_score    : 0.8568
mean_nonlinear_score : 0.7493
mean_probability     : 0.268
mean_terminal_shear  : 1.2676
mean_peak_shear      : 1.5715
```

## Trisomy 21 Summary

```text
n_embryos            : 96.0
prevalence           : 0.2396
closure_fraction     : 0.1354
mean_closure_time    : 3.4394
mean_gap             : 0.8331
median_gap           : 0.4265
std_gap              : 1.7051
mean_linear_score    : 0.7936
mean_nonlinear_score : 0.6899
mean_probability     : 0.3622
mean_terminal_shear  : 1.2338
mean_peak_shear      : 1.3586
```


# Posterior Summary

Acceptance rate: 0.51
Mean score: 0.5258
MAP score: 0.0004
MAP log posterior: 11.2466
Retained samples: 3804
Burn-in: 250
Max R-hat: 1.0546

## Posterior Means and Intervals

- `alpha_DMP`: 0.6246 [0.4345, 0.8795] (R-hat 1.0055)
- `alpha_EMT`: 0.8218 [0.5994, 1.0947] (R-hat 1.0067)
- `k_A`: 0.7037 [0.5072, 0.9246] (R-hat 1.0456)
- `k_E`: 0.7065 [0.4966, 0.9537] (R-hat 1.0367)
- `k_G`: 0.9343 [0.6684, 1.2072] (R-hat 1.0546)
- `threshold`: 0.6847 [0.5899, 0.7717] (R-hat 1.0013)

## MAP Parameters

- `alpha_DMP`: 0.6082
- `alpha_EMT`: 0.7981
- `k_A`: 0.6898
- `k_E`: 0.7196
- `k_G`: 0.9512
- `threshold`: 0.7066


# Validation

Mean score: 0.9389

## Replicate Scores

- replicate 1: 0.0052
- replicate 2: 0.9492
- replicate 3: 0.0083
- replicate 4: 2.7274
- replicate 5: 0.2926
- replicate 6: 0.0041
- replicate 7: 3.5208
- replicate 8: 0.0039

## Pass Rates

- `closure_time`: 1.0
- `euploid_prevalence`: 0.25
- `gap`: 0.375
- `terminal_shear`: 1.0
- `trisomy21_prevalence`: 1.0

## Mean Penalties

- `closure_time`: 0.0
- `euploid_prevalence`: 0.3116
- `gap`: 0.0027
- `terminal_shear`: 0.0
- `trisomy21_prevalence`: 0.0


## Sensitivity

Baseline score: 0.0052

### Local Effects

- `alpha_DMP`: -4.609
- `k_A`: -4.0637
- `alpha_EMT`: -2.8161
- `k_G`: -0.0356
- `k_E`: 0.0
