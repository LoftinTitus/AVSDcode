# AVSD Research Results

## Cohort Summary

```text
n_embryos            : 96.0
prevalence           : 0.0521
closure_fraction     : 0.125
mean_closure_time    : 3.4524
mean_gap             : 0.7267
std_gap              : 1.5117
mean_linear_score    : 0.8191
mean_nonlinear_score : 0.7084
mean_probability     : 0.2475
mean_terminal_shear  : 1.2348
mean_peak_shear      : 1.3281
```

# Calibration

## Objective

Score: 0.5581

## Penalties

- `closure_time`: 0.0073
- `euploid_prevalence`: 0.0
- `gap`: 0.3648
- `terminal_shear`: 0.0
- `trisomy21_prevalence`: 0.0

## Euploid Summary

```text
n_embryos            : 96.0
prevalence           : 0.0312
closure_fraction     : 0.1146
mean_closure_time    : 3.462
mean_gap             : 0.5572
std_gap              : 1.0907
mean_linear_score    : 0.8289
mean_nonlinear_score : 0.7209
mean_probability     : 0.223
mean_terminal_shear  : 1.288
mean_peak_shear      : 1.3417
```

## Trisomy 21 Summary

```text
n_embryos            : 96.0
prevalence           : 0.1771
closure_fraction     : 0.0625
mean_closure_time    : 3.4827
mean_gap             : 0.5655
std_gap              : 0.7822
mean_linear_score    : 0.7684
mean_nonlinear_score : 0.6678
mean_probability     : 0.3163
mean_terminal_shear  : 1.1633
mean_peak_shear      : 1.1906
```


# Posterior Summary

Acceptance rate: 0.5553
Mean score: 1.5761
MAP score: 0.5581
MAP log posterior: 10.218
Retained samples: 420
Burn-in: 32
Max R-hat: 1.5868

## Posterior Means and Intervals

- `alpha_DMP`: 0.6673 [0.4841, 0.8875] (R-hat 1.5007)
- `alpha_EMT`: 0.7722 [0.4529, 0.9443] (R-hat 1.0773)
- `k_A`: 0.5979 [0.4626, 0.7118] (R-hat 1.2521)
- `k_E`: 0.6775 [0.422, 0.917] (R-hat 1.5816)
- `k_G`: 0.8222 [0.6557, 1.0528] (R-hat 1.2448)
- `threshold`: 0.6286 [0.5558, 0.6949] (R-hat 1.0418)

## MAP Parameters

- `alpha_DMP`: 0.5737
- `alpha_EMT`: 0.7619
- `k_A`: 0.7047
- `k_E`: 0.7514
- `k_G`: 0.8187
- `threshold`: 0.6497


# Validation

Mean score: 1.6311

## Replicate Scores

- replicate 1: 1.1921
- replicate 2: 2.0392
- replicate 3: 2.0418
- replicate 4: 2.045
- replicate 5: 0.8372

## Pass Rates

- `closure_time`: 0.0
- `euploid_prevalence`: 0.4
- `gap`: 0.0
- `terminal_shear`: 1.0
- `trisomy21_prevalence`: 0.2

## Mean Penalties

- `closure_time`: 0.0071
- `euploid_prevalence`: 0.0253
- `gap`: 0.8294
- `terminal_shear`: 0.0
- `trisomy21_prevalence`: 0.0751


## Sensitivity

Baseline score: 1.1921

### Local Effects

- `k_A`: 9.1788
- `alpha_EMT`: 7.7187
- `alpha_DMP`: 6.8442
- `k_G`: 5.6692
- `k_E`: 0.7633
