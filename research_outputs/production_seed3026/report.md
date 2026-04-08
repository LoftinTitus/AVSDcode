# AVSD Research Results

## Cohort Summary

```text
n_embryos            : 96.0
prevalence           : 0.0938
closure_fraction     : 0.1458
mean_closure_time    : 3.4495
mean_gap             : 0.7415
std_gap              : 1.522
mean_linear_score    : 0.8419
mean_nonlinear_score : 0.7358
mean_probability     : 0.2654
mean_terminal_shear  : 1.2399
mean_peak_shear      : 1.3335
```

# Calibration

## Objective

Score: 0.2851

## Penalties

- `closure_time`: 0.0077
- `euploid_prevalence`: 0.0
- `gap`: 0.1824
- `terminal_shear`: 0.0
- `trisomy21_prevalence`: 0.0

## Euploid Summary

```text
n_embryos            : 96.0
prevalence           : 0.0417
closure_fraction     : 0.0521
mean_closure_time    : 3.4761
mean_gap             : 0.4991
std_gap              : 0.7784
mean_linear_score    : 0.8481
mean_nonlinear_score : 0.735
mean_probability     : 0.2588
mean_terminal_shear  : 1.2542
mean_peak_shear      : 1.2802
```

## Trisomy 21 Summary

```text
n_embryos            : 96.0
prevalence           : 0.2292
closure_fraction     : 0.0417
mean_closure_time    : 3.4862
mean_gap             : 0.4998
std_gap              : 0.1805
mean_linear_score    : 0.801
mean_nonlinear_score : 0.6941
mean_probability     : 0.3298
mean_terminal_shear  : 1.1554
mean_peak_shear      : 1.1702
```


# Posterior Summary

Acceptance rate: 0.5703
Mean score: 1.514
MAP score: 0.2851
MAP log posterior: 10.3611
Retained samples: 1028
Burn-in: 128
Max R-hat: 1.4169

## Posterior Means and Intervals

- `alpha_DMP`: 0.5319 [0.4368, 0.6731] (R-hat 1.4169)
- `alpha_EMT`: 0.8274 [0.5094, 1.162] (R-hat 1.2351)
- `k_A`: 0.6787 [0.5066, 0.8554] (R-hat 1.3533)
- `k_E`: 0.7178 [0.6291, 0.906] (R-hat 1.1049)
- `k_G`: 0.7941 [0.6931, 0.9241] (R-hat 1.0866)
- `threshold`: 0.643 [0.5838, 0.702] (R-hat 1.3596)

## MAP Parameters

- `alpha_DMP`: 0.542
- `alpha_EMT`: 0.864
- `k_A`: 0.7332
- `k_E`: 0.7043
- `k_G`: 0.7642
- `threshold`: 0.6923


# Validation

Mean score: 3.7433

## Replicate Scores

- replicate 1: 2.9457
- replicate 2: 1.9834
- replicate 3: 3.8431
- replicate 4: 0.617
- replicate 5: 10.0004
- replicate 6: 2.0019
- replicate 7: 1.4295
- replicate 8: 7.1255

## Pass Rates

- `closure_time`: 0.0
- `euploid_prevalence`: 0.125
- `gap`: 0.0
- `terminal_shear`: 1.0
- `trisomy21_prevalence`: 1.0

## Mean Penalties

- `closure_time`: 0.0072
- `euploid_prevalence`: 0.8869
- `gap`: 0.7144
- `terminal_shear`: 0.0
- `trisomy21_prevalence`: 0.0


## Sensitivity

Baseline score: 2.9457

### Local Effects

- `alpha_DMP`: -88.0962
- `alpha_EMT`: -69.6067
- `k_A`: -35.4714
- `k_G`: 13.8502
- `k_E`: -2.0121
