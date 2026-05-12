# AVSD Research Results

## Cohort Summary

```text
n_embryos            : 96.0
prevalence           : 0.0938
closure_fraction     : 0.3438
mean_closure_time    : 3.34
mean_gap             : 1.5056
median_gap           : 0.3595
std_gap              : 2.7589
mean_linear_score    : 0.8478
mean_nonlinear_score : 0.7377
mean_probability     : 0.2694
mean_terminal_shear  : 1.2266
mean_peak_shear      : 1.5572
```

# Calibration

## Objective

Score: 0.0027

## Penalties

- `closure_time`: 0.0
- `euploid_prevalence`: 0.0009
- `gap`: 0.0001
- `terminal_shear`: 0.0
- `trisomy21_prevalence`: 0.0

## Euploid Summary

```text
n_embryos            : 96.0
prevalence           : 0.0781
closure_fraction     : 0.3698
mean_closure_time    : 3.3435
mean_gap             : 1.3296
median_gap           : 0.3424
std_gap              : 2.6142
mean_linear_score    : 0.8662
mean_nonlinear_score : 0.7598
mean_probability     : 0.2444
mean_terminal_shear  : 1.298
mean_peak_shear      : 1.5903
```

## Trisomy 21 Summary

```text
n_embryos            : 96.0
prevalence           : 0.2135
closure_fraction     : 0.1146
mean_closure_time    : 3.4555
mean_gap             : 0.6802
median_gap           : 0.417
std_gap              : 1.313
mean_linear_score    : 0.7828
mean_nonlinear_score : 0.6898
mean_probability     : 0.3563
mean_terminal_shear  : 1.256
mean_peak_shear      : 1.3413
```


# Posterior Summary

Acceptance rate: 0.5392
Mean score: 0.4855
MAP score: 0.0027
MAP log posterior: 11.2136
Retained samples: 3804
Burn-in: 250
Max R-hat: 1.0839

## Posterior Means and Intervals

- `alpha_DMP`: 0.6103 [0.4056, 0.8315] (R-hat 1.0839)
- `alpha_EMT`: 0.8145 [0.5942, 1.1199] (R-hat 1.0379)
- `k_A`: 0.6919 [0.5013, 0.9368] (R-hat 1.0472)
- `k_E`: 0.7027 [0.4741, 0.8961] (R-hat 1.0145)
- `k_G`: 0.965 [0.7545, 1.1916] (R-hat 1.0478)
- `threshold`: 0.675 [0.5856, 0.7732] (R-hat 1.006)

## MAP Parameters

- `alpha_DMP`: 0.6133
- `alpha_EMT`: 0.8058
- `k_A`: 0.6641
- `k_E`: 0.742
- `k_G`: 0.9841
- `threshold`: 0.6935


# Validation

Mean score: 0.2172

## Replicate Scores

- replicate 1: 0.0052
- replicate 2: 0.0047
- replicate 3: 0.0
- replicate 4: 0.2738
- replicate 5: 0.0009
- replicate 6: 1.4413
- replicate 7: 0.0113
- replicate 8: 0.0008

## Pass Rates

- `closure_time`: 1.0
- `euploid_prevalence`: 0.625
- `gap`: 0.375
- `terminal_shear`: 1.0
- `trisomy21_prevalence`: 1.0

## Mean Penalties

- `closure_time`: 0.0
- `euploid_prevalence`: 0.0716
- `gap`: 0.0016
- `terminal_shear`: 0.0
- `trisomy21_prevalence`: 0.0


## Sensitivity

Baseline score: 0.0052

### Local Effects

- `alpha_EMT`: -14.5998
- `k_A`: -14.4726
- `alpha_DMP`: -10.8312
- `k_G`: -0.0071
- `k_E`: 0.0
