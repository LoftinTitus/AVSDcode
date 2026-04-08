# AVSD Research Results

## Cohort Summary

```text
n_embryos            : 96.0
prevalence           : 0.0521
closure_fraction     : 0.1979
mean_closure_time    : 3.4099
mean_gap             : 0.9943
std_gap              : 2.0978
mean_linear_score    : 0.8356
mean_nonlinear_score : 0.713
mean_probability     : 0.2343
mean_terminal_shear  : 1.261
mean_peak_shear      : 1.4421
```

# Calibration

## Objective

Score: 0.3005

## Penalties

- `closure_time`: 0.0068
- `euploid_prevalence`: 0.0017
- `gap`: 0.19
- `terminal_shear`: 0.0
- `trisomy21_prevalence`: 0.0

## Euploid Summary

```text
n_embryos            : 96.0
prevalence           : 0.0521
closure_fraction     : 0.1354
mean_closure_time    : 3.4517
mean_gap             : 0.5312
std_gap              : 1.0831
mean_linear_score    : 0.831
mean_nonlinear_score : 0.7125
mean_probability     : 0.2331
mean_terminal_shear  : 1.3422
mean_peak_shear      : 1.3895
```

## Trisomy 21 Summary

```text
n_embryos            : 96.0
prevalence           : 0.1667
closure_fraction     : 0.0729
mean_closure_time    : 3.477
mean_gap             : 0.4671
std_gap              : 0.2007
mean_linear_score    : 0.7745
mean_nonlinear_score : 0.6666
mean_probability     : 0.3133
mean_terminal_shear  : 1.2104
mean_peak_shear      : 1.233
```


# Posterior Summary

Acceptance rate: 0.6035
Mean score: 1.5077
MAP score: 0.3005
MAP log posterior: 10.4103
Retained samples: 1028
Burn-in: 128
Max R-hat: 1.5718

## Posterior Means and Intervals

- `alpha_DMP`: 0.6577 [0.4896, 0.8759] (R-hat 1.4633)
- `alpha_EMT`: 0.7097 [0.5718, 0.8412] (R-hat 1.5423)
- `k_A`: 0.6926 [0.5463, 0.9237] (R-hat 1.0179)
- `k_E`: 0.6698 [0.5115, 0.8283] (R-hat 1.5718)
- `k_G`: 0.827 [0.6032, 1.0176] (R-hat 1.1832)
- `threshold`: 0.6376 [0.5627, 0.6999] (R-hat 1.3153)

## MAP Parameters

- `alpha_DMP`: 0.5995
- `alpha_EMT`: 0.7642
- `k_A`: 0.7423
- `k_E`: 0.6202
- `k_G`: 0.8791
- `threshold`: 0.6553


# Validation

Mean score: 3.9997

## Replicate Scores

- replicate 1: 2.6465
- replicate 2: 3.8258
- replicate 3: 6.1296
- replicate 4: 4.4101
- replicate 5: 2.4471
- replicate 6: 4.276
- replicate 7: 2.8883
- replicate 8: 5.3743

## Pass Rates

- `closure_time`: 0.0
- `euploid_prevalence`: 0.625
- `gap`: 0.0
- `terminal_shear`: 1.0
- `trisomy21_prevalence`: 0.0

## Mean Penalties

- `closure_time`: 0.0058
- `euploid_prevalence`: 0.0419
- `gap`: 2.1498
- `terminal_shear`: 0.0
- `trisomy21_prevalence`: 0.1602


## Sensitivity

Baseline score: 2.6465

### Local Effects

- `k_G`: 37.3594
- `alpha_EMT`: 22.6772
- `alpha_DMP`: 21.5612
- `k_A`: 17.9488
- `k_E`: -0.5992
