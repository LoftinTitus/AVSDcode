# AVSD Research Results

## Cohort Summary

```text
n_embryos            : 96.0
prevalence           : 0.0625
closure_fraction     : 0.1771
mean_closure_time    : 3.4231
mean_gap             : 0.7695
std_gap              : 1.6762
mean_linear_score    : 0.8424
mean_nonlinear_score : 0.7409
mean_probability     : 0.2336
mean_terminal_shear  : 1.2855
mean_peak_shear      : 1.3977
```

# Calibration

## Objective

Score: 0.4479

## Penalties

- `closure_time`: 0.0074
- `euploid_prevalence`: 0.0
- `gap`: 0.2171
- `terminal_shear`: 0.0
- `trisomy21_prevalence`: 0.0278

## Euploid Summary

```text
n_embryos            : 96.0
prevalence           : 0.0312
closure_fraction     : 0.125
mean_closure_time    : 3.4606
mean_gap             : 0.4787
std_gap              : 0.8073
mean_linear_score    : 0.8428
mean_nonlinear_score : 0.7462
mean_probability     : 0.2271
mean_terminal_shear  : 1.3199
mean_peak_shear      : 1.3679
```

## Trisomy 21 Summary

```text
n_embryos            : 96.0
prevalence           : 0.125
closure_fraction     : 0.0208
mean_closure_time    : 3.4894
mean_gap             : 0.5414
std_gap              : 0.7689
mean_linear_score    : 0.787
mean_nonlinear_score : 0.6953
mean_probability     : 0.2901
mean_terminal_shear  : 1.1943
mean_peak_shear      : 1.2211
```


# Posterior Summary

Acceptance rate: 0.6453
Mean score: 1.5428
MAP score: 0.4479
MAP log posterior: 10.262
Retained samples: 644
Burn-in: 64
Max R-hat: 1.9122

## Posterior Means and Intervals

- `alpha_DMP`: 0.6475 [0.4861, 0.8439] (R-hat 1.0497)
- `alpha_EMT`: 0.7443 [0.5613, 0.9385] (R-hat 1.9122)
- `k_A`: 0.6423 [0.4357, 0.8208] (R-hat 1.189)
- `k_E`: 0.6734 [0.4361, 0.845] (R-hat 1.5145)
- `k_G`: 0.8509 [0.7087, 1.0106] (R-hat 1.1682)
- `threshold`: 0.6317 [0.5393, 0.6919] (R-hat 1.1388)

## MAP Parameters

- `alpha_DMP`: 0.5366
- `alpha_EMT`: 0.8983
- `k_A`: 0.6612
- `k_E`: 0.7855
- `k_G`: 0.8457
- `threshold`: 0.6605


# Validation

Mean score: 2.3525

## Replicate Scores

- replicate 1: 1.9429
- replicate 2: 2.8343
- replicate 3: 3.7509
- replicate 4: 1.768
- replicate 5: 1.4663

## Pass Rates

- `closure_time`: 0.0
- `euploid_prevalence`: 0.6
- `gap`: 0.0
- `terminal_shear`: 1.0
- `trisomy21_prevalence`: 0.2

## Mean Penalties

- `closure_time`: 0.0063
- `euploid_prevalence`: 0.0128
- `gap`: 1.3058
- `terminal_shear`: 0.0
- `trisomy21_prevalence`: 0.0865


## Sensitivity

Baseline score: 1.9429

### Local Effects

- `k_G`: 23.4244
- `k_A`: 16.0856
- `alpha_DMP`: 12.6988
- `alpha_EMT`: 9.4181
- `k_E`: -1.5399
