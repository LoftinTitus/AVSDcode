# AVSD Research Results

## Cohort Summary

```text
n_embryos            : 96.0
prevalence           : 0.0625
closure_fraction     : 0.1146
mean_closure_time    : 3.4566
mean_gap             : 0.6625
std_gap              : 1.3204
mean_linear_score    : 0.7995
mean_nonlinear_score : 0.6857
mean_probability     : 0.2577
mean_terminal_shear  : 1.2356
mean_peak_shear      : 1.3036
```

# Calibration

## Objective

Score: 0.2678

## Penalties

- `closure_time`: 0.0076
- `euploid_prevalence`: 0.0017
- `gap`: 0.1675
- `terminal_shear`: 0.0
- `trisomy21_prevalence`: 0.0

## Euploid Summary

```text
n_embryos            : 96.0
prevalence           : 0.0521
closure_fraction     : 0.1042
mean_closure_time    : 3.4645
mean_gap             : 0.4644
std_gap              : 0.7866
mean_linear_score    : 0.8073
mean_nonlinear_score : 0.6938
mean_probability     : 0.2387
mean_terminal_shear  : 1.3287
mean_peak_shear      : 1.3543
```

## Trisomy 21 Summary

```text
n_embryos            : 96.0
prevalence           : 0.2604
closure_fraction     : 0.0208
mean_closure_time    : 3.4925
mean_gap             : 0.5171
std_gap              : 0.1645
mean_linear_score    : 0.7201
mean_nonlinear_score : 0.6197
mean_probability     : 0.3609
mean_terminal_shear  : 1.1051
mean_peak_shear      : 1.1087
```


# Posterior Summary

Acceptance rate: 0.3197
Mean score: 1.6168
MAP score: 0.2678
MAP log posterior: 10.0192
Retained samples: 420
Burn-in: 32
Max R-hat: 1.5285

## Posterior Means and Intervals

- `alpha_DMP`: 0.6884 [0.5551, 0.8358] (R-hat 1.2562)
- `alpha_EMT`: 0.7019 [0.41, 0.9471] (R-hat 1.5285)
- `k_A`: 0.6126 [0.5221, 0.758] (R-hat 1.2088)
- `k_E`: 0.6985 [0.5098, 0.8602] (R-hat 1.0964)
- `k_G`: 0.8618 [0.7041, 1.0229] (R-hat 1.4103)
- `threshold`: 0.6459 [0.5893, 0.7044] (R-hat 1.0)

## MAP Parameters

- `alpha_DMP`: 0.6721
- `alpha_EMT`: 0.6997
- `k_A`: 0.5508
- `k_E`: 0.7616
- `k_G`: 0.8943
- `threshold`: 0.6367


# Validation

Mean score: 2.3378

## Replicate Scores

- replicate 1: 2.1117
- replicate 2: 1.6428
- replicate 3: 0.6244
- replicate 4: 1.703
- replicate 5: 5.607

## Pass Rates

- `closure_time`: 0.0
- `euploid_prevalence`: 0.2
- `gap`: 0.0
- `terminal_shear`: 1.0
- `trisomy21_prevalence`: 0.8

## Mean Penalties

- `closure_time`: 0.0073
- `euploid_prevalence`: 0.5066
- `gap`: 0.5376
- `terminal_shear`: 0.0
- `trisomy21_prevalence`: 0.0002


## Sensitivity

Baseline score: 2.1117

### Local Effects

- `alpha_DMP`: -44.128
- `k_A`: -32.1597
- `k_G`: 21.9268
- `k_E`: -13.0005
- `alpha_EMT`: -0.5047
