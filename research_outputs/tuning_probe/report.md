# AVSD Research Results

## Cohort Summary

```text
n_embryos            : 96.0
prevalence           : 0.0521
closure_fraction     : 0.1146
mean_closure_time    : 3.4578
mean_gap             : 0.666
std_gap              : 1.3244
mean_linear_score    : 0.8488
mean_nonlinear_score : 0.7147
mean_probability     : 0.2308
mean_terminal_shear  : 1.2289
mean_peak_shear      : 1.2978
```

# Calibration

## Objective

Score: 0.4524

## Penalties

- `closure_time`: 0.008
- `euploid_prevalence`: 0.0625
- `gap`: 0.1686
- `terminal_shear`: 0.0
- `trisomy21_prevalence`: 0.0

## Euploid Summary

```text
n_embryos            : 96.0
prevalence           : 0.0625
closure_fraction     : 0.0625
mean_closure_time    : 3.4768
mean_gap             : 0.4973
std_gap              : 0.7752
mean_linear_score    : 0.8567
mean_nonlinear_score : 0.7139
mean_probability     : 0.2192
mean_terminal_shear  : 1.2504
mean_peak_shear      : 1.276
```

## Trisomy 21 Summary

```text
n_embryos            : 96.0
prevalence           : 0.1667
closure_fraction     : 0.0312
mean_closure_time    : 3.4951
mean_gap             : 0.4901
std_gap              : 0.1759
mean_linear_score    : 0.7762
mean_nonlinear_score : 0.6589
mean_probability     : 0.3256
mean_terminal_shear  : 1.1792
mean_peak_shear      : 1.183
```


# Posterior Summary

Acceptance rate: 0.2321
Mean score: 1.4498
MAP score: 0.4524
MAP log posterior: 10.0849
Retained samples: 228
Burn-in: 24
Max R-hat: 1.6314

## Posterior Means and Intervals

- `alpha_DMP`: 0.6788 [0.5283, 0.827] (R-hat 1.6314)
- `alpha_EMT`: 0.6947 [0.5444, 0.8428] (R-hat 1.2)
- `k_A`: 0.6364 [0.4653, 0.7557] (R-hat 1.3139)
- `k_E`: 0.7068 [0.452, 0.8664] (R-hat 1.2324)
- `k_G`: 0.8651 [0.7336, 0.9926] (R-hat 1.2262)
- `threshold`: 0.642 [0.5852, 0.7002] (R-hat 1.1329)

## MAP Parameters

- `alpha_DMP`: 0.6636
- `alpha_EMT`: 0.7576
- `k_A`: 0.6823
- `k_E`: 0.6317
- `k_G`: 0.783
- `threshold`: 0.6626


# Validation

Mean score: 1.208

## Replicate Scores

- replicate 1: 1.0237
- replicate 2: 1.2512
- replicate 3: 0.7485
- replicate 4: 0.7957
- replicate 5: 2.221

## Pass Rates

- `closure_time`: 0.0
- `euploid_prevalence`: 0.2
- `gap`: 0.0
- `terminal_shear`: 1.0
- `trisomy21_prevalence`: 0.2

## Mean Penalties

- `closure_time`: 0.0074
- `euploid_prevalence`: 0.1142
- `gap`: 0.5148
- `terminal_shear`: 0.0
- `trisomy21_prevalence`: 0.0205


## Sensitivity

Baseline score: 1.0237

### Local Effects

- `k_G`: 26.475
- `k_E`: -11.7716
- `alpha_EMT`: 10.921
- `alpha_DMP`: -10.6582
- `k_A`: 7.7499
