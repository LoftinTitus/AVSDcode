using AVSDModel

params = load_parameters()
targets = default_calibration_targets()

fit = fit_parameters(
    params = params,
    targets = targets,
    parameter_names = [:alpha_EMT, :alpha_DMP, :k_G, :k_E, :k_A],
    n_candidates = 8,
    n_embryos = 32,
    seed = 2026,
)
posterior = posterior_summary(fit; burn_in = 2)

validation = validate_parameters(
    fit.best_parameters;
    targets = targets,
    n_replicates = 3,
    n_embryos = 32,
    seed = 3030,
)

sensitivity = global_sensitivity(
    [:alpha_EMT, :alpha_DMP, :k_G, :k_E, :k_A];
    params = fit.best_parameters,
    targets = targets,
    n_embryos = 24,
    seed = 4040,
)

println("Best calibration score:")
println(fit.best_evaluation.score)

println("\nPosterior summary:")
println(format_posterior_summary(posterior))

println("\nValidation mean score:")
println(validation.mean_score)

println("\nLocal sensitivity ranking:")
println(sort(collect(sensitivity.local_effects); by = x -> abs(last(x)), rev = true))
