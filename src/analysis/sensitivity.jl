function local_sensitivity(
    parameter_names::AbstractVector{Symbol};
    params::ModelParameters = default_parameters(),
    n_embryos::Integer = 64,
    perturbation::Real = 0.05,
    seed::Integer = 42,
    initial_state::PhysicalState = default_initial_state(),
)
    baseline_results = simulate_population(
        n_embryos;
        params = params,
        initial_state = initial_state,
        seed = seed,
    )
    baseline_prevalence = prevalence(baseline_results)

    sensitivities = Dict{Symbol,Float64}()

    for (idx, name) in enumerate(parameter_names)
        haskey(params.rates, name) || continue

        perturbed_rates = copy(params.rates)
        perturbed_rates[name] *= 1 + perturbation
        perturbed_params = with_rates(params; rates = perturbed_rates)

        perturbed_results = simulate_population(
            n_embryos;
            params = perturbed_params,
            initial_state = initial_state,
            seed = seed + idx,
        )

        sensitivities[name] = (prevalence(perturbed_results) - baseline_prevalence) / (params.rates[name] * perturbation)
    end

    return sensitivities
end
