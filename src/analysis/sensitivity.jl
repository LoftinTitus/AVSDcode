function local_sensitivity(
    parameter_names::AbstractVector{Symbol};
    params::ModelParameters = default_parameters(),
    n_embryos::Integer = 64,
    perturbation::Real = 0.05,
    seed::Integer = 42,
    initial_state::PhysicalState = default_initial_state(),
    solver::Symbol = :stochastic_heun,
)
    baseline_results = simulate_population(
        n_embryos;
        params = params,
        initial_state = initial_state,
        seed = seed,
        solver = solver,
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
            solver = solver,
        )

        sensitivities[name] = (prevalence(perturbed_results) - baseline_prevalence) / (params.rates[name] * perturbation)
    end

    return sensitivities
end

struct SensitivityResult{T<:Real}
    baseline_score::T
    baseline_evaluation::CalibrationEvaluation{T}
    local_effects::Dict{Symbol,T}
    global_effects::Dict{Symbol,Dict{Symbol,T}}
    sweeps::Dict{Symbol,Vector{NamedTuple{(:factor, :score, :euploid_prevalence, :trisomy21_prevalence),Tuple{T,T,T,T}}}}
end

function parameter_sweep(
    parameter_name::Symbol;
    params::ModelParameters{T} = default_parameters(),
    targets::Union{Nothing,CalibrationTargets{T}} = nothing,
    factors::Union{Nothing,AbstractVector{T}} = nothing,
    n_embryos::Integer = 64,
    seed::Integer = 42,
    initial_state::Union{Nothing,PhysicalState{T}} = nothing,
    solver::Symbol = :stochastic_heun,
) where {T<:Real}
    targets = isnothing(targets) ? default_calibration_targets(; T = T) : targets
    factors = isnothing(factors) ? T[0.8, 0.9, 1.0, 1.1, 1.2] : factors
    initial_state = isnothing(initial_state) ? default_initial_state() : initial_state

    haskey(params.rates, parameter_name) || error("Unknown rate $(parameter_name)")
    sweep_results = Vector{NamedTuple{(:factor, :score, :euploid_prevalence, :trisomy21_prevalence),Tuple{T,T,T,T}}}()

    for (idx, factor) in enumerate(factors)
        rates = copy(params.rates)
        rates[parameter_name] = max(T(1.0e-6), rates[parameter_name] * factor)
        evaluation = calibration_summary(
            with_rates(params; rates = rates);
            targets = targets,
            n_embryos = n_embryos,
            seed = seed + idx,
            initial_state = initial_state,
            solver = solver,
        )
        push!(sweep_results, (
            factor = factor,
            score = evaluation.score,
            euploid_prevalence = T(evaluation.euploid_summary[:prevalence]),
            trisomy21_prevalence = T(evaluation.trisomy21_summary[:prevalence]),
        ))
    end

    return sweep_results
end

function global_sensitivity(
    parameter_names::AbstractVector{Symbol};
    params::ModelParameters{T} = default_parameters(),
    targets::Union{Nothing,CalibrationTargets{T}} = nothing,
    n_embryos::Integer = 64,
    perturbation::Real = 0.10,
    sweep_factors::Union{Nothing,AbstractVector{T}} = nothing,
    seed::Integer = 42,
    initial_state::Union{Nothing,PhysicalState{T}} = nothing,
    solver::Symbol = :stochastic_heun,
) where {T<:Real}
    targets = isnothing(targets) ? default_calibration_targets(; T = T) : targets
    sweep_factors = isnothing(sweep_factors) ? T[0.8, 0.9, 1.0, 1.1, 1.2] : sweep_factors
    initial_state = isnothing(initial_state) ? default_initial_state() : initial_state

    baseline_evaluation = calibration_summary(
        params;
        targets = targets,
        n_embryos = n_embryos,
        seed = seed,
        initial_state = initial_state,
        solver = solver,
    )

    local_effects = Dict{Symbol,T}()
    global_effects = Dict{Symbol,Dict{Symbol,T}}()
    sweeps = Dict{Symbol,Vector{NamedTuple{(:factor, :score, :euploid_prevalence, :trisomy21_prevalence),Tuple{T,T,T,T}}}}()

    for (idx, parameter_name) in enumerate(parameter_names)
        haskey(params.rates, parameter_name) || continue

        plus_rates = copy(params.rates)
        minus_rates = copy(params.rates)
        plus_rates[parameter_name] *= 1 + perturbation
        minus_rates[parameter_name] = max(T(1.0e-6), minus_rates[parameter_name] * (1 - perturbation))

        plus_eval = calibration_summary(
            with_rates(params; rates = plus_rates);
            targets = targets,
            n_embryos = n_embryos,
            seed = seed + 100 * idx,
            initial_state = initial_state,
            solver = solver,
        )
        minus_eval = calibration_summary(
            with_rates(params; rates = minus_rates);
            targets = targets,
            n_embryos = n_embryos,
            seed = seed + 100 * idx + 1,
            initial_state = initial_state,
            solver = solver,
        )

        scale = max(params.rates[parameter_name] * perturbation, T(1.0e-6))
        local_effects[parameter_name] = (plus_eval.score - minus_eval.score) / (2 * scale)
        global_effects[parameter_name] = Dict(
            :score_span => abs(plus_eval.score - minus_eval.score),
            :euploid_prevalence_shift => abs(T(plus_eval.euploid_summary[:prevalence] - minus_eval.euploid_summary[:prevalence])),
            :trisomy21_prevalence_shift => abs(T(plus_eval.trisomy21_summary[:prevalence] - minus_eval.trisomy21_summary[:prevalence])),
            :closure_time_shift => abs(T(plus_eval.trisomy21_summary[:mean_closure_time] - minus_eval.trisomy21_summary[:mean_closure_time])),
        )
        sweeps[parameter_name] = parameter_sweep(
            parameter_name;
            params = params,
            targets = targets,
            factors = sweep_factors,
            n_embryos = n_embryos,
            seed = seed + 1000 * idx,
            initial_state = initial_state,
            solver = solver,
        )
    end

    return SensitivityResult(baseline_evaluation.score, baseline_evaluation, local_effects, global_effects, sweeps)
end
