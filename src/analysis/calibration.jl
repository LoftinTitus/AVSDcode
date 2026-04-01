struct CalibrationTargets{T<:Real}
    euploid_prevalence::Tuple{T,T}
    trisomy21_prevalence::Tuple{T,T}
    closure_time::Tuple{T,T}
    shear_range::Tuple{T,T}
    max_gap::T
    weights::Dict{Symbol,T}
end

struct CalibrationEvaluation{T<:Real}
    score::T
    penalties::Dict{Symbol,T}
    euploid_summary::Dict{Symbol,Float64}
    trisomy21_summary::Dict{Symbol,Float64}
end

struct CalibrationResult{T<:Real}
    best_parameters::ModelParameters{T}
    best_evaluation::CalibrationEvaluation{T}
    history::Vector{Dict{Symbol,Any}}
end

struct ValidationResult{T<:Real}
    mean_score::T
    replicate_scores::Vector{T}
    pass_rates::Dict{Symbol,T}
    mean_penalties::Dict{Symbol,T}
    euploid_summary::Dict{Symbol,Float64}
    trisomy21_summary::Dict{Symbol,Float64}
end

function default_calibration_targets(; T::Type{<:Real} = Float64)
    weights = Dict{Symbol,T}(
        :euploid_prevalence => T(3.0),
        :trisomy21_prevalence => T(4.0),
        :closure_time => T(1.5),
        :terminal_shear => T(1.0),
        :gap => T(1.5),
    )

    return CalibrationTargets(
        (T(0.00), T(0.05)),
        (T(0.15), T(0.33)),
        (T(1.5), T(3.2)),
        (T(0.5), T(4.0)),
        T(0.35),
        weights,
    )
end

function _range_penalty(value::Real, bounds::Tuple{<:Real,<:Real})
    lower, upper = bounds

    if value < lower
        scale = max(abs(lower), 1.0e-6)
        return ((lower - value) / scale)^2
    elseif value > upper
        scale = max(abs(upper), 1.0e-6)
        return ((value - upper) / scale)^2
    end

    return 0.0
end

function _upper_penalty(value::Real, upper::Real)
    value <= upper && return 0.0
    scale = max(abs(upper), 1.0e-6)
    return ((value - upper) / scale)^2
end

function _mean_summary_dict(dicts::AbstractVector{<:Dict{Symbol,Float64}})
    isempty(dicts) && return Dict{Symbol,Float64}()
    common_keys = collect(keys(first(dicts)))
    return Dict(
        key => mean(dict[key] for dict in dicts)
        for key in common_keys
    )
end

function evaluate_calibration(
    euploid_results::AbstractVector{<:EmbryoResult},
    trisomy21_results::AbstractVector{<:EmbryoResult},
    targets::CalibrationTargets{T} = default_calibration_targets(; T = Float64),
) where {T<:Real}
    euploid_summary = summary_metrics(euploid_results)
    trisomy21_summary = summary_metrics(trisomy21_results)

    penalties = Dict{Symbol,T}(
        :euploid_prevalence => T(_range_penalty(euploid_summary[:prevalence], targets.euploid_prevalence)),
        :trisomy21_prevalence => T(_range_penalty(trisomy21_summary[:prevalence], targets.trisomy21_prevalence)),
        :closure_time => T(0.5 * (
            _range_penalty(euploid_summary[:mean_closure_time], targets.closure_time) +
            _range_penalty(trisomy21_summary[:mean_closure_time], targets.closure_time)
        )),
        :terminal_shear => T(0.5 * (
            _range_penalty(euploid_summary[:mean_terminal_shear], targets.shear_range) +
            _range_penalty(trisomy21_summary[:mean_terminal_shear], targets.shear_range)
        )),
        :gap => T(0.5 * (
            _upper_penalty(euploid_summary[:mean_gap], targets.max_gap) +
            _upper_penalty(trisomy21_summary[:mean_gap], targets.max_gap)
        )),
    )

    score = zero(T)
    for (penalty_name, penalty_value) in penalties
        score += get(targets.weights, penalty_name, one(T)) * penalty_value
    end

    return CalibrationEvaluation(score, penalties, euploid_summary, trisomy21_summary)
end

function calibration_summary(
    params::ModelParameters{T};
    targets::Union{Nothing,CalibrationTargets{T}} = nothing,
    n_embryos::Integer = 96,
    seed::Integer = 42,
    initial_state::Union{Nothing,PhysicalState{T}} = nothing,
    solver::Symbol = :stochastic_heun,
) where {T<:Real}
    targets = isnothing(targets) ? default_calibration_targets(; T = T) : targets
    initial_state = isnothing(initial_state) ? default_initial_state() : initial_state

    euploid_results = simulate_population(
        n_embryos;
        params = params,
        initial_state = initial_state,
        seed = seed,
        trisomy_fraction = 0.0,
        solver = solver,
    )
    trisomy21_results = simulate_population(
        n_embryos;
        params = params,
        initial_state = initial_state,
        seed = seed + 10_000,
        trisomy_fraction = 1.0,
        solver = solver,
    )

    return evaluate_calibration(euploid_results, trisomy21_results, targets)
end

function _prior_group(priors::Dict{String,Any}, parameter_name::Symbol)
    for group_name in ("rates", "noise")
        group = get(priors, group_name, Dict{String,Any}())
        haskey(group, string(parameter_name, "_mean")) && return group
    end

    return Dict{String,Any}()
end

function sample_parameter_candidate(
    rng::AbstractRNG,
    params::ModelParameters{T},
    parameter_names::AbstractVector{Symbol},
    priors::Dict{String,Any} = load_priors(),
) where {T<:Real}
    candidate_rates = copy(params.rates)

    for name in parameter_names
        haskey(candidate_rates, name) || continue
        group = _prior_group(priors, name)
        mean_key = string(name, "_mean")
        sd_key = string(name, "_sd")
        baseline = candidate_rates[name]
        mean_value = T(get(group, mean_key, baseline))
        sd_value = T(get(group, sd_key, max(abs(baseline) * 0.15, 1.0e-3)))
        sampled_value = mean_value + randn(rng) * sd_value
        candidate_rates[name] = max(T(1.0e-6), sampled_value)
    end

    return with_rates(params; rates = candidate_rates)
end

function fit_parameters(
    ;
    params::ModelParameters{T} = default_parameters(),
    targets::Union{Nothing,CalibrationTargets{T}} = nothing,
    parameter_names::Union{Nothing,AbstractVector{Symbol}} = nothing,
    priors::Dict{String,Any} = load_priors(),
    n_candidates::Integer = 24,
    n_embryos::Integer = 64,
    seed::Integer = 42,
    initial_state::Union{Nothing,PhysicalState{T}} = nothing,
    solver::Symbol = :stochastic_heun,
) where {T<:Real}
    targets = isnothing(targets) ? default_calibration_targets(; T = T) : targets
    parameter_names = isnothing(parameter_names) ? collect(keys(params.rates)) : parameter_names
    initial_state = isnothing(initial_state) ? default_initial_state() : initial_state

    rng = MersenneTwister(seed)
    history = Vector{Dict{Symbol,Any}}()

    best_params = params
    best_eval = calibration_summary(
        params;
        targets = targets,
        n_embryos = n_embryos,
        seed = seed,
        initial_state = initial_state,
        solver = solver,
    )

    push!(history, Dict(
        :iteration => 0,
        :score => best_eval.score,
        :rates => copy(best_params.rates),
        :penalties => copy(best_eval.penalties),
    ))

    for iteration in 1:n_candidates
        candidate_params = sample_parameter_candidate(rng, params, parameter_names, priors)
        evaluation = calibration_summary(
            candidate_params;
            targets = targets,
            n_embryos = n_embryos,
            seed = seed + iteration,
            initial_state = initial_state,
            solver = solver,
        )

        push!(history, Dict(
            :iteration => iteration,
            :score => evaluation.score,
            :rates => copy(candidate_params.rates),
            :penalties => copy(evaluation.penalties),
        ))

        if evaluation.score < best_eval.score
            best_params = candidate_params
            best_eval = evaluation
        end
    end

    return CalibrationResult(best_params, best_eval, history)
end

function validate_parameters(
    params::ModelParameters{T};
    targets::Union{Nothing,CalibrationTargets{T}} = nothing,
    n_replicates::Integer = 5,
    n_embryos::Integer = 96,
    seed::Integer = 42,
    initial_state::Union{Nothing,PhysicalState{T}} = nothing,
    solver::Symbol = :stochastic_heun,
) where {T<:Real}
    targets = isnothing(targets) ? default_calibration_targets(; T = T) : targets
    initial_state = isnothing(initial_state) ? default_initial_state() : initial_state

    evaluations = Vector{CalibrationEvaluation{T}}(undef, n_replicates)

    for replicate in 1:n_replicates
        evaluations[replicate] = calibration_summary(
            params;
            targets = targets,
            n_embryos = n_embryos,
            seed = seed + 1000 * replicate,
            initial_state = initial_state,
            solver = solver,
        )
    end

    pass_rates = Dict{Symbol,T}()
    for penalty_name in keys(evaluations[1].penalties)
        pass_rates[penalty_name] = T(mean(evaluation -> iszero(evaluation.penalties[penalty_name]) ? 1.0 : 0.0, evaluations))
    end

    return ValidationResult(
        T(mean(evaluation.score for evaluation in evaluations)),
        T[evaluation.score for evaluation in evaluations],
        pass_rates,
        Dict(
            key => T(mean(evaluation.penalties[key] for evaluation in evaluations))
            for key in keys(evaluations[1].penalties)
        ),
        _mean_summary_dict([evaluation.euploid_summary for evaluation in evaluations]),
        _mean_summary_dict([evaluation.trisomy21_summary for evaluation in evaluations]),
    )
end
