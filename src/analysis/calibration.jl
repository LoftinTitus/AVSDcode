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

struct MultiChainCalibration{T<:Real}
    chains::Vector{CalibrationResult{T}}
    best_chain_index::Int
    best_parameters::ModelParameters{T}
    best_evaluation::CalibrationEvaluation{T}
end

struct PosteriorSummary{T<:Real}
    acceptance_rate::T
    mean_score::T
    map_score::T
    map_log_posterior::T
    posterior_mean::Dict{Symbol,T}
    posterior_interval::Dict{Symbol,Tuple{T,T}}
    map_parameters::Dict{Symbol,T}
    rhat::Dict{Symbol,T}
    n_samples::Int
    burn_in::Int
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
    sample_strategy::Symbol = :stratified,
    calibration_panel::Union{Nothing,CalibrationPanel{T}} = nothing,
) where {T<:Real}
    targets = isnothing(targets) ? default_calibration_targets(; T = T) : targets
    initial_state = isnothing(initial_state) ? default_initial_state() : initial_state

    panel = isnothing(calibration_panel) ? build_calibration_panel(
        n_embryos;
        params = params,
        seed = seed,
        sample_strategy = sample_strategy,
    ) : calibration_panel

    euploid_results = simulate_population(
        panel.euploid;
        params = params,
        initial_state = initial_state,
        solver = solver,
    )
    trisomy21_results = simulate_population(
        panel.trisomy21;
        params = params,
        initial_state = initial_state,
        solver = solver,
    )

    return evaluate_calibration(euploid_results, trisomy21_results, targets)
end

function _prior_group(priors::Dict{String,Any}, parameter_name::Symbol)
    parameter_key = string(parameter_name)

    for group_name in ("rates", "noise", "phenotype")
        group = get(priors, group_name, Dict{String,Any}())
        haskey(group, parameter_key * "_mean") && return group_name, group
    end

    return "", Dict{String,Any}()
end

function _parameter_baseline(params::ModelParameters, parameter_name::Symbol)
    if haskey(params.rates, parameter_name)
        return float(params.rates[parameter_name])
    elseif parameter_name == :threshold
        return float(params.phenotype.threshold)
    elseif parameter_name == :slope
        return float(params.phenotype.slope)
    else
        error("Unknown calibration parameter $(parameter_name)")
    end
end

function _prior_statistics(priors::Dict{String,Any}, params::ModelParameters, parameter_name::Symbol)
    _, group = _prior_group(priors, parameter_name)
    baseline = _parameter_baseline(params, parameter_name)
    mean_value = float(get(group, string(parameter_name, "_mean"), baseline))
    sd_value = float(get(group, string(parameter_name, "_sd"), max(abs(baseline) * 0.15, 1.0e-3)))
    return mean_value, max(sd_value, 1.0e-6)
end

function _parameter_value_dict(params::ModelParameters, parameter_names::AbstractVector{Symbol})
    values = Dict{Symbol,Float64}()
    for name in parameter_names
        values[name] = _parameter_baseline(params, name)
    end
    return values
end

function _with_parameter_values(
    params::ModelParameters{T},
    values::Dict{Symbol,<:Real},
) where {T<:Real}
    rates = copy(params.rates)
    threshold = params.phenotype.threshold
    slope = params.phenotype.slope

    for (name, raw_value) in values
        value = T(max(float(raw_value), 1.0e-6))
        if haskey(rates, name)
            rates[name] = value
        elseif name == :threshold
            threshold = clamp(value, T(1.0e-4), T(1 - 1.0e-4))
        elseif name == :slope
            slope = max(value, T(1.0e-3))
        else
            error("Unknown calibration parameter $(name)")
        end
    end

    phenotype = PhenotypeParameters(
        params.phenotype.w_C,
        params.phenotype.w_D,
        params.phenotype.w_E,
        params.phenotype.w_A,
        threshold,
        slope,
    )

    return ModelParameters(
        params.T,
        params.dt,
        params.windows,
        params.hemo,
        params.activation,
        phenotype,
        rates,
        params.genetics,
    )
end

function _log_prior_density(
    values::Dict{Symbol,<:Real},
    params::ModelParameters,
    parameter_names::AbstractVector{Symbol},
    priors::Dict{String,Any},
)
    log_density = 0.0

    for name in parameter_names
        value = float(values[name])
        value <= 0.0 && return -Inf
        mean_value, sd_value = _prior_statistics(priors, params, name)
        z = (value - mean_value) / sd_value
        log_density += -0.5 * z^2 - log(sd_value)
    end

    return log_density
end

function _log_posterior(
    evaluation::CalibrationEvaluation,
    log_prior::Real;
    penalty_scale::Real = 1.0,
)
    return float(log_prior) - float(penalty_scale) * float(evaluation.score)
end

function sample_parameter_candidate(
    rng::AbstractRNG,
    params::ModelParameters{T},
    parameter_names::AbstractVector{Symbol},
    priors::Dict{String,Any} = load_priors(),
) where {T<:Real}
    candidate_values = Dict{Symbol,Float64}()

    for name in parameter_names
        mean_value, sd_value = _prior_statistics(priors, params, name)
        candidate_values[name] = max(1.0e-6, mean_value + randn(rng) * sd_value)
    end

    return _with_parameter_values(params, candidate_values)
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
    sample_strategy::Symbol = :stratified,
    proposal_scale::Real = 0.20,
    prior_mix::Real = 0.15,
    refinement_rounds::Integer = 0,
    refinement_candidates::Integer = max(cld(Int(n_candidates), 2), 1),
    refinement_scale_factor::Real = 0.5,
    penalty_scale::Real = 1.0,
) where {T<:Real}
    targets = isnothing(targets) ? default_calibration_targets(; T = T) : targets
    parameter_names = isnothing(parameter_names) ? collect(keys(params.rates)) : collect(parameter_names)
    initial_state = isnothing(initial_state) ? default_initial_state() : initial_state
    calibration_panel = build_calibration_panel(
        n_embryos;
        params = params,
        seed = seed,
        sample_strategy = sample_strategy,
    )

    rng = MersenneTwister(seed)
    history = Vector{Dict{Symbol,Any}}()

    current_values = _parameter_value_dict(params, parameter_names)
    current_params = _with_parameter_values(params, current_values)
    current_eval = calibration_summary(
        current_params;
        targets = targets,
        n_embryos = n_embryos,
        seed = seed,
        initial_state = initial_state,
        solver = solver,
        sample_strategy = sample_strategy,
        calibration_panel = calibration_panel,
    )
    current_log_prior = _log_prior_density(current_values, params, parameter_names, priors)
    current_log_posterior = _log_posterior(current_eval, current_log_prior; penalty_scale = penalty_scale)

    best_params = current_params
    best_eval = current_eval
    best_values = copy(current_values)
    best_log_prior = current_log_prior
    best_log_posterior = current_log_posterior

    push!(history, Dict(
        :iteration => 0,
        :score => current_eval.score,
        :log_prior => current_log_prior,
        :log_posterior => current_log_posterior,
        :accepted => true,
        :parameter_values => copy(current_values),
        :rates => copy(current_params.rates),
        :penalties => copy(current_eval.penalties),
    ))

    function run_walk!(
        iteration_offset::Int,
        stage_candidates::Int,
        stage_proposal_scale::Real,
        stage_prior_mix::Real,
    )
        for stage_iteration in 1:stage_candidates
            iteration = iteration_offset + stage_iteration
            proposed_values = copy(current_values)

            for name in parameter_names
                current_value = max(proposed_values[name], 1.0e-6)
                mean_value, sd_value = _prior_statistics(priors, params, name)
                log_step = randn(rng) * stage_proposal_scale
                proposed_value = current_value * exp(log_step)

                if stage_prior_mix > 0
                    proposed_value = (1 - stage_prior_mix) * proposed_value + stage_prior_mix * max(1.0e-6, mean_value + randn(rng) * sd_value)
                end

                proposed_values[name] = max(1.0e-6, proposed_value)
            end

            proposed_params = _with_parameter_values(params, proposed_values)
            proposed_eval = calibration_summary(
                proposed_params;
                targets = targets,
                n_embryos = n_embryos,
                seed = seed + iteration,
                initial_state = initial_state,
                solver = solver,
                sample_strategy = sample_strategy,
                calibration_panel = calibration_panel,
            )
            proposed_log_prior = _log_prior_density(proposed_values, params, parameter_names, priors)
            proposed_log_posterior = _log_posterior(proposed_eval, proposed_log_prior; penalty_scale = penalty_scale)

            log_acceptance = proposed_log_posterior - current_log_posterior
            accepted = log(rand(rng)) < min(0.0, log_acceptance)

            if accepted
                current_values = proposed_values
                current_params = proposed_params
                current_eval = proposed_eval
                current_log_prior = proposed_log_prior
                current_log_posterior = proposed_log_posterior
            end

            if proposed_log_posterior > best_log_posterior
                best_values = copy(proposed_values)
                best_params = proposed_params
                best_eval = proposed_eval
                best_log_prior = proposed_log_prior
                best_log_posterior = proposed_log_posterior
            end

            push!(history, Dict(
                :iteration => iteration,
                :score => current_eval.score,
                :log_prior => current_log_prior,
                :log_posterior => current_log_posterior,
                :accepted => accepted,
                :parameter_values => copy(current_values),
                :rates => copy(current_params.rates),
                :penalties => copy(current_eval.penalties),
            ))
        end

        return iteration_offset + stage_candidates
    end

    last_iteration = run_walk!(0, Int(n_candidates), proposal_scale, prior_mix)

    refined_candidates = max(Int(refinement_candidates), 1)
    for refinement_round in 1:Int(refinement_rounds)
        current_values = copy(best_values)
        current_params = best_params
        current_eval = best_eval
        current_log_prior = best_log_prior
        current_log_posterior = best_log_posterior

        stage_scale = max(float(proposal_scale) * refinement_scale_factor^refinement_round, 1.0e-3)
        stage_prior_mix = max(float(prior_mix) * refinement_scale_factor^refinement_round, 0.0)
        last_iteration = run_walk!(last_iteration, refined_candidates, stage_scale, stage_prior_mix)
    end

    return CalibrationResult(best_params, best_eval, history)
end

function posterior_summary(
    calibration::CalibrationResult{T};
    burn_in::Integer = 0,
    interval_level::Real = 0.95,
) where {T<:Real}
    isempty(calibration.history) && error("Calibration history is empty")

    capped_burn_in = clamp(Int(burn_in), 0, max(length(calibration.history) - 1, 0))
    samples = calibration.history[(capped_burn_in + 1):end]
    isempty(samples) && error("No posterior samples remain after burn-in")

    proposal_samples = calibration.history[min(length(calibration.history), capped_burn_in + 2):end]
    acceptance_rate = isempty(proposal_samples) ? 1.0 : mean(item -> item[:accepted] ? 1.0 : 0.0, proposal_samples)
    mean_score = mean(item -> float(item[:score]), samples)
    map_item = samples[argmax([float(item[:log_posterior]) for item in samples])]
    parameter_names = sort!(collect(keys(samples[1][:parameter_values])); by = string)
    posterior_mean = Dict(
        name => T(mean(item -> float(item[:parameter_values][name]), samples))
        for name in parameter_names
    )
    posterior_interval = Dict(
        name => _credible_interval(
            [float(item[:parameter_values][name]) for item in samples];
            level = interval_level,
            T = T,
        )
        for name in parameter_names
    )
    map_parameters = Dict(
        name => T(float(map_item[:parameter_values][name]))
        for name in parameter_names
    )

    return PosteriorSummary(
        T(acceptance_rate),
        T(mean_score),
        T(float(map_item[:score])),
        T(float(map_item[:log_posterior])),
        posterior_mean,
        posterior_interval,
        map_parameters,
        Dict{Symbol,T}(),
        length(samples),
        capped_burn_in,
    )
end

function _history_log_posterior_peak(calibration::CalibrationResult)
    return maximum(float(item[:log_posterior]) for item in calibration.history)
end

function _quantile(values::AbstractVector{<:Real}, p::Real)
    isempty(values) && error("Cannot compute quantile of an empty collection")
    0.0 <= p <= 1.0 || error("Quantile probability must lie between 0 and 1")

    sorted = sort(float.(collect(values)))
    length(sorted) == 1 && return sorted[1]

    rank = 1 + (length(sorted) - 1) * p
    lower = floor(Int, rank)
    upper = ceil(Int, rank)
    lower == upper && return sorted[lower]

    fraction = rank - lower
    return (1 - fraction) * sorted[lower] + fraction * sorted[upper]
end

function _credible_interval(
    values::AbstractVector{<:Real};
    level::Real = 0.95,
    T::Type{<:Real} = Float64,
)
    0.0 < level < 1.0 || error("Credible interval level must lie between 0 and 1")
    alpha = (1 - level) / 2
    return T(_quantile(values, alpha)), T(_quantile(values, 1 - alpha))
end

function _rhat(chains::AbstractVector{<:AbstractVector{<:Real}})
    length(chains) < 2 && return 1.0
    n = minimum(length(chain) for chain in chains)
    n < 2 && return 1.0

    trimmed = [float.(collect(chain[(end - n + 1):end])) for chain in chains]
    chain_means = [mean(chain) for chain in trimmed]
    chain_vars = [var(chain; corrected = true) for chain in trimmed]
    W = mean(chain_vars)
    W <= 1.0e-12 && return 1.0

    B = n * var(chain_means; corrected = true)
    var_hat = ((n - 1) / n) * W + B / n
    return sqrt(max(var_hat / W, 1.0))
end

function run_calibration_chains(
    ;
    n_chains::Integer = 4,
    chain_seed_stride::Integer = 10_000,
    seed::Integer = 42,
    kwargs...,
)
    n_chains >= 1 || error("n_chains must be at least 1")

    first_chain = fit_parameters(; seed = seed, kwargs...)
    chains = Vector{typeof(first_chain)}(undef, n_chains)
    chains[1] = first_chain

    for chain_idx in 2:n_chains
        chains[chain_idx] = fit_parameters(
            ;
            seed = seed + chain_seed_stride * (chain_idx - 1),
            kwargs...,
        )
    end

    best_scores = [_history_log_posterior_peak(chain) for chain in chains]
    best_chain_index = argmax(best_scores)
    best_chain = chains[best_chain_index]

    return MultiChainCalibration(chains, best_chain_index, best_chain.best_parameters, best_chain.best_evaluation)
end

function posterior_summary(
    calibration::MultiChainCalibration{T};
    burn_in::Integer = 0,
    interval_level::Real = 0.95,
) where {T<:Real}
    isempty(calibration.chains) && error("Multi-chain calibration has no chains")

    chain_samples = Vector{Vector{Dict{Symbol,Any}}}(undef, length(calibration.chains))
    acceptance_rates = Float64[]
    capped_burn_in = 0

    for (idx, chain) in enumerate(calibration.chains)
        isempty(chain.history) && error("Calibration chain $(idx) has empty history")
        capped_burn_in = clamp(Int(burn_in), 0, max(length(chain.history) - 1, 0))
        chain_samples[idx] = chain.history[(capped_burn_in + 1):end]
        proposal_samples = chain.history[min(length(chain.history), capped_burn_in + 2):end]
        push!(acceptance_rates, isempty(proposal_samples) ? 1.0 : mean(item -> item[:accepted] ? 1.0 : 0.0, proposal_samples))
    end

    samples = reduce(vcat, chain_samples)
    isempty(samples) && error("No posterior samples remain after burn-in")

    map_item = samples[argmax([float(item[:log_posterior]) for item in samples])]
    parameter_names = sort!(collect(keys(samples[1][:parameter_values])); by = string)
    posterior_mean = Dict(
        name => T(mean(item -> float(item[:parameter_values][name]), samples))
        for name in parameter_names
    )
    posterior_interval = Dict(
        name => _credible_interval(
            [float(item[:parameter_values][name]) for item in samples];
            level = interval_level,
            T = T,
        )
        for name in parameter_names
    )
    map_parameters = Dict(
        name => T(float(map_item[:parameter_values][name]))
        for name in parameter_names
    )

    rhat = Dict{Symbol,T}()
    for name in parameter_names
        rhat[name] = T(_rhat([[float(item[:parameter_values][name]) for item in chain] for chain in chain_samples]))
    end
    rhat[:score] = T(_rhat([[float(item[:score]) for item in chain] for chain in chain_samples]))

    return PosteriorSummary(
        T(mean(acceptance_rates)),
        T(mean(item -> float(item[:score]), samples)),
        T(float(map_item[:score])),
        T(float(map_item[:log_posterior])),
        posterior_mean,
        posterior_interval,
        map_parameters,
        rhat,
        length(samples),
        capped_burn_in,
    )
end

function validate_parameters(
    params::ModelParameters{T};
    targets::Union{Nothing,CalibrationTargets{T}} = nothing,
    n_replicates::Integer = 5,
    n_embryos::Integer = 96,
    seed::Integer = 42,
    initial_state::Union{Nothing,PhysicalState{T}} = nothing,
    solver::Symbol = :stochastic_heun,
    sample_strategy::Symbol = :stratified,
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
            sample_strategy = sample_strategy,
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
