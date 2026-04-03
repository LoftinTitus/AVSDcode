function _format_float(value::Real; digits::Integer = 4)
    return string(round(float(value); digits = digits))
end

function _summary_key_order(summary::Dict{Symbol,Float64})
    preferred = [
        :n_embryos,
        :prevalence,
        :closure_fraction,
        :mean_closure_time,
        :mean_gap,
        :std_gap,
        :mean_linear_score,
        :mean_nonlinear_score,
        :mean_probability,
        :mean_terminal_shear,
        :mean_peak_shear,
    ]
    present = [key for key in preferred if haskey(summary, key)]
    extras = sort([key for key in keys(summary) if key ∉ present], by = string)
    return vcat(present, extras)
end

function format_summary_table(summary::Dict{Symbol,Float64}; digits::Integer = 4)
    isempty(summary) && return "No summary metrics available.\n"

    ordered_keys = _summary_key_order(summary)
    label_width = maximum(length(string(key)) for key in ordered_keys)
    lines = [
        rpad(string(key), label_width) * " : " * _format_float(summary[key]; digits = digits)
        for key in ordered_keys
    ]
    return join(lines, "\n") * "\n"
end

function format_calibration_report(
    evaluation::CalibrationEvaluation;
    digits::Integer = 4,
    title::AbstractString = "Calibration Report",
)
    sections = String["# " * title, "", "## Objective", "", "Score: " * _format_float(evaluation.score; digits = digits), ""]

    push!(sections, "## Penalties")
    push!(sections, "")
    for key in sort(collect(keys(evaluation.penalties)); by = string)
        push!(sections, "- `" * string(key) * "`: " * _format_float(evaluation.penalties[key]; digits = digits))
    end
    push!(sections, "")

    push!(sections, "## Euploid Summary")
    push!(sections, "")
    push!(sections, "```text")
    push!(sections, chomp(format_summary_table(evaluation.euploid_summary; digits = digits)))
    push!(sections, "```")
    push!(sections, "")

    push!(sections, "## Trisomy 21 Summary")
    push!(sections, "")
    push!(sections, "```text")
    push!(sections, chomp(format_summary_table(evaluation.trisomy21_summary; digits = digits)))
    push!(sections, "```")
    push!(sections, "")

    return join(sections, "\n")
end

function format_validation_report(
    validation::ValidationResult;
    digits::Integer = 4,
    title::AbstractString = "Validation Report",
)
    lines = String[
        "# " * title,
        "",
        "Mean score: " * _format_float(validation.mean_score; digits = digits),
        "",
        "## Replicate Scores",
        "",
    ]

    for (idx, score) in enumerate(validation.replicate_scores)
        push!(lines, "- replicate $(idx): " * _format_float(score; digits = digits))
    end

    push!(lines, "")
    push!(lines, "## Pass Rates")
    push!(lines, "")
    for key in sort(collect(keys(validation.pass_rates)); by = string)
        push!(lines, "- `" * string(key) * "`: " * _format_float(validation.pass_rates[key]; digits = digits))
    end

    push!(lines, "")
    push!(lines, "## Mean Penalties")
    push!(lines, "")
    for key in sort(collect(keys(validation.mean_penalties)); by = string)
        push!(lines, "- `" * string(key) * "`: " * _format_float(validation.mean_penalties[key]; digits = digits))
    end

    return join(lines, "\n") * "\n"
end

function format_posterior_summary(
    summary::PosteriorSummary;
    digits::Integer = 4,
    title::AbstractString = "Posterior Summary",
)
    max_rhat = isempty(summary.rhat) ? nothing : maximum(values(summary.rhat))
    lines = String[
        "# " * title,
        "",
        "Acceptance rate: " * _format_float(summary.acceptance_rate; digits = digits),
        "Mean score: " * _format_float(summary.mean_score; digits = digits),
        "MAP score: " * _format_float(summary.map_score; digits = digits),
        "MAP log posterior: " * _format_float(summary.map_log_posterior; digits = digits),
        "Retained samples: " * string(summary.n_samples),
        "Burn-in: " * string(summary.burn_in),
    ]

    if max_rhat !== nothing
        push!(lines, "Max R-hat: " * _format_float(max_rhat; digits = digits))
    end

    append!(lines, [
        "",
        "## Posterior Means and Intervals",
        "",
    ])

    for key in sort(collect(keys(summary.posterior_mean)); by = string)
        lower, upper = summary.posterior_interval[key]
        interval_text = "[" * _format_float(lower; digits = digits) * ", " * _format_float(upper; digits = digits) * "]"
        rhat_text = haskey(summary.rhat, key) ? " (R-hat " * _format_float(summary.rhat[key]; digits = digits) * ")" : ""
        push!(lines, "- `" * string(key) * "`: " * _format_float(summary.posterior_mean[key]; digits = digits) * " " * interval_text * rhat_text)
    end

    push!(lines, "")
    push!(lines, "## MAP Parameters")
    push!(lines, "")

    for key in sort(collect(keys(summary.map_parameters)); by = string)
        push!(lines, "- `" * string(key) * "`: " * _format_float(summary.map_parameters[key]; digits = digits))
    end

    return join(lines, "\n") * "\n"
end

function write_posterior_summary_csv(path::AbstractString, summary::PosteriorSummary)
    header = ["parameter", "posterior_mean", "lower_credible", "upper_credible", "rhat"]

    open(path, "w") do io
        write(io, join(header, ",") * "\n")
        for key in sort(collect(keys(summary.posterior_mean)); by = string)
            lower, upper = summary.posterior_interval[key]
            row = [
                key,
                summary.posterior_mean[key],
                lower,
                upper,
                get(summary.rhat, key, ""),
            ]
            write(io, join((_csv_escape(value) for value in row), ","))
            write(io, "\n")
        end
    end

    return path
end

function _csv_escape(value)
    text = string(value)
    escaped = replace(text, "\"" => "\"\"")
    if occursin(',', escaped) || occursin('\n', escaped) || occursin('"', escaped)
        return "\"" * escaped * "\""
    end
    return escaped
end

function write_summary_csv(path::AbstractString, summary::Dict{Symbol,Float64}; label::AbstractString = "cohort")
    open(path, "w") do io
        write(io, "label,metric,value\n")
        for key in _summary_key_order(summary)
            write(io, join((_csv_escape(label), _csv_escape(string(key)), _csv_escape(summary[key])), ","))
            write(io, "\n")
        end
    end

    return path
end

function write_population_csv(path::AbstractString, results::AbstractVector{<:EmbryoResult}; gap_threshold::Real = 0.15)
    header = [
        "embryo_id",
        "solver",
        "trisomy21",
        "phenotype",
        "final_C",
        "final_D",
        "final_G",
        "final_E",
        "final_A",
        "linear_score",
        "nonlinear_score",
        "probability",
        "closure_time",
        "closed",
        "terminal_shear",
        "peak_shear",
    ]

    open(path, "w") do io
        write(io, join(header, ",") * "\n")
        for (idx, result) in enumerate(results)
            final_state = result.trajectory.states[end]
            row = [
                idx,
                result.solver,
                result.genotype.trisomy21,
                result.phenotype,
                final_state.C,
                final_state.D,
                final_state.G,
                final_state.E,
                final_state.A,
                result.scores.linear,
                result.scores.nonlinear,
                result.scores.probability,
                closure_time(result.trajectory; gap_threshold = gap_threshold),
                trajectory_closes(result.trajectory; gap_threshold = gap_threshold),
                terminal_shear(result),
                peak_shear(result),
            ]
            write(io, join((_csv_escape(value) for value in row), ","))
            write(io, "\n")
        end
    end

    return path
end

function write_calibration_history_csv(path::AbstractString, calibration::CalibrationResult)
    penalty_names = sort!(collect(keys(calibration.history[1][:penalties])); by = string)
    has_log_prior = haskey(calibration.history[1], :log_prior)
    has_log_posterior = haskey(calibration.history[1], :log_posterior)
    has_accepted = haskey(calibration.history[1], :accepted)

    open(path, "w") do io
        header = ["iteration", "score"]
        has_log_prior && push!(header, "log_prior")
        has_log_posterior && push!(header, "log_posterior")
        has_accepted && push!(header, "accepted")
        append!(header, [string(name) for name in penalty_names])
        write(io, join(header, ",") * "\n")

        for item in calibration.history
            row = Any[item[:iteration], item[:score]]
            has_log_prior && push!(row, item[:log_prior])
            has_log_posterior && push!(row, item[:log_posterior])
            has_accepted && push!(row, item[:accepted])
            append!(row, [item[:penalties][name] for name in penalty_names])
            write(io, join((_csv_escape(value) for value in row), ","))
            write(io, "\n")
        end
    end

    return path
end

function write_markdown_report(
    path::AbstractString;
    title::AbstractString = "AVSD Analysis Report",
    summary::Union{Nothing,Dict{Symbol,Float64}} = nothing,
    calibration::Union{Nothing,<:CalibrationEvaluation} = nothing,
    posterior::Union{Nothing,<:PosteriorSummary} = nothing,
    validation::Union{Nothing,<:ValidationResult} = nothing,
    sensitivity::Union{Nothing,<:SensitivityResult} = nothing,
    digits::Integer = 4,
)
    lines = String["# " * title, ""]

    if summary !== nothing
        push!(lines, "## Cohort Summary")
        push!(lines, "")
        push!(lines, "```text")
        push!(lines, chomp(format_summary_table(summary; digits = digits)))
        push!(lines, "```")
        push!(lines, "")
    end

    if calibration !== nothing
        push!(lines, format_calibration_report(calibration; digits = digits, title = "Calibration"))
        push!(lines, "")
    end

    if posterior !== nothing
        push!(lines, format_posterior_summary(posterior; digits = digits, title = "Posterior Summary"))
        push!(lines, "")
    end

    if validation !== nothing
        push!(lines, format_validation_report(validation; digits = digits, title = "Validation"))
        push!(lines, "")
    end

    if sensitivity !== nothing
        push!(lines, "## Sensitivity")
        push!(lines, "")
        push!(lines, "Baseline score: " * _format_float(sensitivity.baseline_score; digits = digits))
        push!(lines, "")
        push!(lines, "### Local Effects")
        push!(lines, "")
        for (name, value) in sort(collect(sensitivity.local_effects); by = item -> abs(last(item)), rev = true)
            push!(lines, "- `" * string(name) * "`: " * _format_float(value; digits = digits))
        end
    end

    open(path, "w") do io
        write(io, join(lines, "\n"))
        write(io, "\n")
    end

    return path
end

function generate_results_bundle(
    output_dir::AbstractString;
    params::ModelParameters{T} = default_parameters(),
    targets::Union{Nothing,CalibrationTargets{T}} = nothing,
    initial_state::PhysicalState{T} = default_initial_state(),
    parameter_names::AbstractVector{Symbol} = [:alpha_EMT, :alpha_DMP, :k_G, :k_E, :k_A],
    sensitivity_parameters::AbstractVector{Symbol} = collect(parameter_names),
    population_size::Integer = 96,
    calibration_candidates::Integer = 24,
    validation_replicates::Integer = 5,
    trisomy_fraction::Real = 0.25,
    solver::Symbol = :stochastic_heun,
    sample_strategy::Symbol = :stratified,
    burn_in::Integer = 4,
    n_chains::Integer = 4,
    trace_parameters::AbstractVector{Symbol} = Symbol[],
    seed::Integer = 2026,
) where {T<:Real}
    targets = isnothing(targets) ? default_calibration_targets(; T = T) : targets
    mkpath(output_dir)

    calibration = n_chains == 1 ?
        fit_parameters(
            params = params,
            targets = targets,
            parameter_names = parameter_names,
            n_candidates = calibration_candidates,
            n_embryos = population_size,
            seed = seed,
            initial_state = initial_state,
            solver = solver,
        ) :
        run_calibration_chains(
            params = params,
            targets = targets,
            parameter_names = parameter_names,
            n_candidates = calibration_candidates,
            n_embryos = population_size,
            seed = seed,
            initial_state = initial_state,
            solver = solver,
            n_chains = n_chains,
        )

    fit = calibration isa MultiChainCalibration ? calibration.chains[calibration.best_chain_index] : calibration
    best_params = calibration isa MultiChainCalibration ? calibration.best_parameters : fit.best_parameters
    posterior = posterior_summary(calibration; burn_in = burn_in)
    validation = validate_parameters(
        best_params;
        targets = targets,
        n_replicates = validation_replicates,
        n_embryos = population_size,
        seed = seed + 1_000,
        initial_state = initial_state,
        solver = solver,
    )
    sensitivity = global_sensitivity(
        sensitivity_parameters;
        params = best_params,
        targets = targets,
        n_embryos = population_size,
        seed = seed + 2_000,
        initial_state = initial_state,
        solver = solver,
    )

    mixed_results = simulate_population(
        population_size;
        params = best_params,
        initial_state = initial_state,
        seed = seed + 3_000,
        trisomy_fraction = trisomy_fraction,
        solver = solver,
        sample_strategy = sample_strategy,
    )
    euploid_results = simulate_population(
        population_size;
        params = best_params,
        initial_state = initial_state,
        seed = seed + 4_000,
        trisomy_fraction = 0.0,
        solver = solver,
        sample_strategy = sample_strategy,
    )
    trisomy_results = simulate_population(
        population_size;
        params = best_params,
        initial_state = initial_state,
        seed = seed + 5_000,
        trisomy_fraction = 1.0,
        solver = solver,
        sample_strategy = sample_strategy,
    )

    mixed_summary = summary_metrics(mixed_results)
    euploid_summary = summary_metrics(euploid_results)
    trisomy_summary = summary_metrics(trisomy_results)

    paths = Dict{Symbol,String}()
    paths[:population_mixed_csv] = write_population_csv(joinpath(output_dir, "population_mixed.csv"), mixed_results)
    paths[:population_euploid_csv] = write_population_csv(joinpath(output_dir, "population_euploid.csv"), euploid_results)
    paths[:population_trisomy21_csv] = write_population_csv(joinpath(output_dir, "population_trisomy21.csv"), trisomy_results)
    paths[:summary_mixed_csv] = write_summary_csv(joinpath(output_dir, "summary_mixed.csv"), mixed_summary; label = "mixed")
    paths[:summary_euploid_csv] = write_summary_csv(joinpath(output_dir, "summary_euploid.csv"), euploid_summary; label = "euploid")
    paths[:summary_trisomy21_csv] = write_summary_csv(joinpath(output_dir, "summary_trisomy21.csv"), trisomy_summary; label = "trisomy21")
    if calibration isa MultiChainCalibration
        for (idx, chain) in enumerate(calibration.chains)
            paths[Symbol("calibration_history_chain_$(idx)_csv")] = write_calibration_history_csv(
                joinpath(output_dir, "calibration_history_chain_$(idx).csv"),
                chain,
            )
        end
    else
        paths[:calibration_history_csv] = write_calibration_history_csv(joinpath(output_dir, "calibration_history.csv"), fit)
    end
    paths[:posterior_summary_csv] = write_posterior_summary_csv(joinpath(output_dir, "posterior_summary.csv"), posterior)
    paths[:gap_trajectories_svg] = write_trajectory_svg(
        joinpath(output_dir, "gap_trajectories.svg"),
        [result.trajectory for result in mixed_results];
        variable = :G,
        title = "Mixed Cohort Gap Trajectories",
    )
    paths[:ecm_trajectories_svg] = write_trajectory_svg(
        joinpath(output_dir, "ecm_trajectories.svg"),
        [result.trajectory for result in mixed_results];
        variable = :E,
        title = "Mixed Cohort ECM Trajectories",
    )
    paths[:calibration_history_svg] = write_calibration_history_svg(joinpath(output_dir, "calibration_history.svg"), fit)
    paths[:sensitivity_svg] = write_sensitivity_svg(joinpath(output_dir, "sensitivity.svg"), sensitivity)
    paths[:cohort_comparison_svg] = write_cohort_comparison_svg(
        joinpath(output_dir, "cohort_comparison.svg"),
        Dict(
            :euploid => euploid_summary,
            :mixed => mixed_summary,
            :trisomy21 => trisomy_summary,
        ),
    )
    if calibration isa MultiChainCalibration
        paths[:posterior_trace_svg] = write_posterior_trace_svg(
            joinpath(output_dir, "posterior_traces.svg"),
            calibration;
            parameter_names = isempty(trace_parameters) ? Symbol[] : collect(trace_parameters),
        )
    end
    paths[:report_md] = write_markdown_report(
        joinpath(output_dir, "report.md");
        title = "AVSD Research Results",
        summary = mixed_summary,
        calibration = calibration isa MultiChainCalibration ? calibration.best_evaluation : fit.best_evaluation,
        posterior = posterior,
        validation = validation,
        sensitivity = sensitivity,
    )

    open(joinpath(output_dir, "posterior_summary.md"), "w") do io
        write(io, format_posterior_summary(posterior))
    end
    paths[:posterior_summary_md] = joinpath(output_dir, "posterior_summary.md")

    open(joinpath(output_dir, "run_manifest.txt"), "w") do io
        write(io, "seed=$(seed)\n")
        write(io, "solver=$(solver)\n")
        write(io, "sample_strategy=$(sample_strategy)\n")
        write(io, "population_size=$(population_size)\n")
        write(io, "calibration_candidates=$(calibration_candidates)\n")
        write(io, "validation_replicates=$(validation_replicates)\n")
        write(io, "n_chains=$(n_chains)\n")
        write(io, "trisomy_fraction=$(trisomy_fraction)\n")
        write(io, "burn_in=$(burn_in)\n")
        write(io, "parameter_names=$(join(string.(parameter_names), ','))\n")
    end
    paths[:run_manifest] = joinpath(output_dir, "run_manifest.txt")

    return paths
end

function run_research_pipeline(
    config::ResearchRunConfig;
    params::ModelParameters = load_parameters(),
    initial_state::PhysicalState = load_initial_state(),
)
    targets = load_calibration_targets(config.calibration_targets_path)
    return generate_results_bundle(
        config.output_dir;
        params = params,
        targets = targets,
        initial_state = initial_state,
        parameter_names = config.parameter_names,
        sensitivity_parameters = config.sensitivity_parameters,
        population_size = config.population_size,
        calibration_candidates = config.calibration_candidates,
        validation_replicates = config.validation_replicates,
        trisomy_fraction = config.trisomy_fraction,
        solver = config.solver,
        sample_strategy = config.sample_strategy,
        burn_in = config.burn_in,
        n_chains = config.n_chains,
        trace_parameters = config.trace_parameters,
        seed = config.seed,
    )
end

function run_research_pipeline(path::AbstractString = joinpath(project_root(), "config", "research_run.yaml"))
    return run_research_pipeline(load_research_run_config(path))
end
