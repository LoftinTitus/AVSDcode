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

    open(path, "w") do io
        header = vcat(["iteration", "score"], [string(name) for name in penalty_names])
        write(io, join(header, ",") * "\n")

        for item in calibration.history
            row = Any[item[:iteration], item[:score]]
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
