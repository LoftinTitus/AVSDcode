function parse_scalar(token::AbstractString)
    stripped = strip(token)
    isempty(stripped) && return ""

    lowered = lowercase(stripped)
    lowered == "true" && return true
    lowered == "false" && return false

    if startswith(stripped, "[") && endswith(stripped, "]")
        inner = strip(stripped[2:end-1])
        isempty(inner) && return Any[]
        return [parse_scalar(piece) for piece in split(inner, ",")]
    end

    if (startswith(stripped, "\"") && endswith(stripped, "\"")) || (startswith(stripped, "'") && endswith(stripped, "'"))
        return stripped[2:end-1]
    end

    try
        return parse(Int, stripped)
    catch
    end

    try
        return parse(Float64, stripped)
    catch
    end

    return stripped
end

struct ResearchRunConfig
    output_dir::String
    calibration_targets_path::String
    population_size::Int
    calibration_candidates::Int
    proposal_scale::Float64
    prior_mix::Float64
    refinement_rounds::Int
    refinement_candidates::Int
    validation_replicates::Int
    n_chains::Int
    burn_in::Int
    seed::Int
    trisomy_fraction::Float64
    solver::Symbol
    sample_strategy::Symbol
    parameter_names::Vector{Symbol}
    sensitivity_parameters::Vector{Symbol}
    trace_parameters::Vector{Symbol}
end

function load_yaml_like(path::AbstractString)
    root = Dict{String,Any}()
    dict_stack = [root]
    indent_stack = [-1]

    for raw_line in eachline(path)
        line = replace(raw_line, '\t' => "    ")
        stripped = strip(line)

        if isempty(stripped) || startswith(stripped, "#")
            continue
        end

        indent = length(line) - length(lstrip(line))

        while indent <= last(indent_stack) && length(dict_stack) > 1
            pop!(dict_stack)
            pop!(indent_stack)
        end

        current = last(dict_stack)
        parts = split(stripped, ":"; limit = 2)
        length(parts) == 2 || error("Invalid YAML-like line: $raw_line")

        key = strip(parts[1])
        value = strip(parts[2])

        if isempty(value)
            child = Dict{String,Any}()
            current[key] = child
            push!(dict_stack, child)
            push!(indent_stack, indent)
        else
            current[key] = parse_scalar(value)
        end
    end

    return root
end

project_root() = dirname(dirname(@__DIR__))

function _resolve_project_path(path::AbstractString)
    return isabspath(path) ? path : normpath(joinpath(project_root(), path))
end

function _float_dict(source::Dict{String,Any})
    return Dict(Symbol(key) => float(value) for (key, value) in source)
end

function _vector_float_dict(source::Dict{String,Any})
    return Dict(Symbol(key) => Float64[float(v) for v in value] for (key, value) in source)
end

function _symbol_vector(values)
    return Symbol[Symbol(string(value)) for value in values]
end

function load_parameters(path::AbstractString = joinpath(project_root(), "config", "parameters.yaml"))
    raw = load_yaml_like(path)

    windows = (
        emt = TimingWindow(
            float(raw["windows"]["emt"]["t_on"]),
            float(raw["windows"]["emt"]["t_off"]),
            float(raw["windows"]["emt"]["smoothness"]),
        ),
        dmp = TimingWindow(
            float(raw["windows"]["dmp"]["t_on"]),
            float(raw["windows"]["dmp"]["t_off"]),
            float(raw["windows"]["dmp"]["smoothness"]),
        ),
        close = TimingWindow(
            float(raw["windows"]["close"]["t_on"]),
            float(raw["windows"]["close"]["t_off"]),
            float(raw["windows"]["close"]["smoothness"]),
        ),
        avc = TimingWindow(
            float(raw["windows"]["avc"]["t_on"]),
            float(raw["windows"]["avc"]["t_off"]),
            float(raw["windows"]["avc"]["smoothness"]),
        ),
    )

    hemo = HemodynamicsParameters(
        float(raw["hemodynamics"]["Q0"]),
        float(raw["hemodynamics"]["k_Q"]),
        float(raw["hemodynamics"]["R0"]),
        float(raw["hemodynamics"]["kappa_R"]),
        float(raw["hemodynamics"]["mu"]),
        float(raw["hemodynamics"]["tau_star"]),
    )

    activation = ActivationParameters(
        float(raw["activation"]["hill_n"]),
        float(raw["activation"]["hill_K"]),
        float(raw["activation"]["N_crit"]),
    )

    phenotype = PhenotypeParameters(
        float(raw["phenotype"]["w_C"]),
        float(raw["phenotype"]["w_D"]),
        float(raw["phenotype"]["w_E"]),
        float(raw["phenotype"]["w_A"]),
        float(raw["phenotype"]["threshold"]),
        float(raw["phenotype"]["slope"]),
    )

    genetics = GeneticParameters(
        Int(raw["genetics"]["latent_dim"]),
        float(raw["genetics"]["latent_sd"]),
        float(raw["genetics"]["parameter_noise_sd"]),
        _float_dict(raw["genetics"]["trisomy_effects"]),
        _vector_float_dict(raw["genetics"]["loadings"]),
    )

    return ModelParameters(
        float(raw["simulation"]["T"]),
        float(raw["simulation"]["dt"]),
        windows,
        hemo,
        activation,
        phenotype,
        _float_dict(raw["dynamics"]),
        genetics,
    )
end

function load_initial_state(path::AbstractString = joinpath(project_root(), "config", "parameters.yaml"))
    raw = load_yaml_like(path)
    return PhysicalState(
        float(raw["initial_state"]["C"]),
        float(raw["initial_state"]["D"]),
        float(raw["initial_state"]["G"]),
        float(raw["initial_state"]["E"]),
        float(raw["initial_state"]["A"]),
    )
end

load_priors(path::AbstractString = joinpath(project_root(), "config", "priors.yaml")) = load_yaml_like(path)

function load_calibration_targets(
    path::AbstractString = joinpath(project_root(), "data", "calibration_targets.csv");
    T::Type{<:Real} = Float64,
)
    euploid_prevalence = nothing
    trisomy21_prevalence = nothing
    closure_time = nothing
    shear_range = nothing
    max_gap = nothing
    weights = Dict{Symbol,T}()

    for raw_line in eachline(path)
        stripped = strip(raw_line)
        if isempty(stripped) || startswith(stripped, "#") || lowercase(stripped) == "metric,cohort,lower,upper,weight"
            continue
        end

        fields = split(stripped, ",")
        length(fields) == 5 || error("Invalid calibration target row: $raw_line")

        metric = strip(fields[1])
        cohort = strip(fields[2])
        lower = T(parse(Float64, strip(fields[3])))
        upper = T(parse(Float64, strip(fields[4])))
        weight = T(parse(Float64, strip(fields[5])))

        if metric == "prevalence" && cohort == "euploid"
            euploid_prevalence = (lower, upper)
            weights[:euploid_prevalence] = weight
        elseif metric == "prevalence" && cohort == "trisomy21"
            trisomy21_prevalence = (lower, upper)
            weights[:trisomy21_prevalence] = weight
        elseif metric == "closure_time"
            closure_time = (lower, upper)
            weights[:closure_time] = weight
        elseif metric == "terminal_shear"
            shear_range = (lower, upper)
            weights[:terminal_shear] = weight
        elseif metric == "gap"
            max_gap = upper
            weights[:gap] = weight
        else
            error("Unsupported calibration target metric/cohort combination: $(metric), $(cohort)")
        end
    end

    isnothing(euploid_prevalence) && error("Missing euploid prevalence target in $(path)")
    isnothing(trisomy21_prevalence) && error("Missing trisomy21 prevalence target in $(path)")
    isnothing(closure_time) && error("Missing closure time target in $(path)")
    isnothing(shear_range) && error("Missing terminal shear target in $(path)")
    isnothing(max_gap) && error("Missing gap target in $(path)")

    return CalibrationTargets(
        euploid_prevalence,
        trisomy21_prevalence,
        closure_time,
        shear_range,
        max_gap,
        weights,
    )
end

function load_research_run_config(path::AbstractString = joinpath(project_root(), "config", "research_run.yaml"))
    raw = load_yaml_like(path)
    run = haskey(raw, "run") ? raw["run"] : raw

    parameter_names = _symbol_vector(get(run, "parameter_names", Any["alpha_EMT", "alpha_DMP", "k_G", "k_E", "k_A"]))
    sensitivity_parameters = _symbol_vector(get(run, "sensitivity_parameters", Any[string(name) for name in parameter_names]))
    trace_parameters = _symbol_vector(get(run, "trace_parameters", Any[string(name) for name in parameter_names[1:min(end, 3)]]))

    return ResearchRunConfig(
        _resolve_project_path(string(get(run, "output_dir", "research_outputs/default_run"))),
        _resolve_project_path(string(get(run, "calibration_targets_path", joinpath("data", "calibration_targets.csv")))),
        Int(get(run, "population_size", 96)),
        Int(get(run, "calibration_candidates", 24)),
        float(get(run, "proposal_scale", 0.20)),
        float(get(run, "prior_mix", 0.15)),
        Int(get(run, "refinement_rounds", 0)),
        Int(get(run, "refinement_candidates", max(cld(Int(get(run, "calibration_candidates", 24)), 2), 1))),
        Int(get(run, "validation_replicates", 5)),
        Int(get(run, "n_chains", 4)),
        Int(get(run, "burn_in", 8)),
        Int(get(run, "seed", 2026)),
        float(get(run, "trisomy_fraction", 0.25)),
        Symbol(string(get(run, "solver", "stochastic_heun"))),
        Symbol(string(get(run, "sample_strategy", "stratified"))),
        parameter_names,
        sensitivity_parameters,
        trace_parameters,
    )
end
