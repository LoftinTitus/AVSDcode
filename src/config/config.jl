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

function _float_dict(source::Dict{String,Any})
    return Dict(Symbol(key) => float(value) for (key, value) in source)
end

function _vector_float_dict(source::Dict{String,Any})
    return Dict(Symbol(key) => Float64[float(v) for v in value] for (key, value) in source)
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
