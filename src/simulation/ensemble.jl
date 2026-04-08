struct EmbryoScenario{T<:Real}
    seed::Int
    trisomy21::Bool
    genotype::GenotypeSample{T}
end

struct CalibrationPanel{T<:Real}
    euploid::Vector{EmbryoScenario{T}}
    trisomy21::Vector{EmbryoScenario{T}}
end

function _population_assignment(
    rng::AbstractRNG,
    n_embryos::Integer,
    trisomy_fraction::Real,
    sample_strategy::Symbol,
)
    0.0 <= trisomy_fraction <= 1.0 || error("trisomy_fraction must lie between 0 and 1")

    if sample_strategy == :random
        return [rand(rng) < trisomy_fraction for _ in 1:n_embryos]
    elseif sample_strategy == :stratified
        n_trisomy = round(Int, n_embryos * trisomy_fraction)
        assignments = fill(false, n_embryos)
        if n_trisomy > 0
            for idx in randperm(rng, n_embryos)[1:n_trisomy]
                assignments[idx] = true
            end
        end
        return assignments
    else
        error("Unsupported sample strategy $(sample_strategy). Choose :random or :stratified.")
    end
end

function _latin_hypercube_latent_axes(
    rng::AbstractRNG,
    n_embryos::Integer,
    genetics::GeneticParameters{T},
) where {T<:Real}
    axes = Matrix{T}(undef, n_embryos, genetics.latent_dim)
    n_embryos == 0 && return axes

    for dim in 1:genetics.latent_dim
        strata = randperm(rng, n_embryos)
        for idx in 1:n_embryos
            u = (strata[idx] - rand(rng)) / n_embryos
            u = clamp(u, eps(Float64), 1 - eps(Float64))
            axes[idx, dim] = T(_inverse_standard_normal_cdf(u) * genetics.latent_sd)
        end
    end

    return axes
end

function build_population_panel(
    n_embryos::Integer;
    params::ModelParameters{T} = default_parameters(),
    seed::Integer = 42,
    trisomy_fraction::Real = 0.0,
    sample_strategy::Symbol = :stratified,
) where {T<:Real}
    rng = MersenneTwister(seed)
    scenarios = Vector{EmbryoScenario{T}}(undef, n_embryos)
    trisomy_assignments = _population_assignment(rng, n_embryos, trisomy_fraction, sample_strategy)

    latent_samples = if sample_strategy == :stratified
        _latin_hypercube_latent_axes(rng, n_embryos, params.genetics)
    else
        Matrix{T}(undef, 0, 0)
    end

    for idx in 1:n_embryos
        genotype = if sample_strategy == :stratified
            sample_genotype(
                rng,
                params.genetics;
                trisomy21 = trisomy_assignments[idx],
                latent_axes = vec(latent_samples[idx, :]),
            )
        else
            sample_genotype(
                rng,
                params.genetics;
                trisomy21 = trisomy_assignments[idx],
            )
        end

        scenarios[idx] = EmbryoScenario{T}(
            rand(rng, 1:typemax(Int)),
            trisomy_assignments[idx],
            genotype,
        )
    end

    return scenarios
end

function build_calibration_panel(
    n_embryos::Integer;
    params::ModelParameters{T} = default_parameters(),
    seed::Integer = 42,
    sample_strategy::Symbol = :stratified,
) where {T<:Real}
    return CalibrationPanel(
        build_population_panel(
            n_embryos;
            params = params,
            seed = seed,
            trisomy_fraction = 0.0,
            sample_strategy = sample_strategy,
        ),
        build_population_panel(
            n_embryos;
            params = params,
            seed = seed + 10_000,
            trisomy_fraction = 1.0,
            sample_strategy = sample_strategy,
        ),
    )
end

function simulate_population(
    scenarios::AbstractVector{<:EmbryoScenario{T}};
    params::ModelParameters{T} = default_parameters(),
    initial_state::PhysicalState{T} = default_initial_state(),
    solver::Symbol = :stochastic_heun,
) where {T<:Real}
    results = Vector{EmbryoResult{T}}(undef, length(scenarios))

    for (idx, scenario) in enumerate(scenarios)
        results[idx] = solve_embryo(
            MersenneTwister(scenario.seed);
            params = params,
            initial_state = initial_state,
            trisomy21 = scenario.trisomy21,
            genotype = scenario.genotype,
            solver = solver,
        )
    end

    return results
end

function simulate_population(
    n_embryos::Integer;
    params::ModelParameters{T} = default_parameters(),
    initial_state::PhysicalState{T} = default_initial_state(),
    seed::Integer = 42,
    trisomy_fraction::Real = 0.0,
    solver::Symbol = :stochastic_heun,
    sample_strategy::Symbol = :stratified,
) where {T<:Real}
    scenarios = build_population_panel(
        n_embryos;
        params = params,
        seed = seed,
        trisomy_fraction = trisomy_fraction,
        sample_strategy = sample_strategy,
    )
    return simulate_population(
        scenarios;
        params = params,
        initial_state = initial_state,
        solver = solver,
    )
end
