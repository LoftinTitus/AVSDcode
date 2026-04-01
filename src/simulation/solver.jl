struct EmbryoResult{T<:Real}
    trajectory::Trajectory{T}
    parameters::ModelParameters{T}
    sampled_rates::Dict{Symbol,T}
    genotype::GenotypeSample{T}
    solver::Symbol
    scores::NamedTuple{(:linear, :nonlinear, :probability, :gap),Tuple{T,T,T,T}}
    phenotype::Bool
end

function simulate_trajectory(
    rng::AbstractRNG,
    initial_state::PhysicalState{T},
    params::ModelParameters{T};
    solver::Symbol = :stochastic_heun,
) where {T<:Real}
    if solver == :stochastic_heun
        return stochastic_heun(rng, initial_state, params)
    elseif solver == :euler_maruyama
        return euler_maruyama(rng, initial_state, params)
    elseif solver == :bdf1
        return bdf1_stiff(initial_state, params)
    else
        error("Unsupported solver $(solver). Choose :bdf1, :euler_maruyama, or :stochastic_heun.")
    end
end

function solve_embryo(
    rng::AbstractRNG;
    params::ModelParameters{T} = default_parameters(),
    initial_state::PhysicalState{T} = default_initial_state(),
    trisomy21::Bool = false,
    solver::Symbol = :stochastic_heun,
) where {T<:Real}
    genotype = sample_genotype(rng, params.genetics; trisomy21 = trisomy21)
    sampled_params, sampled_rates = apply_genotype(rng, params, genotype)
    trajectory = simulate_trajectory(rng, initial_state, sampled_params; solver = solver)
    final_state = trajectory.states[end]
    scores = phenotype_scores(final_state, sampled_params)
    phenotype = phenotype_label(final_state, sampled_params)
    return EmbryoResult(trajectory, sampled_params, sampled_rates, genotype, solver, scores, phenotype)
end
