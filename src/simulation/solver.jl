struct EmbryoResult{T<:Real}
    trajectory::Trajectory{T}
    parameters::ModelParameters{T}
    sampled_rates::Dict{Symbol,T}
    genotype::GenotypeSample{T}
    scores::NamedTuple{(:linear, :nonlinear, :probability, :gap),Tuple{T,T,T,T}}
    phenotype::Bool
end

function solve_embryo(
    rng::AbstractRNG;
    params::ModelParameters{T} = default_parameters(),
    initial_state::PhysicalState{T} = default_initial_state(),
    trisomy21::Bool = false,
) where {T<:Real}
    genotype = sample_genotype(rng, params.genetics; trisomy21 = trisomy21)
    sampled_params, sampled_rates = apply_genotype(rng, params, genotype)
    trajectory = euler_maruyama(rng, initial_state, sampled_params)
    final_state = trajectory.states[end]
    scores = phenotype_scores(final_state, sampled_params)
    phenotype = phenotype_label(final_state, sampled_params)
    return EmbryoResult(trajectory, sampled_params, sampled_rates, genotype, scores, phenotype)
end
