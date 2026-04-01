struct GenotypeSample{T<:Real}
    trisomy21::Bool
    latent_axes::Vector{T}
end

function sample_genotype(
    rng::AbstractRNG,
    genetics::GeneticParameters{T};
    trisomy21::Bool = false,
) where {T<:Real}
    z = [randn(rng) * genetics.latent_sd for _ in 1:genetics.latent_dim]
    return GenotypeSample{T}(trisomy21, z)
end

function apply_genotype(
    rng::AbstractRNG,
    params::ModelParameters{T},
    genotype::GenotypeSample{T},
) where {T<:Real}
    sampled_rates = copy(params.rates)

    for (key, baseline_value) in params.rates
        loading = get(params.genetics.loadings, key, fill(zero(T), params.genetics.latent_dim))
        n_shared = min(length(loading), length(genotype.latent_axes))
        latent_shift = n_shared == 0 ? zero(T) : dot(view(loading, 1:n_shared), view(genotype.latent_axes, 1:n_shared))
        trisomy_shift = genotype.trisomy21 ? get(params.genetics.trisomy_effects, key, zero(T)) : zero(T)
        stochastic_shift = randn(rng) * params.genetics.parameter_noise_sd
        sampled_rates[key] = baseline_value * exp(trisomy_shift + latent_shift + stochastic_shift)
    end

    return with_rates(params; rates = sampled_rates), sampled_rates
end
