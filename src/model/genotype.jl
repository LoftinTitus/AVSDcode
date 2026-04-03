struct GenotypeSample{T<:Real}
    trisomy21::Bool
    latent_axes::Vector{T}
end

function _inverse_standard_normal_cdf(p::Real)
    0.0 < p < 1.0 || error("Probability must lie strictly between 0 and 1")

    a = (
        -3.969683028665376e+01,
        2.209460984245205e+02,
        -2.759285104469687e+02,
        1.383577518672690e+02,
        -3.066479806614716e+01,
        2.506628277459239e+00,
    )
    b = (
        -5.447609879822406e+01,
        1.615858368580409e+02,
        -1.556989798598866e+02,
        6.680131188771972e+01,
        -1.328068155288572e+01,
    )
    c = (
        -7.784894002430293e-03,
        -3.223964580411365e-01,
        -2.400758277161838e+00,
        -2.549732539343734e+00,
        4.374664141464968e+00,
        2.938163982698783e+00,
    )
    d = (
        7.784695709041462e-03,
        3.224671290700398e-01,
        2.445134137142996e+00,
        3.754408661907416e+00,
    )

    plow = 0.02425
    phigh = 1 - plow

    if p < plow
        q = sqrt(-2 * log(p))
        return (((((c[1] * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) * q + c[6]) /
               ((((d[1] * q + d[2]) * q + d[3]) * q + d[4]) * q + 1)
    elseif p <= phigh
        q = p - 0.5
        r = q * q
        return (((((a[1] * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * r + a[6]) * q /
               (((((b[1] * r + b[2]) * r + b[3]) * r + b[4]) * r + b[5]) * r + 1)
    else
        q = sqrt(-2 * log(1 - p))
        return -(((((c[1] * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) * q + c[6]) /
                ((((d[1] * q + d[2]) * q + d[3]) * q + d[4]) * q + 1)
    end
end

function sample_genotype(
    rng::AbstractRNG,
    genetics::GeneticParameters{T};
    trisomy21::Bool = false,
    latent_axes::Union{Nothing,AbstractVector{<:Real}} = nothing,
) where {T<:Real}
    z = if isnothing(latent_axes)
        [randn(rng) * genetics.latent_sd for _ in 1:genetics.latent_dim]
    else
        [T(value) for value in latent_axes]
    end
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
