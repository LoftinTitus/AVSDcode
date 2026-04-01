struct TimingWindow{T<:Real}
    t_on::T
    t_off::T
    smoothness::T
end

struct HemodynamicsParameters{T<:Real}
    Q0::T
    k_Q::T
    R0::T
    kappa_R::T
    mu::T
    tau_star::T
end

struct ActivationParameters{T<:Real}
    hill_n::T
    hill_K::T
    N_crit::T
end

struct PhenotypeParameters{T<:Real}
    w_C::T
    w_D::T
    w_E::T
    w_A::T
    threshold::T
    slope::T
end

struct GeneticParameters{T<:Real}
    latent_dim::Int
    latent_sd::T
    parameter_noise_sd::T
    trisomy_effects::Dict{Symbol,T}
    loadings::Dict{Symbol,Vector{T}}
end

struct ModelParameters{T<:Real}
    T::T
    dt::T
    windows::NamedTuple{(:emt, :dmp, :close, :avc),Tuple{TimingWindow{T},TimingWindow{T},TimingWindow{T},TimingWindow{T}}}
    hemo::HemodynamicsParameters{T}
    activation::ActivationParameters{T}
    phenotype::PhenotypeParameters{T}
    rates::Dict{Symbol,T}
    genetics::GeneticParameters{T}
end

function default_parameters()
    windows = (
        emt = TimingWindow(0.45, 1.55, 0.12),
        dmp = TimingWindow(0.45, 1.55, 0.12),
        close = TimingWindow(1.35, 3.20, 0.18),
        avc = TimingWindow(0.15, 2.10, 0.15),
    )

    hemo = HemodynamicsParameters(1.0, 0.25, 1.1, 0.6, 1.0, 3.0)
    activation = ActivationParameters(3.0, 0.45, 1.05)
    phenotype = PhenotypeParameters(0.35, 0.25, 0.20, 0.20, 0.70, 8.0)

    rates = Dict{Symbol,Float64}(
        :alpha_EMT => 0.80,
        :r_C => 0.55,
        :delta_C => 0.14,
        :sigma_C => 0.08,
        :alpha_DMP => 0.60,
        :r_D => 0.42,
        :delta_D => 0.12,
        :sigma_D => 0.07,
        :k_G => 0.95,
        :omega_C => 0.55,
        :omega_D => 0.45,
        :omega_E => 0.50,
        :omega_A => 0.50,
        :k_N => 0.18,
        :sigma_G => 0.05,
        :k_E => 0.72,
        :gamma_E => 0.10,
        :sigma_E => 0.02,
        :k_A => 0.68,
        :gamma_A => 0.12,
        :sigma_A => 0.02,
        :late_decay => 1.30,
        :C_scale => 1.00,
        :D_scale => 0.80,
        :gap_scale => 1.00,
    )

    genetics = GeneticParameters(
        2,
        0.55,
        0.12,
        Dict(
            :alpha_EMT => -0.08,
            :alpha_DMP => -0.10,
            :k_G => -0.08,
            :k_E => -0.05,
            :k_A => -0.05,
            :k_N => 0.06,
            :sigma_G => 0.05,
        ),
        Dict(
            :alpha_EMT => [-0.10, -0.05],
            :alpha_DMP => [-0.07, -0.03],
            :k_G => [-0.12, -0.04],
            :k_E => [-0.05, 0.04],
            :k_A => [-0.05, -0.02],
            :k_N => [0.08, 0.05],
            :sigma_G => [0.06, 0.03],
        ),
    )

    return ModelParameters(3.5, 0.01, windows, hemo, activation, phenotype, rates, genetics)
end

function with_rates(
    params::ModelParameters{T};
    rates::Dict{Symbol,T} = params.rates,
) where {T<:Real}
    return ModelParameters(params.T, params.dt, params.windows, params.hemo, params.activation, params.phenotype, rates, params.genetics)
end

window_value(t::Real, window::TimingWindow) = timing_window(t, window.t_on, window.t_off, window.smoothness)

rate(params::ModelParameters, key::Symbol) = params.rates[key]

function normalized_components(state::PhysicalState, params::ModelParameters)
    c = state.C / max(rate(params, :C_scale), STATE_EPS)
    d = state.D / max(rate(params, :D_scale), STATE_EPS)
    g = state.G / max(rate(params, :gap_scale), STATE_EPS)
    return (c = c, d = d, g = g)
end

function drift(state::PhysicalState, t::Real, params::ModelParameters)
    comps = normalized_components(state, params)
    N = normalized_shear(t, state, params)

    phi_A = hill_activation(state.A, params.activation.hill_n, params.activation.hill_K)
    phi_E = hill_activation(state.E, params.activation.hill_n, params.activation.hill_K)
    psi_N = hill_activation(N, params.activation.hill_n, params.activation.hill_K)
    psi_C = hill_activation(comps.c, params.activation.hill_n, params.activation.hill_K)

    dC = (
        rate(params, :alpha_EMT) * window_value(t, params.windows.emt) +
        rate(params, :r_C) * state.C * phi_A * phi_E -
        rate(params, :delta_C) * state.C
    )

    dD = (
        rate(params, :alpha_DMP) * window_value(t, params.windows.dmp) * phi_A * late_penalty(t, params.windows.dmp.t_off, rate(params, :late_decay)) +
        rate(params, :r_D) * state.D -
        rate(params, :delta_D) * state.D
    )

    closure_force = (
        rate(params, :omega_C) * comps.c +
        rate(params, :omega_D) * comps.d
    ) * (
        rate(params, :omega_E) * state.E +
        rate(params, :omega_A) * state.A
    )

    dG = (
        -rate(params, :k_G) * window_value(t, params.windows.close) * closure_force +
        rate(params, :k_N) * excess_over_threshold(N, params.activation.N_crit)
    )

    dE = (
        rate(params, :k_E) * (1 - state.E) * psi_N * psi_C -
        rate(params, :gamma_E) * state.E
    )

    dA = (
        rate(params, :k_A) * window_value(t, params.windows.avc) * (1 - state.A) -
        rate(params, :gamma_A) * state.A
    )

    return PhysicalState(dC, dD, dG, dE, dA)
end

function diffusion(state::PhysicalState, ::Real, params::ModelParameters)
    return PhysicalState(
        rate(params, :sigma_C) * state.C,
        rate(params, :sigma_D) * state.D,
        rate(params, :sigma_G),
        rate(params, :sigma_E),
        rate(params, :sigma_A),
    )
end

function transformed_drift(state::PhysicalState, t::Real, params::ModelParameters)
    mu = drift(state, t, params)
    E = clamp(state.E, STATE_EPS, 1 - STATE_EPS)
    A = clamp(state.A, STATE_EPS, 1 - STATE_EPS)
    return TransformedState(
        mu.C / max(state.C, STATE_EPS),
        mu.D / max(state.D, STATE_EPS),
        mu.G / max(state.G, STATE_EPS),
        mu.E / (E * (1 - E)),
        mu.A / (A * (1 - A)),
    )
end

function transformed_diffusion(state::PhysicalState, t::Real, params::ModelParameters)
    sigma = diffusion(state, t, params)
    E = clamp(state.E, STATE_EPS, 1 - STATE_EPS)
    A = clamp(state.A, STATE_EPS, 1 - STATE_EPS)
    return TransformedState(
        sigma.C / max(state.C, STATE_EPS),
        sigma.D / max(state.D, STATE_EPS),
        sigma.G / max(state.G, STATE_EPS),
        sigma.E / (E * (1 - E)),
        sigma.A / (A * (1 - A)),
    )
end

function drift!(du::AbstractVector, state::PhysicalState, t::Real, params::ModelParameters)
    mu = drift(state, t, params)
    du[1] = mu.C
    du[2] = mu.D
    du[3] = mu.G
    du[4] = mu.E
    du[5] = mu.A
    return du
end

function diffusion!(du::AbstractVector, state::PhysicalState, t::Real, params::ModelParameters)
    sigma = diffusion(state, t, params)
    du[1] = sigma.C
    du[2] = sigma.D
    du[3] = sigma.G
    du[4] = sigma.E
    du[5] = sigma.A
    return du
end
