flow_rate(t::Real, params::ModelParameters) = params.hemo.Q0 * exp(params.hemo.k_Q * t)

function effective_radius(state::PhysicalState, params::ModelParameters)
    gap_term = normalized_components(state, params).g
    return max(params.hemo.R0 + params.hemo.kappa_R * gap_term, STATE_EPS)
end

function shear_stress(t::Real, state::PhysicalState, params::ModelParameters)
    numerator = 4 * params.hemo.mu * flow_rate(t, params)
    denominator = pi * effective_radius(state, params)^3
    return numerator / denominator
end

normalized_shear(t::Real, state::PhysicalState, params::ModelParameters) = shear_stress(t, state, params) / params.hemo.tau_star
