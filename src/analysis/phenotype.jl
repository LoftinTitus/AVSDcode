function phenotype_scores(state::PhysicalState, params::ModelParameters)
    comps = normalized_components(state, params)

    linear = (
        params.phenotype.w_C * comps.c +
        params.phenotype.w_D * comps.d +
        params.phenotype.w_E * state.E +
        params.phenotype.w_A * state.A
    )

    nonlinear = (
        max(comps.c, STATE_EPS)^params.phenotype.w_C *
        max(comps.d, STATE_EPS)^params.phenotype.w_D *
        max(state.E, STATE_EPS)^params.phenotype.w_E *
        max(state.A, STATE_EPS)^params.phenotype.w_A
    )

    probability = logistic(params.phenotype.slope * (params.phenotype.threshold - linear))

    return (
        linear = linear,
        nonlinear = nonlinear,
        probability = probability,
        gap = state.G,
    )
end

function phenotype_probability(state::PhysicalState, params::ModelParameters; model::Symbol = :linear)
    scores = phenotype_scores(state, params)
    score = model == :nonlinear ? scores.nonlinear : scores.linear
    return logistic(params.phenotype.slope * (params.phenotype.threshold - score))
end

function phenotype_label(state::PhysicalState, params::ModelParameters; model::Symbol = :linear)
    scores = phenotype_scores(state, params)
    score = model == :nonlinear ? scores.nonlinear : scores.linear
    return score <= params.phenotype.threshold
end
