struct PhysicalState{T<:Real}
    C::T
    D::T
    G::T
    E::T
    A::T
end

struct TransformedState{T<:Real}
    xC::T
    xD::T
    xG::T
    xE::T
    xA::T
end

const POSITIVE_STATE_MIN = log(STATE_EPS)
const POSITIVE_STATE_MAX = log(8.0)
const LOGIT_STATE_MIN = -12.0
const LOGIT_STATE_MAX = 12.0

default_initial_state() = PhysicalState(0.22, 0.12, 1.00, 0.18, 0.24)

function bounded_transformed_value(
    value::Real;
    lower::Real = POSITIVE_STATE_MIN,
    upper::Real = POSITIVE_STATE_MAX,
    fallback::Real = 0.0,
)
    raw = float(value)
    isnan(raw) && return float(fallback)
    raw == Inf && return float(upper)
    raw == -Inf && return float(lower)
    return clamp(raw, float(lower), float(upper))
end

function bounded_transformed_state(state::TransformedState)
    return TransformedState(
        bounded_transformed_value(state.xC; lower = POSITIVE_STATE_MIN, upper = POSITIVE_STATE_MAX),
        bounded_transformed_value(state.xD; lower = POSITIVE_STATE_MIN, upper = POSITIVE_STATE_MAX),
        bounded_transformed_value(state.xG; lower = POSITIVE_STATE_MIN, upper = POSITIVE_STATE_MAX),
        bounded_transformed_value(state.xE; lower = LOGIT_STATE_MIN, upper = LOGIT_STATE_MAX),
        bounded_transformed_value(state.xA; lower = LOGIT_STATE_MIN, upper = LOGIT_STATE_MAX),
    )
end

function transformed_state(state::PhysicalState)
    return TransformedState(
        log(max(float(state.C), STATE_EPS)),
        log(max(float(state.D), STATE_EPS)),
        log(max(float(state.G), STATE_EPS)),
        logit(state.E),
        logit(state.A),
    )
end

function physical_state(state::TransformedState)
    bounded = bounded_transformed_state(state)
    return PhysicalState(
        exp(bounded.xC),
        exp(bounded.xD),
        exp(bounded.xG),
        logistic(bounded.xE),
        logistic(bounded.xA),
    )
end

physical_state(values::AbstractVector{<:Real}) = PhysicalState(values[1], values[2], values[3], values[4], values[5])

transformed_state(values::AbstractVector{<:Real}) = TransformedState(values[1], values[2], values[3], values[4], values[5])

state_vector(state::PhysicalState) = [state.C, state.D, state.G, state.E, state.A]

state_vector(state::TransformedState) = [state.xC, state.xD, state.xG, state.xE, state.xA]
