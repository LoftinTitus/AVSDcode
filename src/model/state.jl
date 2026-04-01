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

default_initial_state() = PhysicalState(0.22, 0.12, 1.00, 0.18, 0.24)

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
    return PhysicalState(
        exp(state.xC),
        exp(state.xD),
        exp(state.xG),
        logistic(state.xE),
        logistic(state.xA),
    )
end

physical_state(values::AbstractVector{<:Real}) = PhysicalState(values[1], values[2], values[3], values[4], values[5])

transformed_state(values::AbstractVector{<:Real}) = TransformedState(values[1], values[2], values[3], values[4], values[5])

state_vector(state::PhysicalState) = [state.C, state.D, state.G, state.E, state.A]

state_vector(state::TransformedState) = [state.xC, state.xD, state.xG, state.xE, state.xA]
