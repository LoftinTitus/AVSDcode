function _time_grid(dt::T, T_final::T) where {T<:Real}
    n_steps = Int(floor(T_final / dt)) + 1
    return collect(range(zero(T), step = dt, length = n_steps))
end

function _finite_difference_jacobian(f, x::Vector{T}; epsilon::T = T(1.0e-6)) where {T<:Real}
    fx = f(x)
    n = length(x)
    jacobian = Matrix{T}(undef, n, n)

    for j in 1:n
        perturbation = zeros(T, n)
        perturbation[j] = epsilon
        jacobian[:, j] = (f(x .+ perturbation) .- fx) ./ epsilon
    end

    return jacobian, fx
end

function _transformed_state_vector(state::TransformedState{T}) where {T<:Real}
    return T[state.xC, state.xD, state.xG, state.xE, state.xA]
end

function _transformed_state_from_vector(values::AbstractVector{T}) where {T<:Real}
    return TransformedState(values[1], values[2], values[3], values[4], values[5])
end

function bdf1_stiff(
    initial_state::PhysicalState{T},
    params::ModelParameters{T};
    dt::T = params.dt,
    T_final::T = params.T,
    max_newton_iters::Int = 12,
    tolerance::T = T(1.0e-8),
) where {T<:Real}
    times = _time_grid(dt, T_final)
    states = Vector{PhysicalState{T}}(undef, length(times))
    states[1] = initial_state

    x_prev = transformed_state(initial_state)

    for i in 2:length(times)
        t_next = times[i]
        f_prev = transformed_drift(physical_state(x_prev), times[i - 1], params)
        x_guess = TransformedState(
            x_prev.xC + dt * f_prev.xC,
            x_prev.xD + dt * f_prev.xD,
            x_prev.xG + dt * f_prev.xG,
            x_prev.xE + dt * f_prev.xE,
            x_prev.xA + dt * f_prev.xA,
        )

        x_prev_vec = _transformed_state_vector(x_prev)
        x_vec = _transformed_state_vector(x_guess)

        function residual(candidate_vec::Vector{T})
            candidate_state = physical_state(_transformed_state_from_vector(candidate_vec))
            candidate_drift = transformed_drift(candidate_state, t_next, params)
            drift_vec = _transformed_state_vector(candidate_drift)
            return candidate_vec .- x_prev_vec .- dt .* drift_vec
        end

        converged = false

        for _ in 1:max_newton_iters
            jacobian, residual_value = _finite_difference_jacobian(residual, x_vec)

            if norm(residual_value, Inf) <= tolerance
                converged = true
                break
            end

            newton_step = jacobian \ residual_value
            x_vec .-= newton_step

            if norm(newton_step, Inf) <= tolerance
                converged = true
                break
            end
        end

        if !converged
            error("BDF1 stiff solver failed to converge at time $t_next")
        end

        x_prev = _transformed_state_from_vector(x_vec)
        states[i] = physical_state(x_prev)
    end

    return Trajectory(times, states)
end
