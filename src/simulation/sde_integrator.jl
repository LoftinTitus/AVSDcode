struct Trajectory{T<:Real}
    times::Vector{T}
    states::Vector{PhysicalState{T}}
end

function stochastic_heun(
    rng::AbstractRNG,
    initial_state::PhysicalState{T},
    params::ModelParameters{T};
    dt::T = params.dt,
    T_final::T = params.T,
) where {T<:Real}
    times = _time_grid(dt, T_final)
    states = Vector{PhysicalState{T}}(undef, length(times))
    states[1] = initial_state

    x = transformed_state(initial_state)
    sqrt_dt = sqrt(dt)

    for i in 2:length(times)
        t = times[i - 1]
        state_now = physical_state(x)
        drift_now = transformed_drift(state_now, t, params)
        diffusion_now = transformed_diffusion(state_now, t, params)
        noise = randn(rng, 5)
        dW = sqrt_dt .* noise

        predictor = TransformedState(
            x.xC + drift_now.xC * dt + diffusion_now.xC * dW[1],
            x.xD + drift_now.xD * dt + diffusion_now.xD * dW[2],
            x.xG + drift_now.xG * dt + diffusion_now.xG * dW[3],
            x.xE + drift_now.xE * dt + diffusion_now.xE * dW[4],
            x.xA + drift_now.xA * dt + diffusion_now.xA * dW[5],
        )

        state_predictor = physical_state(predictor)
        drift_predictor = transformed_drift(state_predictor, times[i], params)
        diffusion_predictor = transformed_diffusion(state_predictor, times[i], params)

        x = TransformedState(
            x.xC + 0.5 * (drift_now.xC + drift_predictor.xC) * dt + 0.5 * (diffusion_now.xC + diffusion_predictor.xC) * dW[1],
            x.xD + 0.5 * (drift_now.xD + drift_predictor.xD) * dt + 0.5 * (diffusion_now.xD + diffusion_predictor.xD) * dW[2],
            x.xG + 0.5 * (drift_now.xG + drift_predictor.xG) * dt + 0.5 * (diffusion_now.xG + diffusion_predictor.xG) * dW[3],
            x.xE + 0.5 * (drift_now.xE + drift_predictor.xE) * dt + 0.5 * (diffusion_now.xE + diffusion_predictor.xE) * dW[4],
            x.xA + 0.5 * (drift_now.xA + drift_predictor.xA) * dt + 0.5 * (diffusion_now.xA + diffusion_predictor.xA) * dW[5],
        )

        states[i] = physical_state(x)
    end

    return Trajectory(times, states)
end

function euler_maruyama(
    rng::AbstractRNG,
    initial_state::PhysicalState{T},
    params::ModelParameters{T};
    dt::T = params.dt,
    T_final::T = params.T,
) where {T<:Real}
    times = _time_grid(dt, T_final)
    states = Vector{PhysicalState{T}}(undef, length(times))
    states[1] = initial_state

    x = transformed_state(initial_state)
    sqrt_dt = sqrt(dt)

    for i in 2:length(times)
        t = times[i - 1]
        current_state = physical_state(x)
        mu = transformed_drift(current_state, t, params)
        sigma = transformed_diffusion(current_state, t, params)
        noise = randn(rng, 5)

        x = TransformedState(
            x.xC + mu.xC * dt + sigma.xC * sqrt_dt * noise[1],
            x.xD + mu.xD * dt + sigma.xD * sqrt_dt * noise[2],
            x.xG + mu.xG * dt + sigma.xG * sqrt_dt * noise[3],
            x.xE + mu.xE * dt + sigma.xE * sqrt_dt * noise[4],
            x.xA + mu.xA * dt + sigma.xA * sqrt_dt * noise[5],
        )

        states[i] = physical_state(x)
    end

    return Trajectory(times, states)
end
