struct Trajectory{T<:Real}
    times::Vector{T}
    states::Vector{PhysicalState{T}}
end

function euler_maruyama(
    rng::AbstractRNG,
    initial_state::PhysicalState{T},
    params::ModelParameters{T};
    dt::T = params.dt,
    T_final::T = params.T,
) where {T<:Real}
    n_steps = Int(floor(T_final / dt)) + 1
    times = collect(range(zero(T), step = dt, length = n_steps))
    states = Vector{PhysicalState{T}}(undef, n_steps)
    states[1] = initial_state

    x = transformed_state(initial_state)
    sqrt_dt = sqrt(dt)

    for i in 2:n_steps
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
