function simulate_population(
    n_embryos::Integer;
    params::ModelParameters{T} = default_parameters(),
    initial_state::PhysicalState{T} = default_initial_state(),
    seed::Integer = 42,
    trisomy_fraction::Real = 0.0,
) where {T<:Real}
    rng = MersenneTwister(seed)
    results = Vector{EmbryoResult{T}}(undef, n_embryos)

    for idx in 1:n_embryos
        trisomy21 = rand(rng) < trisomy_fraction
        results[idx] = solve_embryo(
            rng;
            params = params,
            initial_state = initial_state,
            trisomy21 = trisomy21,
        )
    end

    return results
end
