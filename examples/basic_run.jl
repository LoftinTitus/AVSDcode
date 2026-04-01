using AVSDModel

params = load_parameters()
initial_state = load_initial_state()
results = simulate_population(100; params = params, initial_state = initial_state, seed = 2026, trisomy_fraction = 0.2)

println("Population summary:")
println(summary_metrics(results))
