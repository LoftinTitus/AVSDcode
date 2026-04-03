using AVSDModel

params = load_parameters()
initial_state = load_initial_state()
targets = default_calibration_targets()

paths = generate_results_bundle(
    joinpath(@__DIR__, "generated_results");
    params = params,
    targets = targets,
    initial_state = initial_state,
    population_size = 48,
    calibration_candidates = 12,
    validation_replicates = 3,
    burn_in = 3,
    seed = 2026,
)

println("Generated research outputs:")
for key in sort(collect(keys(paths)); by = string)
    println(" - ", key, ": ", paths[key])
end
