using AVSDModel

params = load_parameters()
results = simulate_population(24; params = params, seed = 2027, trisomy_fraction = 0.25)
summary = summary_metrics(results)

println(format_summary_table(summary))

write_population_csv("population_outputs.csv", results)
write_summary_csv("population_summary.csv", summary; label = "mixed_cohort")
write_trajectory_svg("gap_trajectories.svg", [result.trajectory for result in results]; variable = :G, title = "Gap Trajectories")

println("Wrote CSV and SVG outputs to the current directory.")
