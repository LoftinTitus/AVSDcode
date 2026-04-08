using AVSDModel

const PRODUCTION_CONFIG_PATH = normpath(joinpath(@__DIR__, "..", "config", "research_run_production.yaml"))
const PRODUCTION_SEEDS = [2026, 3026, 4026]

base_config = load_research_run_config(PRODUCTION_CONFIG_PATH)

println("Running production research pipeline for seeds: ", join(PRODUCTION_SEEDS, ", "))

for seed in PRODUCTION_SEEDS
    output_dir = joinpath(dirname(base_config.output_dir), "production_seed$(seed)")
    config = ResearchRunConfig(
        output_dir,
        base_config.calibration_targets_path,
        base_config.population_size,
        base_config.calibration_candidates,
        base_config.proposal_scale,
        base_config.prior_mix,
        base_config.refinement_rounds,
        base_config.refinement_candidates,
        base_config.validation_replicates,
        base_config.n_chains,
        base_config.burn_in,
        seed,
        base_config.trisomy_fraction,
        base_config.solver,
        base_config.sample_strategy,
        base_config.parameter_names,
        base_config.sensitivity_parameters,
        base_config.trace_parameters,
    )

    paths = run_research_pipeline(config)
    println("Seed $(seed): ", paths[:report_md])
end
