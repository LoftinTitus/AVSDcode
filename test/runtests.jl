using Random
using Test

using AVSDModel

@testset "AVSDModel" begin
    @testset "State transforms are reversible" begin
        state = PhysicalState(0.35, 0.18, 0.82, 0.44, 0.63)
        recovered = physical_state(transformed_state(state))

        @test isapprox(recovered.C, state.C; atol = 1e-10)
        @test isapprox(recovered.D, state.D; atol = 1e-10)
        @test isapprox(recovered.G, state.G; atol = 1e-10)
        @test isapprox(recovered.E, state.E; atol = 1e-10)
        @test isapprox(recovered.A, state.A; atol = 1e-10)
        bounded = physical_state(TransformedState(Inf, -Inf, NaN, Inf, -Inf))
        @test all(isfinite, [bounded.C, bounded.D, bounded.G, bounded.E, bounded.A])
        @test 0.0 < bounded.E < 1.0
        @test 0.0 < bounded.A < 1.0
    end

    @testset "Supporting functions behave sensibly" begin
        @test 0.0 < logistic(0.0) < 1.0
        @test timing_window(1.0, 0.5, 1.5, 0.1) > timing_window(0.0, 0.5, 1.5, 0.1)
        @test late_penalty(1.0, 0.5, 1.0) < 1.0
        @test hill_activation(1.0, 3.0, 0.5) > hill_activation(0.2, 3.0, 0.5)
        @test excess_over_threshold(0.8, 1.0) == 0.0
    end

    @testset "Config loader returns a valid parameter set" begin
        params = load_parameters()
        targets = load_calibration_targets()
        run_config = load_research_run_config()
        production_config = load_research_run_config(joinpath(pwd(), "config", "research_run_production.yaml"))
        @test params.T == 3.5
        @test params.dt == 0.01
        @test haskey(params.rates, :alpha_EMT)
        @test params.genetics.latent_dim == 2
        @test targets.trisomy21_prevalence == (0.15, 0.33)
        @test run_config.n_chains == 4
        @test :threshold in run_config.parameter_names
        @test run_config.proposal_scale > 0.0
        @test run_config.refinement_candidates >= 1
        @test production_config.calibration_candidates > run_config.calibration_candidates
        @test production_config.burn_in > run_config.burn_in
    end

    @testset "Single embryo simulation stays in bounds" begin
        rng = MersenneTwister(7)
        params = default_parameters()
        result = solve_embryo(rng; params = params)
        final_state = result.trajectory.states[end]

        @test length(result.trajectory.times) > 10
        @test final_state.C > 0.0
        @test final_state.D > 0.0
        @test final_state.G > 0.0
        @test 0.0 < final_state.E < 1.0
        @test 0.0 < final_state.A < 1.0
        @test 0.0 <= result.scores.probability <= 1.0
        @test result.solver == :stochastic_heun
    end

    @testset "Stiff and higher-order solvers run" begin
        params = default_parameters()
        initial_state = default_initial_state()

        stiff_trajectory = bdf1_stiff(initial_state, params; dt = 0.02)
        srk_trajectory = stochastic_heun(MersenneTwister(11), initial_state, params; dt = 0.02)
        stiff_final = stiff_trajectory.states[end]
        srk_final = srk_trajectory.states[end]

        @test length(stiff_trajectory.times) == length(stiff_trajectory.states)
        @test length(srk_trajectory.times) == length(srk_trajectory.states)
        @test stiff_final.C > 0.0
        @test stiff_final.D > 0.0
        @test stiff_final.G > 0.0
        @test 0.0 < stiff_final.E < 1.0
        @test 0.0 < stiff_final.A < 1.0
        @test srk_final.C > 0.0
        @test srk_final.D > 0.0
        @test srk_final.G > 0.0
        @test 0.0 < srk_final.E < 1.0
        @test 0.0 < srk_final.A < 1.0
    end

    @testset "Solver selection flows through embryo and population simulations" begin
        params = default_parameters()
        rng = MersenneTwister(17)
        embryo = solve_embryo(rng; params = params, solver = :bdf1)
        population = simulate_population(6; params = params, seed = 123, solver = :euler_maruyama)

        @test embryo.solver == :bdf1
        @test all(result -> result.solver == :euler_maruyama, population)
    end

    @testset "Stratified population sampling improves cohort coverage" begin
        params = default_parameters()
        population = simulate_population(
            12;
            params = params,
            seed = 321,
            trisomy_fraction = 0.25,
            sample_strategy = :stratified,
        )

        @test count(result -> result.genotype.trisomy21, population) == 3
        @test all(length(result.genotype.latent_axes) == params.genetics.latent_dim for result in population)
        @test abs(sum(result.genotype.latent_axes[1] for result in population) / length(population)) < params.genetics.latent_sd
    end

    @testset "Population summaries and sensitivity run" begin
        params = default_parameters()
        results = simulate_population(12; params = params, seed = 21, trisomy_fraction = 0.25)
        metrics = summary_metrics(results)
        sens = local_sensitivity([:alpha_EMT, :k_G]; params = params, n_embryos = 10, seed = 99)

        @test length(results) == 12
        @test 0.0 <= prevalence(results) <= 1.0
        @test haskey(metrics, :mean_gap)
        @test haskey(metrics, :mean_closure_time)
        @test haskey(metrics, :mean_terminal_shear)
        @test haskey(sens, :alpha_EMT)
        @test haskey(sens, :k_G)
    end

    @testset "Lower structural scores map to AVSD phenotype risk" begin
        params = default_parameters()
        robust_state = PhysicalState(1.0, 0.8, 0.1, 0.9, 0.9)
        weak_state = PhysicalState(0.1, 0.1, 1.0, 0.1, 0.1)

        @test phenotype_scores(weak_state, params).linear < phenotype_scores(robust_state, params).linear
        @test phenotype_probability(weak_state, params) > phenotype_probability(robust_state, params)
        @test phenotype_label(weak_state, params)
        @test !phenotype_label(robust_state, params)
    end

    @testset "Calibration, fitting, validation, and global sensitivity run" begin
        params = default_parameters()
        targets = default_calibration_targets()
        panel = AVSDModel.build_calibration_panel(8; params = params, seed = 5)

        evaluation = calibration_summary(params; targets = targets, n_embryos = 8, seed = 5)
        fixed_panel_eval_a = calibration_summary(
            params;
            targets = targets,
            n_embryos = 8,
            seed = 5,
            calibration_panel = panel,
        )
        fixed_panel_eval_b = calibration_summary(
            params;
            targets = targets,
            n_embryos = 8,
            seed = 500,
            calibration_panel = panel,
        )
        fit = fit_parameters(
            params = params,
            targets = targets,
            parameter_names = [:alpha_EMT, :alpha_DMP, :k_G],
            n_candidates = 3,
            n_embryos = 8,
            seed = 12,
            refinement_rounds = 1,
            refinement_candidates = 2,
        )
        posterior = posterior_summary(fit; burn_in = 1)
        multichain = run_calibration_chains(
            params = params,
            targets = targets,
            parameter_names = [:alpha_EMT, :alpha_DMP, :k_G],
            n_candidates = 2,
            n_embryos = 8,
            n_chains = 2,
            seed = 1012,
        )
        multichain_posterior = posterior_summary(multichain; burn_in = 1)
        validation = validate_parameters(
            fit.best_parameters;
            targets = targets,
            n_replicates = 2,
            n_embryos = 8,
            seed = 13,
        )
        sweep = parameter_sweep(
            :alpha_EMT;
            params = params,
            targets = targets,
            factors = [0.9, 1.0, 1.1],
            n_embryos = 8,
            seed = 14,
        )
        gsens = global_sensitivity(
            [:alpha_EMT, :k_G];
            params = params,
            targets = targets,
            n_embryos = 8,
            sweep_factors = [0.9, 1.0, 1.1],
            seed = 15,
        )

        @test evaluation.score >= 0.0
        @test isfinite(evaluation.score)
        @test fixed_panel_eval_a.score == fixed_panel_eval_b.score
        @test haskey(evaluation.penalties, :trisomy21_prevalence)
        @test length(fit.history) == 6
        @test fit.best_evaluation.score >= 0.0
        @test isfinite(fit.best_evaluation.score)
        @test haskey(fit.history[1], :log_prior)
        @test haskey(fit.history[1], :log_posterior)
        @test haskey(fit.history[1], :accepted)
        @test haskey(fit.history[1], :parameter_values)
        @test 0.0 <= posterior.acceptance_rate <= 1.0
        @test posterior.n_samples == 5
        @test haskey(posterior.posterior_mean, :alpha_EMT)
        @test haskey(posterior.posterior_interval, :alpha_EMT)
        @test haskey(posterior.map_parameters, :k_G)
        @test length(multichain.chains) == 2
        @test haskey(multichain_posterior.rhat, :alpha_EMT)
        @test haskey(multichain_posterior.rhat, :score)
        @test length(validation.replicate_scores) == 2
        @test haskey(validation.pass_rates, :gap)
        @test length(sweep) == 3
        @test all(hasproperty(point, :score) for point in sweep)
        @test haskey(gsens.local_effects, :alpha_EMT)
        @test haskey(gsens.global_effects, :k_G)
        @test haskey(gsens.sweeps, :alpha_EMT)
    end

    @testset "Output and visualization tooling writes files" begin
        params = default_parameters()
        results = simulate_population(6; params = params, seed = 222)
        summary = summary_metrics(results)
        fit = fit_parameters(
            params = params,
            parameter_names = [:alpha_EMT, :k_G],
            n_candidates = 2,
            n_embryos = 6,
            seed = 223,
        )
        posterior = posterior_summary(fit; burn_in = 1)
        multichain = run_calibration_chains(
            params = params,
            parameter_names = [:alpha_EMT, :k_G],
            n_candidates = 2,
            n_embryos = 6,
            n_chains = 2,
            seed = 323,
        )
        validation = validate_parameters(fit.best_parameters; n_replicates = 2, n_embryos = 6, seed = 224)
        sensitivity = global_sensitivity(
            [:alpha_EMT, :k_G];
            params = fit.best_parameters,
            n_embryos = 6,
            sweep_factors = [0.9, 1.0, 1.1],
            seed = 225,
        )

        mktempdir() do dir
            pop_csv = write_population_csv(joinpath(dir, "population.csv"), results)
            sum_csv = write_summary_csv(joinpath(dir, "summary.csv"), summary; label = "baseline")
            hist_csv = write_calibration_history_csv(joinpath(dir, "history.csv"), fit)
            report_md = write_markdown_report(
                joinpath(dir, "report.md");
                summary = summary,
                calibration = fit.best_evaluation,
                posterior = posterior,
                validation = validation,
                sensitivity = sensitivity,
            )
            traj_svg = write_trajectory_svg(joinpath(dir, "gap.svg"), [result.trajectory for result in results]; variable = :G)
            hist_svg = write_calibration_history_svg(joinpath(dir, "history.svg"), fit)
            sens_svg = write_sensitivity_svg(joinpath(dir, "sensitivity.svg"), sensitivity)
            posterior_csv = write_posterior_summary_csv(joinpath(dir, "posterior.csv"), posterior)
            comparison_svg = write_cohort_comparison_svg(
                joinpath(dir, "cohorts.svg"),
                Dict(:euploid => summary, :mixed => summary, :trisomy21 => summary),
            )
            trace_svg = write_posterior_trace_svg(joinpath(dir, "trace.svg"), multichain; parameter_names = [:alpha_EMT, :k_G])
            bundle = generate_results_bundle(
                joinpath(dir, "bundle");
                params = params,
                initial_state = default_initial_state(),
                population_size = 6,
                calibration_candidates = 2,
                refinement_rounds = 1,
                refinement_candidates = 1,
                validation_replicates = 2,
                burn_in = 1,
                n_chains = 2,
                seed = 500,
            )
            cfg = ResearchRunConfig(
                joinpath(dir, "config_bundle"),
                joinpath(pwd(), "data", "calibration_targets.csv"),
                6,
                2,
                0.2,
                0.15,
                1,
                1,
                2,
                2,
                1,
                700,
                0.25,
                :stochastic_heun,
                :stratified,
                [:alpha_EMT, :k_G],
                [:alpha_EMT, :k_G],
                [:alpha_EMT, :k_G],
            )
            config_bundle = run_research_pipeline(cfg; params = params, initial_state = default_initial_state())

            @test isfile(pop_csv)
            @test isfile(sum_csv)
            @test isfile(hist_csv)
            @test isfile(posterior_csv)
            @test isfile(report_md)
            @test isfile(traj_svg)
            @test isfile(hist_svg)
            @test isfile(sens_svg)
            @test isfile(comparison_svg)
            @test isfile(trace_svg)
            @test occursin("embryo_id", read(pop_csv, String))
            @test occursin("baseline,prevalence", read(sum_csv, String))
            @test occursin("log_posterior", read(hist_csv, String))
            @test occursin("lower_credible", read(posterior_csv, String))
            @test occursin("MAP Parameters", format_posterior_summary(posterior))
            @test occursin("R-hat", format_posterior_summary(posterior_summary(multichain; burn_in = 1)))
            @test occursin("Posterior Summary", read(report_md, String))
            @test occursin("# AVSD Analysis Report", read(report_md, String))
            @test occursin("<svg", read(traj_svg, String))
            @test occursin("<svg", read(hist_svg, String))
            @test occursin("<svg", read(sens_svg, String))
            @test occursin("<svg", read(comparison_svg, String))
            @test occursin("<svg", read(trace_svg, String))
            @test haskey(bundle, :report_md)
            @test haskey(bundle, :gap_trajectories_svg)
            @test haskey(bundle, :cohort_comparison_svg)
            @test haskey(bundle, :posterior_trace_svg)
            @test isfile(bundle[:report_md])
            @test isfile(bundle[:posterior_summary_md])
            @test isfile(bundle[:run_manifest])
            @test occursin("sample_strategy=", read(bundle[:run_manifest], String))
            @test isfile(config_bundle[:report_md])
        end
    end
end
