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
        @test params.T == 3.5
        @test params.dt == 0.01
        @test haskey(params.rates, :alpha_EMT)
        @test params.genetics.latent_dim == 2
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

    @testset "Calibration, fitting, validation, and global sensitivity run" begin
        params = default_parameters()
        targets = default_calibration_targets()

        evaluation = calibration_summary(params; targets = targets, n_embryos = 8, seed = 5)
        fit = fit_parameters(
            params = params,
            targets = targets,
            parameter_names = [:alpha_EMT, :alpha_DMP, :k_G],
            n_candidates = 3,
            n_embryos = 8,
            seed = 12,
        )
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
        @test haskey(evaluation.penalties, :trisomy21_prevalence)
        @test length(fit.history) == 4
        @test fit.best_evaluation.score >= 0.0
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
                validation = validation,
                sensitivity = sensitivity,
            )
            traj_svg = write_trajectory_svg(joinpath(dir, "gap.svg"), [result.trajectory for result in results]; variable = :G)
            hist_svg = write_calibration_history_svg(joinpath(dir, "history.svg"), fit)
            sens_svg = write_sensitivity_svg(joinpath(dir, "sensitivity.svg"), sensitivity)

            @test isfile(pop_csv)
            @test isfile(sum_csv)
            @test isfile(hist_csv)
            @test isfile(report_md)
            @test isfile(traj_svg)
            @test isfile(hist_svg)
            @test isfile(sens_svg)
            @test occursin("embryo_id", read(pop_csv, String))
            @test occursin("baseline,prevalence", read(sum_csv, String))
            @test occursin("# AVSD Analysis Report", read(report_md, String))
            @test occursin("<svg", read(traj_svg, String))
            @test occursin("<svg", read(hist_svg, String))
            @test occursin("<svg", read(sens_svg, String))
        end
    end
end
