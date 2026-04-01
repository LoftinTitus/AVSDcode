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
    end

    @testset "Population summaries and sensitivity run" begin
        params = default_parameters()
        results = simulate_population(12; params = params, seed = 21, trisomy_fraction = 0.25)
        metrics = summary_metrics(results)
        sens = local_sensitivity([:alpha_EMT, :k_G]; params = params, n_embryos = 10, seed = 99)

        @test length(results) == 12
        @test 0.0 <= prevalence(results) <= 1.0
        @test haskey(metrics, :mean_gap)
        @test haskey(sens, :alpha_EMT)
        @test haskey(sens, :k_G)
    end
end
