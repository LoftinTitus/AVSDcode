function trajectory_closes(trajectory::Trajectory; gap_threshold::Real = 0.15)
    return any(state -> state.G <= gap_threshold, trajectory.states)
end

function closure_time(trajectory::Trajectory; gap_threshold::Real = 0.15)
    for (time, state) in zip(trajectory.times, trajectory.states)
        state.G <= gap_threshold && return time
    end
    return trajectory.times[end]
end

function terminal_shear(result::EmbryoResult)
    return shear_stress(result.trajectory.times[end], result.trajectory.states[end], result.parameters)
end

function peak_shear(result::EmbryoResult)
    return maximum(
        shear_stress(time, state, result.parameters) for (time, state) in zip(result.trajectory.times, result.trajectory.states)
    )
end

function prevalence(results::AbstractVector{<:EmbryoResult})
    isempty(results) && return 0.0
    return mean(result -> result.phenotype ? 1.0 : 0.0, results)
end

function summary_metrics(results::AbstractVector{<:EmbryoResult}; gap_threshold::Real = 0.15)
    isempty(results) && return Dict{Symbol,Float64}()

    gaps = [result.scores.gap for result in results]
    linear_scores = [result.scores.linear for result in results]
    nonlinear_scores = [result.scores.nonlinear for result in results]
    probabilities = [result.scores.probability for result in results]
    closure_times = [closure_time(result.trajectory; gap_threshold = gap_threshold) for result in results]
    closure_fraction = mean(result -> trajectory_closes(result.trajectory; gap_threshold = gap_threshold) ? 1.0 : 0.0, results)
    terminal_shears = [terminal_shear(result) for result in results]
    peak_shears = [peak_shear(result) for result in results]

    return Dict{Symbol,Float64}(
        :n_embryos => Float64(length(results)),
        :prevalence => Float64(prevalence(results)),
        :mean_gap => Float64(mean(gaps)),
        :std_gap => Float64(length(gaps) > 1 ? std(gaps) : 0.0),
        :mean_linear_score => Float64(mean(linear_scores)),
        :mean_nonlinear_score => Float64(mean(nonlinear_scores)),
        :mean_probability => Float64(mean(probabilities)),
        :mean_closure_time => Float64(mean(closure_times)),
        :closure_fraction => Float64(closure_fraction),
        :mean_terminal_shear => Float64(mean(terminal_shears)),
        :mean_peak_shear => Float64(mean(peak_shears)),
    )
end
