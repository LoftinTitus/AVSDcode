function prevalence(results::AbstractVector{<:EmbryoResult})
    isempty(results) && return 0.0
    return mean(result -> result.phenotype ? 1.0 : 0.0, results)
end

function summary_metrics(results::AbstractVector{<:EmbryoResult})
    isempty(results) && return Dict{Symbol,Float64}()

    gaps = [result.scores.gap for result in results]
    linear_scores = [result.scores.linear for result in results]
    nonlinear_scores = [result.scores.nonlinear for result in results]
    probabilities = [result.scores.probability for result in results]

    return Dict(
        :n_embryos => length(results),
        :prevalence => prevalence(results),
        :mean_gap => mean(gaps),
        :std_gap => length(gaps) > 1 ? std(gaps) : 0.0,
        :mean_linear_score => mean(linear_scores),
        :mean_nonlinear_score => mean(nonlinear_scores),
        :mean_probability => mean(probabilities),
    )
end
