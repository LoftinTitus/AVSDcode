const STATE_EPS = 1.0e-6

logistic(x::Real) = inv(1 + exp(-x))

function logit(x::Real; eps::Float64 = STATE_EPS)
    clipped = clamp(float(x), eps, 1 - eps)
    return log(clipped / (1 - clipped))
end

function timing_window(
    t::Real,
    t_on::Real,
    t_off::Real,
    smoothness::Real,
)
    return logistic((t - t_on) / smoothness) * logistic((t_off - t) / smoothness)
end

late_penalty(t::Real, t_dead::Real, lambda::Real) = exp(-lambda * max(t - t_dead, 0.0))

function hill_activation(x::Real, n::Real, K::Real)
    x_pos = max(float(x), 0.0)
    numerator = x_pos^n
    denominator = numerator + K^n
    return iszero(denominator) ? 0.0 : numerator / denominator
end

excess_over_threshold(x::Real, threshold::Real) = max(float(x) - float(threshold), 0.0)
