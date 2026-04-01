function _safe_bounds(values::AbstractVector{<:Real})
    min_value = minimum(values)
    max_value = maximum(values)

    if min_value == max_value
        min_value -= 1.0
        max_value += 1.0
    end

    return float(min_value), float(max_value)
end

function _scale_x(value::Real, min_value::Real, max_value::Real, left::Real, width::Real)
    return left + width * (float(value) - min_value) / (max_value - min_value)
end

function _scale_y(value::Real, min_value::Real, max_value::Real, top::Real, height::Real)
    return top + height * (1 - (float(value) - min_value) / (max_value - min_value))
end

function _svg_header(width::Integer, height::Integer)
    return "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"$(width)\" height=\"$(height)\" viewBox=\"0 0 $(width) $(height)\">"
end

_svg_footer() = "</svg>"

function _state_accessor(variable::Symbol)
    if variable == :C
        return state -> state.C
    elseif variable == :D
        return state -> state.D
    elseif variable == :G
        return state -> state.G
    elseif variable == :E
        return state -> state.E
    elseif variable == :A
        return state -> state.A
    else
        error("Unsupported state variable $(variable)")
    end
end

function write_trajectory_svg(
    path::AbstractString,
    trajectories::AbstractVector{<:Trajectory};
    variable::Symbol = :G,
    width::Integer = 900,
    height::Integer = 420,
    title::AbstractString = "Trajectory Plot",
)
    isempty(trajectories) && error("At least one trajectory is required")

    accessor = _state_accessor(variable)
    all_times = vcat([trajectory.times for trajectory in trajectories]...)
    all_values = vcat([[accessor(state) for state in trajectory.states] for trajectory in trajectories]...)
    xmin, xmax = _safe_bounds(all_times)
    ymin, ymax = _safe_bounds(all_values)

    left = 70.0
    right = 20.0
    top = 45.0
    bottom = 45.0
    plot_width = width - left - right
    plot_height = height - top - bottom
    colors = ["#1f77b4", "#e45756", "#4c9f70", "#f3a712", "#6c5ce7", "#008b8b"]

    lines = String[
        _svg_header(width, height),
        "<rect x=\"0\" y=\"0\" width=\"$(width)\" height=\"$(height)\" fill=\"#ffffff\"/>",
        "<text x=\"$(width ÷ 2)\" y=\"26\" text-anchor=\"middle\" font-family=\"monospace\" font-size=\"18\">$(title)</text>",
        "<line x1=\"$(left)\" y1=\"$(top + plot_height)\" x2=\"$(left + plot_width)\" y2=\"$(top + plot_height)\" stroke=\"#444\" stroke-width=\"1.5\"/>",
        "<line x1=\"$(left)\" y1=\"$(top)\" x2=\"$(left)\" y2=\"$(top + plot_height)\" stroke=\"#444\" stroke-width=\"1.5\"/>",
        "<text x=\"$(width ÷ 2)\" y=\"$(height - 10)\" text-anchor=\"middle\" font-family=\"monospace\" font-size=\"12\">time</text>",
        "<text x=\"18\" y=\"$(height ÷ 2)\" text-anchor=\"middle\" transform=\"rotate(-90 18 $(height ÷ 2))\" font-family=\"monospace\" font-size=\"12\">$(variable)</text>",
    ]

    for (idx, trajectory) in enumerate(trajectories)
        points = String[]
        for (time, state) in zip(trajectory.times, trajectory.states)
            x = _scale_x(time, xmin, xmax, left, plot_width)
            y = _scale_y(accessor(state), ymin, ymax, top, plot_height)
            push!(points, _format_float(x; digits = 2) * "," * _format_float(y; digits = 2))
        end
        color = colors[mod1(idx, length(colors))]
        push!(lines, "<polyline fill=\"none\" stroke=\"$(color)\" stroke-width=\"1.5\" points=\"" * join(points, " ") * "\" opacity=\"0.8\"/>")
    end

    push!(lines, _svg_footer())

    open(path, "w") do io
        write(io, join(lines, "\n"))
        write(io, "\n")
    end

    return path
end

function write_calibration_history_svg(
    path::AbstractString,
    calibration::CalibrationResult;
    width::Integer = 900,
    height::Integer = 420,
    title::AbstractString = "Calibration Score History",
)
    iterations = [float(item[:iteration]) for item in calibration.history]
    scores = [float(item[:score]) for item in calibration.history]
    xmin, xmax = _safe_bounds(iterations)
    ymin, ymax = _safe_bounds(scores)

    left = 70.0
    right = 20.0
    top = 45.0
    bottom = 45.0
    plot_width = width - left - right
    plot_height = height - top - bottom

    points = String[]
    for (iteration, score) in zip(iterations, scores)
        x = _scale_x(iteration, xmin, xmax, left, plot_width)
        y = _scale_y(score, ymin, ymax, top, plot_height)
        push!(points, _format_float(x; digits = 2) * "," * _format_float(y; digits = 2))
    end

    lines = String[
        _svg_header(width, height),
        "<rect x=\"0\" y=\"0\" width=\"$(width)\" height=\"$(height)\" fill=\"#ffffff\"/>",
        "<text x=\"$(width ÷ 2)\" y=\"26\" text-anchor=\"middle\" font-family=\"monospace\" font-size=\"18\">$(title)</text>",
        "<line x1=\"$(left)\" y1=\"$(top + plot_height)\" x2=\"$(left + plot_width)\" y2=\"$(top + plot_height)\" stroke=\"#444\" stroke-width=\"1.5\"/>",
        "<line x1=\"$(left)\" y1=\"$(top)\" x2=\"$(left)\" y2=\"$(top + plot_height)\" stroke=\"#444\" stroke-width=\"1.5\"/>",
        "<polyline fill=\"none\" stroke=\"#1f77b4\" stroke-width=\"2\" points=\"" * join(points, " ") * "\"/>",
        _svg_footer(),
    ]

    open(path, "w") do io
        write(io, join(lines, "\n"))
        write(io, "\n")
    end

    return path
end

function write_sensitivity_svg(
    path::AbstractString,
    sensitivity::SensitivityResult;
    metric::Symbol = :score_span,
    width::Integer = 900,
    height::Integer = 420,
    title::AbstractString = "Sensitivity Ranking",
)
    items = sort(
        collect(sensitivity.global_effects);
        by = item -> get(last(item), metric, 0.0),
        rev = true,
    )
    isempty(items) && error("Sensitivity result has no global effects to plot")

    top = 40.0
    left = 180.0
    bar_height = 28.0
    gap = 16.0
    usable_width = width - left - 40.0
    max_value = maximum(get(values, metric, 0.0) for (_, values) in items)
    max_value = max(max_value, 1.0e-6)

    lines = String[
        _svg_header(width, height),
        "<rect x=\"0\" y=\"0\" width=\"$(width)\" height=\"$(height)\" fill=\"#ffffff\"/>",
        "<text x=\"$(width ÷ 2)\" y=\"24\" text-anchor=\"middle\" font-family=\"monospace\" font-size=\"18\">$(title)</text>",
    ]

    for (idx, (name, values)) in enumerate(items)
        value = get(values, metric, 0.0)
        y = top + (idx - 1) * (bar_height + gap)
        bar_width = usable_width * value / max_value
        push!(lines, "<text x=\"$(left - 10)\" y=\"$(y + 18)\" text-anchor=\"end\" font-family=\"monospace\" font-size=\"12\">$(name)</text>")
        push!(lines, "<rect x=\"$(left)\" y=\"$(y)\" width=\"$(bar_width)\" height=\"$(bar_height)\" fill=\"#4c9f70\"/>")
        push!(lines, "<text x=\"$(left + bar_width + 8)\" y=\"$(y + 18)\" font-family=\"monospace\" font-size=\"12\">$(_format_float(value; digits = 4))</text>")
    end

    push!(lines, _svg_footer())

    open(path, "w") do io
        write(io, join(lines, "\n"))
        write(io, "\n")
    end

    return path
end
