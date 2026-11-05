"""
    logsumexp(xs::AbstractArray; α = 1.0)

Return the smoothed maximum of an array, where `α` controls the steepness of the smoothing.
"""
logsumexp(xs::AbstractArray; α = 1.0) = 1 / α * log(sum([exp(α * x) for x in xs]))
"""
    logsumexp(x; α = 1.0)

Return the smoothed maximum of the argument and zero.
"""
logsumexp(x; α = 1.0) = logsumexp([x, zero(x)]; α)

"""
    smooth_daynight(t, t_sunrise, t_sunset, ymin = 0, ymax = 1)

Return a value that smoothly switches between `ymax` at day and `ymin` at night.
"""
function smooth_daynight(t, t_sunrise, t_sunset, ymin = 0, ymax = 1; smoothing = 0.01)
    return daywave(t, t_sunrise, t_sunset) |>
        x -> smooth_sign(x, smoothing) |>
        x -> rebound_sign(x, ymin, ymax)
end

get_offset(d) = -sin((1 - 2 * d) / 2 * pi) # get value of `a` so that the roots of sin(x) + a = 0 have a maximal distance of `d`
offset_sin(x, d) = sin(x + asin(-get_offset(d))) + get_offset(d) # offset a sine a vertical height of `d`, and offset it horizontally so f(0) = 0
daywave(x, t_sunrise, t_sunset) = offset_sin(2 * pi * (x - t_sunrise) / 24, (t_sunset - t_sunrise) / 24) # get an offset sine wave in hours with values >0 for the day
smooth_sign(x, smoothing) = x / (sqrt(x^2 + smoothing^2)) # smooth approximation of the `sign` function
rebound_sign(x, ymin, ymax) = 1 / 2 * ((ymax + ymin) + (ymax - ymin) * x) # for a function with range [-1, 1], change its range to [ymin, ymax]
