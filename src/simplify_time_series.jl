# =============================================
# Algorithm that tries to fit a line between
# start and end point and splits at the point
# of largest distance (given a custom function)
# =============================================

function simplify_ts_dist(ts::TimeSeries, dist, δ::Number=SPLIT_DISTANCE_ALT)
    # linear interpolation of first and last point
    x0 = ts.data[1]
    x1 = ts.data[end]
    t0 = ts.time[1]
    t1 = ts.time[end]
    Δt = t1 .- t0
    slope = (x1 .- x0) ./ Δt

    # compute distance to all data points
    dists = Vector{Float64}(undef, length(ts) - 2)
    @inbounds for i in 2:length(ts)-1
        Δt = ts.time[i] .- t0
        y = x0 .+ slope .* Δt
        dists[i-1] = dist(y, ts.data[i])
    end

    # check if the interpolation is sufficiently precise
    if !isempty(dists)
        max_dist = maximum(dists)
    end
    if isempty(dists) || max_dist <= δ
        time = [t0, t1]
        data = [x0, x1]
        switches = Int[]
        return switches, TimeSeries(time, data)
    end

    # find first point with largest distance
    j = findfirst(dists .== max_dist) + 1

    # recursively split in two
    ts1 = TimeSeries(ts.time[1:j], ts.data[1:j])
    switches1, ts1_simp = simplify_ts_dist(ts1, dist, δ)

    ts2 = TimeSeries(ts.time[j:end], ts.data[j:end])
    switches2, ts2_simp = simplify_ts_dist(ts2, dist, δ)

    # combine results
    switches = vcat(switches1, j, switches2 .+ (j - 1))
    time_simplified = vcat(ts1_simp.time, ts2_simp.time[2:end])
    data_simplified = vcat(ts1_simp.data, ts2_simp.data[2:end])
    ts_simplified = TimeSeries(time_simplified, data_simplified)

    return switches, ts_simplified
end

# ================================================
# Functions that take a time series and return the
# change points and the corresponding slopes
# ================================================

# identity
function cpds_id(ts)
    switches = 2:length(ts)-1
    slopes = get_slopes(ts)
    return switches, slopes
end

# RDP algorithm
function cpds_rdp(ts::TimeSeries, δ::Number=SPLIT_DISTANCE_RDP)
    switches, ts_simplified = simplify_ts_rdp(ts, δ)
    slopes = get_slopes(ts_simplified)
    return switches, slopes
end

# dist algorithm
function cpds_dist(ts::TimeSeries, dist, δ::Number=SPLIT_DISTANCE_ALT)
    switches, ts_simplified = simplify_ts_dist(ts, dist, δ)
    slopes = get_slopes(ts_simplified)
    return switches, slopes
end

# all pieces
function cpds_arg_all_pieces(ts::TimeSeries)
    switches = 2:length(ts)-1
    slopes = get_slopes(ts)
    return switches, slopes
end

# =================================================
# Functions that return another function which take
# only a time series and return the change points
# and the corresponding slopes
# =================================================

# RDP algorithm
function cpds_arg_rdp(δ::Number=SPLIT_DISTANCE_RDP)
    return ts -> cpds_rdp(ts, δ)
end

# dist algorithm
function cpds_arg_dist(dist, δ::Number=SPLIT_DISTANCE_ALT)
    return ts -> cpds_dist(ts, dist, δ)
end

function cpds_arg_all_pieces()
    return ts -> cpds_arg_all_pieces(ts)
end

# ======================================================
# Function to remove points in a time series that are
# close to the linear interpolation of its two neighbors
# ======================================================

function simplify_ts_linear_interpolation!(ts::TimeSeries; δ=0.005, debug=0)
    i = 1
    n_old = length(ts)
    n = n_old - 2
    @inbounds while i <= n
        x = ts.data[i]
        y = ts.data[i+1]
        z = ts.data[i+2]
        tx = ts.time[i]
        ty = ts.time[i+1]
        slope = z - x
        Δt = ty - tx
        y_li = x + slope * Δt
        Δy = relative_distance(y, y_li)
        if Δy < δ
            deleteat!(ts.data, i+1)
            deleteat!(ts.time, i+1)
            n -= 1
        else
            i += 1
        end
    end
    n_new = length(ts)
    Δts = n_old - n_new
    debug > 0 && println("removed $Δts / $n_old points in time series #$debug")
    return ts
end

function relative_distance(x, y)
    Δ = abs.(x - y)
    d = -Inf
    for i in 1:length(Δ)
        if x[i] != 0
            di = abs(Δ[i] / x[i])
        elseif y[i] != 0
            di = abs(Δ[i] / y[i])
        else
            di = zero(Δ)
        end
        d = max(d, di)
    end
    return d
end
