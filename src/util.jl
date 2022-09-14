import DelimitedFiles, Random
using DelimitedFiles: readdlm
using Plots
import ReachabilityAnalysis: solve

# time series
struct TimeSeries
    time::Vector{Float64}  # vector of time points
    data::Vector{Vector{Float64}}  # vector of n-dimensional data points

    # constructor
    function TimeSeries(time, data)
        @assert length(time) > 1
        @assert length(time) == length(data)
        if data isa Vector{Float64}  # auto-convert 1D time series
            data = [[xi] for xi in data]
        end
        return new(time, data)
    end
end

struct Dataset
    traj::Vector{Int}  # vector of trajectory indices
    time::Vector{Float64}  # vector of time points
    data::Vector{Vector{Float64}}  # vector of n-dimensional data points
end

# constructor from a dataset with a trajectory column and a trajectory number
function TimeSeries(d::Dataset, trajno::Int)
    time_trajectory = d.time[d.traj[:] .== trajno]
    data_trajectory = d.data[d.traj[:] .== trajno]
    return TimeSeries(time_trajectory, data_trajectory)
end

function LazySets.dim(ts::TimeSeries)
    return length(ts.data[1])
end

function Base.length(ts::TimeSeries)
    return length(ts.time)
end

function read_dataset(path, c_traj, c_time, cs_data; has_header)
    dataset_raw = readdlm(path, ',')
    if has_header
        traj = dataset_raw[2:end, c_traj] * 1
        time = dataset_raw[2:end, c_time] * 1.0
        data = dataset_raw[2:end, cs_data] * 1.0
    else
        traj = dataset_raw[:, c_traj] * 1
        time = dataset_raw[:, c_time] * 1.0
        data = dataset_raw[:, cs_data] * 1.0
    end

    # convert matrix to vector of data points
    data = [data[i, :] for i in 1:size(data, 1)]

    return Dataset(traj, time, data)
end

abstract type EpsilonMode end

struct SingleEpsilonMode <: EpsilonMode end

mutable struct MultiEpsilonMode <: EpsilonMode end

function plot_ts!(ts::TimeSeries; dim=1, c=:black, alpha=1, m=:+,
                  ylab="x$dim", ms=5)
    n = length(ts.time)
    x = Vector{Float64}(undef, n)
    y = Vector{Float64}(undef, n)
    for i in 1:n
        x[i] = ts.time[i]
        y[i] = ts.data[i][dim]
    end

    scatter!(x, y, seriestype=:shape, alpha=alpha, c=c, lab="", xlab="t",
             m=m, ylab=ylab, ms=ms)
end

# plotting a tube around a time series
function plot_tube!(ts::TimeSeries, ε; dim=1, c=:green, alpha=0.5, ylab="x$dim")
    m = length(ts.time)
    n = 2 * m
    x = Vector{Float64}(undef, n)
    y = Vector{Float64}(undef, n)
    for i in 1:m
        x[i] = ts.time[i]
        x[n - i + 1] = ts.time[i]
        y[i] = ts.data[i][dim] + ε
        y[n - i + 1] = ts.data[i][dim] - ε
    end

    plot!(x, y, seriestype=:shape, alpha=alpha, c=c, lab="", xlab="t",
          ylab=ylab)
end

# infinity-norm ball
ball(x, ε) = BallInf(x, ε)

# solution of a constant-dynamics system x'(t) = m with x(0) = x0 at time point t
function solve(m::AbstractVector{Float64}, x0, t)
    return x0 .+ m * t
end

function get_slopes(ts::TimeSeries)
    m = length(ts.data)
    return @inbounds [(ts.data[i] .- ts.data[i-1]) ./ (ts.time[i] .- ts.time[i-1])
                      for i in 2:m]
end

function get_slopes(all_ts::AbstractVector{TimeSeries})
    vcat([get_slopes(ts) for ts in all_ts]...)
end

function get_delays(ts::TimeSeries, switches)
    delays = zeros(length(switches) + 1)
    t0 = ts.time[1]
    @inbounds for (i, si) in enumerate(switches)
        t1 = ts.time[si]
        delays[i] = t1 - t0
        t0 = t1
    end
    @inbounds delays[end] = ts.time[end] - t0
    @assert sum(delays) ≈ ts.time[end] - ts.time[1]
    return delays
end

# overwrite a Lazysets implementation as a workaround
function LazySets.vertices_list(P::HPolytope{N};
                                backend=nothing, prune::Bool=true) where {N}
    if length(P.constraints) == 0
        return Vector{N}(Vector{N}(undef, 0))
    end

    if dim(P) == 2 && backend == nothing
        return vertices_list(convert(HPolygon, P, prune=prune))
    else
        if backend == nothing
            backend = LazySets.default_polyhedra_backend(P)
        end
        Q = polyhedron(P; backend=backend)
        return collect(Polyhedra.points(Q))
    end
end

function reset_seed(seed)
    Random.seed!(seed)
end

function extract_mode_sequences(switchesv, mode_sequence_all)
    mode_sequencev = Vector{Int}[]
    i = 1
    for switches in switchesv
        j = i + length(switches)
        mode_sequence = mode_sequence_all[i:j]
        push!(mode_sequencev, mode_sequence)
        i = j + 1
    end
    @assert i - 1 == length(mode_sequence_all) "inconsistent vector length"
    return mode_sequencev
end

function extract_vectors(concatenated, n)
    vectors = Vector{Float64}[]
    i = 1
    while i <= length(concatenated)
        j = i + n - 1
        v = concatenated[i:j]
        push!(vectors, v)
        i = j + 1
    end
    @assert i - 1 == length(concatenated) "inconsistent vector length"
    return vectors
end

function plot_results_together(tsv, tsv2, ε; d::Int=1, ylab="none")
    fig = plot()
    [plot_tube!(ts1, ε, c=:green, dim=d, ylab=(ylab == "none" ? "x$d" : ylab)) for ts1 in tsv]
    [plot_tube!(ts2, TUBE_RADIUS, c=:red, dim=d, ylab=(ylab == "none" ? "x$d" : ylab)) for ts2 in tsv2]
    return fig
end

function plot_results_individual(tsv, tsv2, ε; d::Int=1, ylab="none", fig=nothing)
    figs = []
    for (ts1, ts2) in zip(tsv, tsv2)
        if isnothing(fig)
            fig = plot()
        end
        plot_tube!(ts1, ε, c=:green, dim=d, ylab=(ylab == "none" ? "x$d" : ylab))
        plot_tube!(ts2, TUBE_RADIUS, c=:red, dim=d, ylab=(ylab == "none" ? "x$d" : ylab))
        push!(figs, fig)
    end
    fig = plot(figs...)
    return fig
end

function plot_results(tsv, tsv2, ε; d::Int=1)
    fig1 = plot_results_together(tsv, tsv2, ε, d=d)
    fig2 = plot_results_individual(tsv, tsv2, ε, d=d)
    fig = plot(fig1, fig2, layout = (2, 1))
    return fig
end
