using ReachabilityAnalysis, Symbolics, DifferentialEquations
using ReachabilityAnalysis: HalfSpace
import DataStructures
using Plots
using Plots.Measures
using JLD2, FileIO  # for storing results

import DataStructures: OrderedDict

include("../models/thermostat/multi_thermostat.jl")

ns = 1:2:7  # data dimension
rs = vcat(1, 20:20:60)  # number of time series
points = 50:50:200  # number of points of each time series
λs = vcat(1, 5:5:15)  # number of locations

function multi_thermostat_X0(n)
    X0 = Hyperrectangle(low=fill(68.0, n), high=fill(69.0, n))
    return [(1, X0)]
end

function run_multi_thermostat(; n=2,
                             r=1,
                             T=10.0,
                             seed=0,
                             change_points_and_slopes=cpds_arg_dist((x, y) -> norm(x - y), 1.0),
                             ε_mode=SingleEpsilonMode(),
                             debug::Bool=false)
    !isnothing(seed) && reset_seed(seed)
    H = multi_thermostat(n)
    X0 = multi_thermostat_X0(n)
    prob = IVP(H, X0)

    debug && println("obtaining $r simulations from the original system...")
    @time tsv = simulate(prob, T, r)
    debug && println("...done")

    H, ε, x0v, switchesv, mode_sequencev =
        analyze(tsv; change_points_and_slopes=change_points_and_slopes,
                ε_mode=ε_mode, debug=debug)

    return H, tsv, ε, x0v, switchesv, mode_sequencev
end

function run_multi_thermostat_all(; n=2,
                                  r=1,
                                  T=10.0,
                                  seed=0,
                                  plot_dim::Int=0,
                                  ε_mode=SingleEpsilonMode(),
                                  debug::Bool=true)
    H, tsv, ε, x0v, switchesv, mode_sequencev =
        run_multi_thermostat(T=T, r=r, seed=seed, ε_mode=ε_mode, debug=debug,
                             n=n)

    tsv2 = simulate_imitate(H, tsv, x0v, switchesv, mode_sequencev)

    if plot_dim != 0
        fig = plot_results(tsv, tsv2, ε, d=plot_dim)
        return H, ε, tsv, tsv2, fig
    end

    return H, ε, tsv, tsv2
end

function benchmark_multi_thermostat_one(; seed, n, r, λ, tsv, ε_mode)
    reset_seed(seed)
    kmeans_args = (λ, SqEuclidean(), CLUSTERING_THRESHOLD)

    println("running synthesis from $r time series in $n dimensions...")
    res = @timed analyze(tsv;
                         change_points_and_slopes=cpds_arg_all_pieces,
                         kmeans_args=kmeans_args, ε_mode=ε_mode, debug=false)
    H, ε, x0v, switchesv, mode_sequencev = res.value
    runtime = res.time
    println("...done")

    return H, ε, x0v, switchesv, mode_sequencev, runtime
end

# 52 experiments
function benchmark_multi_thermostat!(results=Dict(); ε_mode=SingleEpsilonMode(),
                                     seed=0, warmup::Bool=false)
    T = 40.0

    i = 0
    filter_dict = get_filter_dict()
    for n in ns
        H = multi_thermostat(n)
        X0 = multi_thermostat_X0(n)
        prob = IVP(H, X0)
        for r in rs
            for n_points in points
                reset_seed(seed)
                tsv = simulate(prob, T, r; n_points=n_points)

                for λ in λs
                    tuple = (n, r, n_points, λ)
                    if tuple ∉ filter_dict
                        continue
                    end
                    λ = max(λ, 1)
                    λ >= n_points && break  # does not happen with the given data
                    i += 1
                    H, ε, x0v, switchesv, mode_sequencev, runtime =
                        benchmark_multi_thermostat_one(seed=seed, n=n, r=r, λ=λ,
                                                       tsv=tsv, ε_mode=ε_mode)
                    if warmup
                        return  # only run algorithm once in `warmup` mode
                    end
                    println("$i: n = $n, n_points = $n_points, r = $r, λ = $λ: $runtime seconds")
                    results[tuple] = runtime
                end
            end
        end
    end

    # store results for later
    file = File{format"JLD2"}("scalability_results.jld2")
    save(file, "results", results)
    # load via:
    # data = load(file)
    # results = data["results"]

    return sort!(OrderedDict(results))
end

function filter_data(results, j1, j2, j3)
    dims = [ns, rs, points, λs]
    return [filter(x -> all(
            x[1][j1] == dims[j1][i] &&
            x[1][j2] == dims[j2][i] &&
            x[1][j3] == dims[j3][i]),
        pairs(results)) for i in 1:4]
end

function plot_label(i, j1, j2, j3)
    function lab(idx)
        idx == 1 && return "n"
        idx == 2 && return "r"
        idx == 3 && return "p"
        idx == 4 && return "\\lambda"
        error("invalid input $idx")
    end

    function val(idx, i)
        idx == 1 && return "$(ns[i])"
        idx == 2 && return "$(rs[i])"
        idx == 3 && return "$(points[i])"
        idx == 4 && return "$(λs[i])"
        error("invalid input $idx")
    end

    l1 = lab(j1)
    l2 = lab(j2)
    l3 = lab(j3)
    v1 = val(j1, i)
    v2 = val(j2, i)
    v3 = val(j3, i)
    return "$l1 = $v1, $l2 = $v2, $l3 = $v3"
end

function plot_scalability_results(results)
    j1, j2, j3 = 2, 3, 4
    n_filtered = filter_data(results, j1, j2, j3)
    fig_n = plot(size=(900, 400), tickfont=font(20, "Times"), guidefontsize=20,
                 legendfont=font(15, "Times"), bottom_margin=8mm, left_margin=6mm,
                 xlab="dimension (n)", ylab="time [sec]", xticks=collect(ns), leg=:topleft)
    [plot!(fig_n, ns, collect(values(sort!(OrderedDict(p)))), m=:rect,
           lab=plot_label(i, j1, j2, j3))
        for (i, p) in enumerate(n_filtered)]
    savefig(fig_n, "scalability_n.pdf")

    j1, j2, j3 = 1, 3, 4
    r_filtered = filter_data(results, j1, j2, j3)
    fig_r = plot(size=(900, 400), tickfont=font(20, "Times"), guidefontsize=20,
                 legendfont=font(15, "Times"), bottom_margin=8mm, left_margin=6mm,
                 xlab="number of time series (r)", ylab="time [sec]", xticks=collect(rs), leg=:topleft)
    [plot!(fig_r, rs, collect(values(sort!(OrderedDict(p)))), m=:circle,
           lab=plot_label(i, j1, j2, j3))
        for (i, p) in enumerate(r_filtered)]
    savefig(fig_r, "scalability_r.pdf")

    j1, j2, j3 = 1, 2, 4
    p_filtered = filter_data(results, j1, j2, j3)
    fig_p = plot(size=(900, 400), tickfont=font(20, "Times"), guidefontsize=20,
                 legendfont=font(15, "Times"), bottom_margin=8mm, left_margin=6mm,
                 xlab="number of points in time series (p)", ylab="time [sec]", xticks=collect(points), leg=:topleft)
    [plot!(fig_p, points, collect(values(sort!(OrderedDict(p)))), m=:diamond,
           lab=plot_label(i, j1, j2, j3))
        for (i, p) in enumerate(p_filtered)]
    savefig(fig_p, "scalability_p.pdf")

    j1, j2, j3 = 1, 2, 3
    λ_filtered = filter_data(results, j1, j2, j3)
    fig_λ = plot(size=(900, 400), tickfont=font(20, "Times"), guidefontsize=20,
                 legendfont=font(15, "Times"), bottom_margin=8mm, left_margin=6mm,
                 xlab="number of locations (λ)", ylab="time [sec]", xticks=collect(λs), leg=:right)
    [plot!(fig_λ, λs, collect(values(sort!(OrderedDict(p)))), m=:star5,
           lab=plot_label(i, j1, j2, j3))
        for (i, p) in enumerate(λ_filtered)]
    savefig(fig_λ, "scalability_lambda.pdf")
end

function get_filter_dict()
    results = Dict()
    for n in ns
        for r in rs
            for n_points in points
                for λ in λs
                    results[(n, r, n_points, λ)] = 1
                end
            end
        end
    end

    T = Tuple{Int, Int, Int, Int}
    j1, j2, j3 = 2, 3, 4
    n_filtered = filter_data(results, j1, j2, j3)
    tuples = Set(reduce(append!, keys.(n_filtered[i] for i in 1:length(n_filtered)), init=T[]))
    j1, j2, j3 = 1, 3, 4
    r_filtered = filter_data(results, j1, j2, j3)
    tuples = union(tuples, Set(reduce(append!, keys.(r_filtered[i] for i in 1:length(r_filtered)), init=T[])))
    j1, j2, j3 = 1, 2, 4
    p_filtered = filter_data(results, j1, j2, j3)
    tuples = union(tuples, Set(reduce(append!, keys.(p_filtered[i] for i in 1:length(p_filtered)), init=T[])))
    j1, j2, j3 = 1, 2, 3
    λ_filtered = filter_data(results, j1, j2, j3)
    tuples = union(tuples, Set(reduce(append!, keys.(λ_filtered[i] for i in 1:length(λ_filtered)), init=T[])))

    return tuples
end
