using ReachabilityAnalysis, DifferentialEquations, Plots, Plots.Measures
using ReachabilityAnalysis: Interval

include("../models/cell_cycle/cell_cycle.jl")

function cell_cycle_X0(; projection)
    X0 = Hyperrectangle(low=[1., 1., 1.1, 3., 0.], high=[1.1, 1.1, 1.1, 3.1, 0.])
    if !isnothing(projection)
        X0 = project(X0, projection)
    end
    return [(1, X0)]
end

function simulate_cell_cycle(; k=1, T=30.0)
    H = cell()
    X0 = cell_cycle_X0(projection=nothing)
    prob = IVP(H, X0)
    proj = 1:4  # project to state variables, i.e., forget auxiliary variable t
    tsv = simulate(prob, T, k, projection=proj)
    return tsv
end

function run_cell_cycle(; k=1,
                       T=30.0,
                       seed=0,
                       change_points_and_slopes=cpds_arg_dist((x, y) -> norm(x - y), 0.4),
                       ε_mode=SingleEpsilonMode(),
                       debug::Bool=false)
    reset_seed(seed)

    println("obtaining $k simulations from the original system...")
    @time tsv = simulate_cell_cycle(k=k, T=T)
    n_datapoints = sum(length(ts) for ts in tsv)
    println("...done; $n_datapoints data points in total")

    println("running synthesis from $k time series...")
    res = @timed analyze(tsv; change_points_and_slopes=change_points_and_slopes,
                         ε_mode=ε_mode, debug=debug)
    H, ε, x0v, switchesv, mode_sequencev = res.value
    runtime = res.time
    println("...done after $runtime seconds")

    return H, tsv, ε, x0v, switchesv, mode_sequencev, runtime
end

function run_cell_cycle_all(; k=1,
                           T=30.0,
                           seed=0,
                           plot_dim::Int=0,
                           ε_mode=SingleEpsilonMode(),
                           debug::Bool=true)
    H, tsv, ε, x0v, switchesv, mode_sequencev, runtime =
        run_cell_cycle(T=T, k=k, seed=seed, ε_mode=ε_mode, debug=debug)

    tsv2 = simulate_imitate(H, tsv, x0v, switchesv, mode_sequencev)

    if plot_dim != 0
        fig = plot_results(tsv, tsv2, ε, d=plot_dim)
        return H, ε, tsv, tsv2, fig
    end

    return H, ε, tsv, tsv2
end

function benchmark_cell_cycle(; k=20,
                              T=30.0,
                              seed=0,
                              ε_mode=SingleEpsilonMode(),
                              debug::Bool=false)
    # synthesis
    H, tsv, ε, x0v, switchesv, mode_sequencev, runtime =
        run_cell_cycle(T=T, k=k, seed=seed, ε_mode=ε_mode, debug=debug)
    println("resulting ε: $ε; number of locations: $(length(H.modes))")

    # imitation
    tsv2 = simulate_imitate(H, tsv, x0v, switchesv, mode_sequencev)

    # simulation of new model (search for simulations until the time horizon)
    X0 = [(mode_sequence[1], BallInf(x0, 0.01)) for (mode_sequence, x0) in zip(mode_sequencev, x0v)]
    prob = IVP(H, X0)
    tsv3 = []
    while length(tsv3) < 3
        ts = simulate(prob, T, 1)[1]
        if ts.time[end] ≈ T
            push!(tsv3, ts)
        end
    end

    # plot
    figs = []
    for (d, tube_radius, ylab) in [(1, 0.5, "CycA"), (2, 0.2, "CycB"), (3, 0.2, "CycE")]
        fig = plot_results_individual([tsv[1]], [tsv2[1]], ε, d=d)
        plot_tube!(tsv[1], ε, c=:green, dim=d)
        [plot_tube!(ts3, tube_radius, c=:orange, dim=d, alpha=0.8) for ts3 in tsv3]
        plot_tube!(tsv2[1], tube_radius, c=:red, dim=d, alpha=0.8)
        plot!(size=(900, 400), tickfont=font(20, "Times"), guidefontsize=20,
              bottom_margin=8mm, left_margin=5mm, ylab=ylab)
        push!(figs, fig)
        savefig(fig, "cell_cycle_$d.pdf")
    end

    return H, ε, tsv, tsv2, tsv3, figs
end
