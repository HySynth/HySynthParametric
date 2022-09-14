using ReachabilityAnalysis, DifferentialEquations, Plots, Plots.Measures

include("../models/diauxic_shift/diauxic_shift.jl")

function diauxic_shift_X0(; C1=nothing)
    x0 = zeros(8)
    # good initial states:
    # (C1, 1000, 20, 0.1, 0.01, 0, 0, 0)
    # C1 ∈ {50, 500}
    x0[1] = isnothing(C1) ? 50 : C1
    x0[2] = 1000
    x0[3] = 20
    x0[4] = 0.1
    x0[5] = 0.01
    return BallInf(x0, 0.005)
end

function diauxic_shift_X0_mode(; C1=nothing)
    X0 = diauxic_shift_X0(C1=C1)
    return [(2, X0)]  # mode 2 is "(on, on)"
end

function run_diauxic_shift(; k=1,
                          T=60.0,
                          seed=0,
                          change_points_and_slopes=cpds_arg_dist((x, y) -> norm(x - y), 0.5),
                          ε_mode=SingleEpsilonMode(),
                          C1=nothing,
                          debug::Bool=false)
    reset_seed(seed)
    H = diauxic()
    X0 = diauxic_shift_X0_mode(C1=C1)
    prob = IVP(H, X0)

    println("obtaining $k simulations from the original system...")
    @time tsv = simulate(prob, T, k)
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

function run_diauxic_shift_all(; k=1,
                              T=60.0,
                              seed=0,
                              plot_dim::Int=0,
                              ε_mode=SingleEpsilonMode(),
                              debug::Bool=true)
    H, tsv, ε, x0v, switchesv, mode_sequencev =
        run_diauxic_shift(T=T, k=k, seed=seed, ε_mode=ε_mode, debug=debug)

    tsv2 = simulate_imitate(H, tsv, x0v, switchesv, mode_sequencev)

    if plot_dim != 0
        fig = plot_results(tsv, tsv2, ε, d=plot_dim)
        return H, ε, tsv, tsv2, fig
    end

    return H, ε, tsv, tsv2
end

function benchmark_diauxic_shift(; k=5,
                                 T=60.0,
                                 seed=0,
                                 ε_mode=SingleEpsilonMode(),
                                 debug::Bool=false)
    Hs = []
    figs = []
    for C1 ∈ [50, 500]
        println("run with C1 = $C1")

        # synthesis
        H, tsv, ε, x0v, switchesv, mode_sequencev, runtime =
            run_diauxic_shift(T=T, k=k, seed=seed, C1=C1, ε_mode=ε_mode, debug=debug)
        println("resulting ε: $ε; number of locations: $(length(H.modes))")
        push!(Hs, H)

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
        for (d, tube_radius, ylab) in [(1, C1 == 50 ? 0.8 : 4, "C₁"),
                                       (2, 8, "C₂"),
                                       (3, C1 == 50 ? 1.5 : 3.5, "M")]
            fig = plot_results_individual([tsv[1]], [tsv2[1]], ε, d=d)
            plot_tube!(tsv[1], ε, c=:green, dim=d)
            [plot_tube!(ts3, tube_radius, c=:orange, dim=d, alpha=0.8) for ts3 in tsv3]
            plot_tube!(tsv2[1], tube_radius, c=:red, dim=d, alpha=0.8)
            plot!(size=(900, 400), tickfont=font(20, "Times"), guidefontsize=20,
                  bottom_margin=8mm, left_margin=5mm, ylab=ylab)
            push!(figs, fig)
            savefig(fig, "diauxic_shift_C1-$(C1)_$d.pdf")
        end
    end

    return Hs, figs
end
