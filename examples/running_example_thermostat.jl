using ReachabilityAnalysis, Symbolics, DifferentialEquations, Plots,
      Plots.PlotMeasures
using ReachabilityAnalysis: Interval

include("../models/thermostat/single_thermostat.jl")

function thermostat_X0()
    X0 = Interval(68.0, 69.0)
    return [(1, X0)]
end

function run_running_example_thermostat(; seed=0, debug=true)
    reset_seed(seed)
    H = single_thermostat()
    X0 = thermostat_X0()
    prob = IVP(H, X0)
    T = 5.0
    n_ts = 2
    dtmax = 2.0
    simplify_δ = 0.007

    debug && println("obtaining $n_ts simulations from the original system...")
    @time tsv = simulate(prob, T, n_ts; dtmax=dtmax)
    debug && println("...done")

    # run with dynamic clustering
    tsv_copy = deepcopy(tsv)
    analyze(tsv_copy;
            simplify_δ=simplify_δ,
            change_points_and_slopes=cpds_arg_all_pieces,
            debug=debug)

    debug && println()

    # run with fixed clustering for 2 modes
    kmeans_args = (2, SqEuclidean(), CLUSTERING_THRESHOLD)
    tsv_copy = deepcopy(tsv)
    H, ε, x0v, switchesv, mode_sequencev =
        analyze(tsv;
                simplify_δ=simplify_δ,
                change_points_and_slopes=cpds_arg_all_pieces,
                kmeans_args=kmeans_args,
                debug=debug)

    tsv2 = simulate_imitate(H, tsv, x0v, switchesv, mode_sequencev)

    tsv11 = [tsv[1]]
    tsv12 = [tsv[2]]
    tsv21 = [tsv2[1]]
    tsv22 = [tsv2[2]]
    colors = [:blue, :orange]
    markers = [:utriangle, :dtriangle]
    ms = 8

    println()
    data1 = [d[1] for d in tsv[1].data]
    println("time series 1:\nt = $(round.(tsv[1].time, digits=2))\nd = $(round.(data1, digits=2))")
    println("Δt:$(round.([tsv[1].time[i+1] - tsv[1].time[i] for i in 1:length(tsv[1])-1], digits=2))")
    data2 = [d[1] for d in tsv[2].data]
    println("time series 2:\nt = $(round.(tsv[2].time, digits=2))\nd = $(round.(data2, digits=2))")
    println("Δt:$(round.([tsv[2].time[i+1] - tsv[2].time[i] for i in 1:length(tsv[2])-1], digits=2))")
    println()
    println("loc 1:\nx' = $(H.modes[1].c), x ∈ $(convert(Interval, H.modes[1].X).dat)")
    println("loc 2:\nx' = $(H.modes[2].c), x ∈ $(convert(Interval, H.modes[2].X).dat)")
    println("trans 1-2:\nx ∈ $(convert(Interval, H.resetmaps[1].X).dat)")
    println("trans 2-1:\nx ∈ $(convert(Interval, H.resetmaps[2].X).dat)")

    fig1 = plot(size=(900, 400), tickfont=font(20, "Times"), guidefontsize=20,
                bottom_margin=8mm, left_margin=4mm)
    [plot_ts!(ts; c=c, m=m, ms=ms, ylab="x") for (ts, c, m) in zip(tsv, colors, markers)]
    [plot_tube!(ts, TUBE_RADIUS; c=c, ylab="x") for (ts, c) in zip(tsv, colors)]
    savefig(fig1, "thermostat_execution_original.pdf")

    fig2 = plot(size=(900, 400), tickfont=font(20, "Times"), guidefontsize=20,
                bottom_margin=8mm, left_margin=4mm)
    fig2 = plot_results_individual(tsv11, tsv21, ε, d=1, ylab="x", fig=fig2)
    plot_ts!(tsv11[1], ylab="x", m=markers[1], ms=ms, c=colors[1])
    savefig(fig2, "thermostat_execution_synthesized_1.pdf")

    fig3 = plot(size=(900, 400), tickfont=font(20, "Times"), guidefontsize=20,
                bottom_margin=8mm, left_margin=4mm)
    fig3 = plot_results_individual(tsv12, tsv22, ε, d=1, ylab="x", fig=fig3)
    plot_ts!(tsv12[1], ylab="x", m=markers[2], ms=ms, c=colors[2])
    savefig(fig3, "thermostat_execution_synthesized_2.pdf")

    # run with fixed clustering for 4 and 6 modes
    println()
    for λ in [4, 6]
        kmeans_args = (λ, SqEuclidean(), CLUSTERING_THRESHOLD)
        tsv_copy2 = deepcopy(tsv_copy)
        H2, ε2, x0v2, switchesv2, mode_sequencev2 =
            analyze(tsv_copy2;
                    simplify_δ=simplify_δ,
                    change_points_and_slopes=cpds_arg_all_pieces,
                    kmeans_args=kmeans_args,
                    debug=false)

        println("ε for λ = $λ locations: $ε2")

        tsv2_2 = simulate_imitate(H2, tsv_copy2, x0v2, switchesv2, mode_sequencev2)

        tsv11 = [tsv_copy2[1]]
        tsv12 = [tsv_copy2[2]]
        tsv21 = [tsv2_2[1]]
        tsv22 = [tsv2_2[2]]

        fig4 = plot(size=(900, 400), tickfont=font(20, "Times"), guidefontsize=20,
                    bottom_margin=8mm, left_margin=4mm)
        fig4 = plot_results_individual(tsv11, tsv21, ε2, d=1, ylab="x", fig=fig4)
        plot_ts!(tsv11[1], ylab="x", m=markers[1], ms=ms, c=colors[1])
        savefig(fig4, "thermostat_execution_synthesized_$(λ)_locs.pdf")
    end

    return H, ε, tsv, tsv2, fig1, fig2, fig3
end
