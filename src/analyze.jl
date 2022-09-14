# inputs:
# - a vector of time series
# - a function to compute change points and corresponding slopes
# - a function to cluster the slopes
#
# output:
# the parameter polyhedra for each time series
function analyze(tsv::AbstractVector{TimeSeries};
                 change_points_and_slopes,
                 kmeans_args="default",
                 cluster_predicate=nothing,
                 use_shared_x0::Bool=false,
                 simplify_δ=SIMPLIFY_THRESHOLD,
                 normalize_slopes::Bool=false,
                 cluster_delays::Bool=false,
                 ε_mode::EpsilonMode=SingleEpsilonMode(),
                 debug::Bool=true)
    # get switches and corresponding preliminary slopes
    switchesv = Vector{Int}[]
    slopesv = Vector{Vector{Float64}}[]
    delaysv = Vector{Float64}[]
    for (i, ts) in enumerate(tsv)
        # eliminate linear parts in time series
        simplify_ts_linear_interpolation!(ts; δ=simplify_δ, debug=(debug ? i : 0))

        switches, slopes = change_points_and_slopes(ts)

        # optionally normalize preliminary slopes
        if normalize_slopes
            slopes = normalize.(slopes)
        end
        push!(switchesv, switches)
        push!(slopesv, slopes)

        if cluster_delays
            delays = get_delays(ts, switches)
            push!(delaysv, delays)
	end
    end
    debug && println("switchesv: $switchesv")

    # flatten nested vector of slopes
    slopes = vcat(slopesv...)
    debug && println("slopes: $slopes")
    delays = vcat(delaysv...)

    # cluster slopes and get corresponding mode assignments
    if isnothing(kmeans_args)
        @warn "delay clustering is not supported in this mode"
        clusters, mode_sequence_all = cluster(slopes, cluster_predicate)
    else
        clusters, mode_sequence_all = cluster_kmeans(slopes, kmeans_args;
                                                     delays=delays, debug=debug)
    end
    debug && println("clusters: $clusters")

    mode_sequencev = extract_mode_sequences(switchesv, mode_sequence_all)
    debug && println("mode_sequencev: $mode_sequencev")

    # merge consecutive pieces with the same mode
    for (switches, mode_sequence) in zip(switchesv, mode_sequencev)
        k = length(mode_sequence)
        merge_same_slopes!(switches, mode_sequence)
        debug && println("saved $(k - length(mode_sequence))/$k redundant modes")
    end

    # data dimension
    n = dim(tsv[1])

    # obtain parameter polyhedra
    Ps = LazySet[]
    n_modes = length(clusters)
    n_slopes = n * n_modes  # dimensions used for the slopes
    n_x0 = n * (use_shared_x0 ? 1 : length(tsv))  # dimensions used for x0
    n_without_ε = n_slopes + n_x0  # dimensions used for slopes and x0
    n_total = _get_n_total(ε_mode, n_without_ε, tsv)  # all dimensions
    ε_offset = 0
    for (i, (ts, switches, mode_sequence)) in enumerate(zip(tsv, switchesv, mode_sequencev))
        idx = use_shared_x0 ? 1 : i
        P = parameter_range(ts, switches, mode_sequence, n_without_ε, n_total;
                            p=n_modes, idx=idx, ε_mode=ε_mode, ε_offset=ε_offset)
        push!(Ps, P)
        ε_offset = _increase_ε_offset(ε_mode, ts, ε_offset)
#         debug && println("P$i: $P")
    end


    # intersect parameter polyhedra
    debug && println("computing Pcap...")
    if use_shared_x0
        @time Pcap = intersect_parameters_shared_x0(Ps, n_slopes)
    else
        @time Pcap = intersect_parameters_disjoint_x0(Ps)
    end
    debug && println("...done")

    # obtain slopes together with minimizing ε
    debug && println("minimizing ε and extracting solution...")
    @time if use_shared_x0
        ε, all_slopes = minimizing_slopes_projected(Pcap, n_slopes)
        x0v = minimizing_x0_by_substitution(Ps, ε, all_slopes; ε_mode=ε_mode)
    else
        ε, all_slopes, x0s = minimizing_slopes_x0(Pcap, n_slopes, n_without_ε;
                                                  ε_mode=ε_mode)
        x0v = extract_vectors(x0s, n)
    end
    debug && println("...done")
    slopesv = extract_vectors(all_slopes, n)
    debug && println("ε: $ε")
    debug && println("slopesv: $slopesv")
    debug && println("x0v: $x0v")

    H = construct_automaton(tsv, switchesv, mode_sequencev, ε, slopesv)

    return H, ε, x0v, switchesv, mode_sequencev
end

function _get_n_total(ε_mode::SingleEpsilonMode, n_without_ε::Int,
                      tsv::AbstractVector{TimeSeries})
    return n_without_ε + 1
end

function _get_n_total(ε_mode::MultiEpsilonMode, n_without_ε::Int,
                      tsv::AbstractVector{TimeSeries})
    res = n_without_ε
    for ts in tsv
        res += length(ts)
    end
    return res
end

function _increase_ε_offset(ε_mode::SingleEpsilonMode, ts::TimeSeries, ε_offset::Int)
    return 1
end

function _increase_ε_offset(ε_mode::MultiEpsilonMode, ts::TimeSeries, ε_offset::Int)
    return ε_offset + length(ts.data)
end
