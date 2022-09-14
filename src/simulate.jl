using ReachabilityAnalysis
import Polyhedra, CDDLib
using LazySets: center, dim
import DifferentialEquations
using ReachabilityAnalysis: HalfSpace

function LazySets.sample(X0::Vector{<:Tuple{Int, <:LazySet}}, k::Int=1)
    samples = []
    for _ in 1:k
        mode, X0_mode = rand(X0)
        x0 = sample(X0_mode)
        push!(samples, (mode, x0))
    end
    return samples
end

# k simulations of an IVP with a hybrid automaton H until time horizon T
#
# the result is a vector of TimeSeries objects
#
# `projection` can be used to eliminate some dimensions
function simulate(ivp::IVP{<:HybridSystem}, T::Number, k::Int=1;
                  projection=nothing, dtmax=DTMAX, n_points=Inf)
    sims = ReachabilityAnalysis._solve_ensemble(ivp, T=T, trajectories=k,
                                                dtmax=dtmax)

    # postprocessing: remove jump entries with (almost) zero time distance
    # (these exist because in general a hybrid automaton may have discontinuous
    # jumps)
    tsv = Vector{TimeSeries}(undef, length(sims))
    for (j, s) in enumerate(sims)
        u_proj = _project(s.u, projection)
        t = Vector{eltype(s.t)}()
        u = Vector{eltype(u_proj)}()
        n_points_used = 0
        for i in eachindex(s.t)
            if i > 1 && s.t[i] ≈ s.t[i-1]
#                 if !(u_proj[i] ≈ u_proj[i-1])
#                     @warn "trajectory contains discontinuous jumps; this " *
#                         "requires special treatment"
#                 end
                continue
            end
            push!(t, s.t[i])
            push!(u, u_proj[i])
            n_points_used += 1
            if n_points_used == n_points
                break
            end
        end
        tsv[j] = TimeSeries(t, u)
    end
    return tsv
end

function _project(v, proj)
    return [vi[proj] for vi in v]
end

function _project(v, proj::Nothing)
    return v
end

# "simulate" a hybrid automaton given the initial states, the mode sequence, and
# the switching times (i.e., the result is deterministic)
function simulate_imitate(H::HybridSystem, tsv, x0v, switchesv, mode_sequencev)
    tsv2 = Vector{TimeSeries}(undef, length(tsv))
    @inbounds for (i, (ts, x0, switches, mode_sequence)) in enumerate(zip(tsv, x0v, switchesv, mode_sequencev))
        data = Vector{Vector{Float64}}()
        time = Vector{Float64}()
        t0 = 0.0
        @assert length(switches) == length(mode_sequence) - 1
        for (j, mode_idx) in enumerate(mode_sequence)
            if j < length(mode_sequence)
                mode_idx_next = mode_sequence[j+1]
                if mode_idx_next == mode_idx
                    # stayed in the same mode
                    continue
                end
                switch = switches[j]
                t1 = ts.time[switch]
            else
                t1 = ts.time[end]
            end

            m = mode(H, mode_idx)
            slope = affine_term(m)
            x1 = x0 + slope * (t1 - t0)

            # sanity checks
            inv = stateset(m)
            @assert x0 ∈ inv && x1 ∈ inv "left the invariant during imitation"
            if j < length(mode_sequence)
                trs = transitions(H, mode_idx, mode_idx_next)
                @assert length(trs) == 1
                grd = guard(H, first(trs))
                @assert x1 ∈ grd "did not hit the guard during imitation"
            end

            push!(data, x0)
            push!(time, t0)
            x0 = x1
            t0 = t1
        end
        push!(data, x0)
        push!(time, t0)
        tsv2[i] = TimeSeries(time, data)
    end
    return tsv2
end

######################################################
# old simulation code (which does not simulate well) #
######################################################

# k simulations of an IVP with a hybrid automaton H until time horizon T
#
# the result is a vector of TimeSeries objects
function simulate2(ivp::IVP{<:HybridSystem}, T::Number, k::Int=1; δ::Number=TIME_STEP_BOX_REACH)
    return simulate2(ivp.s, ivp.x0, T, k; δ=δ)
end

# k simulations of a hybrid automaton H from initial set X0 until time horizon T
#
# the result is a vector of TimeSeries objects
function simulate2(H::HybridSystem,
                  X0::Union{LazySet, Vector{<:Tuple{Int, <:LazySet}}},
                  T::Number,
                  k::Int=1;
                  δ::Number=TIME_STEP_BOX_REACH)
    samples = sample(X0, k)
    tsv = [simulate2(H, x0, T; δ=δ) for x0 in samples]
    return tsv
end

# one simulation of a hybrid automaton H from initial set X0 until time horizon T
#
# the result is a TimeSeries object
function simulate2(H::HybridSystem, x0::Tuple{Int, Vector{N}}, T::Number;
                  δ::Number=TIME_STEP_BOX_REACH) where {N<:Number}
    # compute flowpipe
    X0 = [(x0[1], Singleton(x0[2]))]
    prob = IVP(H, X0)
    alg = BOX(δ=δ)
    approx_model = NoBloating()
    Fv = solve(prob, alg=alg, approx_model=approx_model, tspan=(0.0 .. T))

    time = Vector{N}()
    data = Vector{Vector{N}}()
    for (i, F) in enumerate(Fv)
        # extract time
        time_i = [tstart(tspan(R)) for R in F]
        if i > 1
            # remove elements from time series after start of current flowpipe
            skip = length(findall(ti -> ti <= time[end], time_i))
            time = time[1:end-skip]
            data = data[1:end-skip]
        end
        append!(time, time_i)

        # extract data
        sets = [set(R) for R in F]
        if dim(sets[1]) == 1
            # (chebyshev_center of intervals crashes, hence the workaround)
            sets = [convert(Interval, X) for X in sets]
        end
        if i == 1
            # (chebyshev_center of singletons crashes, hence the workaround)
            data_i = vcat([x0[2]], [chebyshev_center(X) for X in sets[2:end]])
        else
            data_i = [chebyshev_center(X) for X in sets]
        end
        append!(data, data_i)
    end

    # construct TimeSeries object
    return TimeSeries(time, data)
end
