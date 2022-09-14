using HybridSystems, LinearAlgebra
using HybridSystems: GraphTransition

# construct an automaton from the following inputs:
# - a vector of time series
# - a vector of sequences of switches
# - a vector of mode sequences
# - a bloating parameter for invariants and guards
# - the computed slope parameters
function construct_automaton(tsv::Vector{TimeSeries}, switchesv, mode_sequencev,
                             δ::Float64, slopesv)
    # dimension
    n = dim(tsv[1])

    # number of modes
    m = length(slopesv)

    # create automaton (discrete structure)
    automaton = GraphAutomaton(m)

    # add modes
    # collect all data points for each mode
    mode2points = [Vector{Float64}[] for _ in 1:m]
    for (ts, mode_sequence, switches) in zip(tsv, mode_sequencev, switchesv)
        K = length(mode_sequence)
        from = 1
        for k in 1:K
            mode = mode_sequence[k]
            to = k == K ? length(ts) : switches[k]
            append!(mode2points[mode], ts.data[from:to])
            from = to
        end
    end
    modes = []
    A = zeros(n, n)
    for (j, slope) in enumerate(slopesv)
        # create invariant
        points = mode2points[j]
        invariant = box_approximation(VPolytope(points))
        invariant = _bloat(invariant, δ)

        # create an affine ODE
        b = slope
        mode = @system(x' = A*x + b, x ∈ invariant)
        push!(modes, mode)
    end

    # add transitions
    transitions = Dict{Tuple{Int, Int}, GraphTransition}()
    resetmaps = []
    symb = 1
    for (ts, mode_sequence, switches) in zip(tsv, mode_sequencev, switchesv)
        source = mode_sequence[1]
        for (switch, target) in enumerate(mode_sequence[2:end])
            if source != target  # do not add self-loops
                guard = _get_constraint(ts.data[switches[switch]], δ)
                transition = get(transitions, (source, target), nothing)
                @assert has_transition(automaton, source, target) == !isnothing(transition)
                if !isnothing(transition)
                    # transition exists => find index
                    idx = symbol(automaton, transition)

                    # update guard
                    resetmap = resetmaps[idx]
                    old_guard = stateset(resetmap)
                    guard = _extend_constraint(old_guard, guard)
                    resetmap = ConstrainedIdentityMap(n, guard)
                    resetmaps[idx] = resetmap
                else
                    # create transition
                    transition = add_transition!(automaton, source, target, symb)
                    transitions[(source, target)] = transition

                    # create guard
                    resetmap = ConstrainedIdentityMap(n, guard)
                    push!(resetmaps, resetmap)
                    symb += 1
                end
            end
            source = target
        end
    end

    # create hybrid automaton
    H = HybridSystem(automaton, modes, resetmaps, [AutonomousSwitching()])

    return H
end

function _get_constraint(p, δ)
    return BallInf(p, δ)
end

function _extend_constraint(X, Y)
    return box_approximation(ConvexHull(X, Y))
end

function _bloat(X, δ)
    return minkowski_sum(X, BallInf(zeros(dim(X)), δ))
end
