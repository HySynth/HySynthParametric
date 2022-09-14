using LazySets: project

# input:
# - a parameter polyhedron
#
# output:
# - the minimizing value ε
function minimize(P::HPolyhedron)
    # minimize P in ε direction
    n_total = dim(P)
    d = _objective(SingleEpsilonMode(), n_total, 0)
    ε = -ρ(d, P)
    return ε
end

# input:
# - a parameter polyhedron
# - number of slope parameters
# - number of slope and x0 parameters
# - mode for ε
#
# output:
# - the minimizing value ε
# - some minimizing slope values
# - some minimizing x0
function minimize(P::HPolyhedron, n_slopes, n_without_ε;
                  ε_mode=SingleEpsilonMode())
    # minimize P in ε direction
    n_total = dim(P)
    d = _objective(ε_mode, n_total, n_without_ε)
    dims = dim.(constraints_list(P))
    slopes_x0_ε = σ(d, P)

    slopes = slopes_x0_ε[1:n_slopes]
    x0 = slopes_x0_ε[n_slopes+1:n_without_ε]
    ε = _get_ε(ε_mode, slopes_x0_ε, n_without_ε)
    return ε, slopes, x0
end

# input:
# - a parameter polyhedron
# - the number of dimensions representing the slopes
# - a value ε
#
# output:
# - the polyhedron representing the feasible slopes corresponding to ε
function feasible_slopes(P::HPolyhedron, m::Int, ε::Number)
    n = dim(P)

    # select region with minimal ε
    a = LazySets.SingleEntryVector(n, n, 1.0)
    Pε = intersection(P, HalfSpace(a, ε))

    # project to slopes
    slopes = project(Pε, 1:m)

    return slopes
end

# input:
# - a parameter polyhedron
# - the number of dimensions representing the slopes
#
# output:
# - the minimizing value ε
# - the polyhedron representing the feasible slopes corresponding to ε
function minimizing_slopes(P::HPolyhedron, m::Int)
    # this code, although technically almost equivalent, does not work for the
    # cellCycle model
    # it computes the minimal ε and then fixes that value in P using the
    # constraint 'x_ε <= min(ε)'
    ε = minimize(P)
    slopes = feasible_slopes(P, m, ε)

    return ε, slopes
end

function minimizing_slopes_x0(P::HPolyhedron, n_slopes::Int,
                              n_without_ε::Int; ε_mode)
    ε, slopes, x0 = minimize(P, n_slopes, n_without_ε; ε_mode=ε_mode)
    return ε, slopes, x0
end

function minimizing_slopes_projected(P::HPolyhedron, n_slopes::Int; ε_mode)
    ε, slopes, x0 = minimize(P, n_slopes, n_slopes; ε_mode=ε_mode)
    return ε, slopes
end

# the `use_hyperplanes` option uses an alternative but slower implementation
function minimizing_x0_by_substitution(Ps, ε, all_slopes;
                                       use_hyperplanes::Bool=false)
    x0v = Vector{Vector{Float64}}(undef, length(Ps))
    n = dim(Ps[1])
    m = length(all_slopes)
    if !use_hyperplanes
        vec = vcat(all_slopes, zeros(n - m - 1), ε)
        proj = m+1:n-1
    end
    for (i, Pi) in enumerate(Ps)
        if use_hyperplanes
            # substitute values by intersecting with corresponding hyperplanes
            Q = copy(Pi)
            for (j, slope) in enumerate(all_slopes)
                Q = intersection(Q, Hyperplane(SingleEntryVector(j, n, 1.0), slope))
            end
            Q = intersection(Q, Hyperplane(SingleEntryVector(n, n, 1.0), ε))
            x0v[i] = an_element(Q)[m+1:n-1]
        else
            # substitute values by extracting matrix and multiplying with vector
            A, b = tosimplehrep(Pi)
            A2 = A[:, proj]
            b2 = b - A * vec
            Q = HPolyhedron(A2, b2)
            x0v[i] = an_element(Q)
        end
    end
    return x0v
end

function _objective(ε_mode::SingleEpsilonMode, n_total::Int, n::Int)
    return LazySets.SingleEntryVector(n_total, n_total, -1.0)
end

function _objective(ε_mode::MultiEpsilonMode, n_total::Int, n::Int)
    d = zeros(n_total)
    @inbounds for i in (n+1):n_total
        d[i] = -1.0
    end
    return d
end

function _get_ε(ε_mode::SingleEpsilonMode, slopes_x0_ε, n)
    return slopes_x0_ε[end]
end

function _get_ε(ε_mode::MultiEpsilonMode, slopes_x0_ε, n)
    return maximum(slopes_x0_ε[n+1:end])
end
