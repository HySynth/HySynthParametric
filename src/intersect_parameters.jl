using LazySets: SingleEntryVector
const SEV = SingleEntryVector

function intersect_parameters_shared_x0(Ps, mPs; remove_redundant_constraints::Bool=false)
    n = dim(Ps[1])
    Q = nothing
    dims = vcat(1:mPs, n)
    @inbounds for (i, Pi) in enumerate(Ps)
        # many constraints are redundant and removing them for the individual
        # P's is cheaper (but may still be too expensive)
        if remove_redundant_constraints
            LazySets.remove_redundant_constraints!(Pi)
        end

        Pp = project(Pi, dims)
        if isnothing(Q)
            Q = Pp
        else
            # removing redundant constraints here is expensive
            Q = intersection(Q, Pp; prune=false)
        end
    end
    return Q
end

function intersect_parameters_disjoint_x0(Ps; remove_redundant_constraints::Bool=false)
    Q = nothing
    @inbounds for P in Ps
        # many constraints are redundant and removing them for the individual
        # P's is cheaper (but may still be too expensive)
        if remove_redundant_constraints
            LazySets.remove_redundant_constraints!(P)
        end

        if isnothing(Q)
            Q = P
        else
            # there are no redundant constraints after taking the intersection
            # because the x0's differ
            Q = intersection(Q, P; prune=false)
        end
    end
    return Q
end
