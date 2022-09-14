import Clustering
using Clustering: kmeans, SqEuclidean
# import Distances
# using Distances: SphericalAngle

function cluster_kmeans(slopes, kmeans_args="default";
                        delays=[], debug::Bool=false)
    # `kmeans_args` must be a tuple with
    # - the number of clusters, or `0` for a refinement procedure
    # - the distance function
    # - the threshold for relative improvement in case of refinement
    #   (ignored if the number of clusters is fixed)
    if kmeans_args == "default"
        kmeans_args = (0, SqEuclidean(), CLUSTERING_THRESHOLD)  # default values
    end
    k, distance, threshold = kmeans_args

    if !isempty(delays)
        # add delay dimension
        @assert length(slopes) == length(delays)
        slopes = [vcat(slopes[i], 10*delays[i]) for i in 1:length(slopes)]
    end

    X = hcat(slopes...)  # align slopes as columns in a matrix
    if k == 0
        # refinement loop
        prev_cost = Inf
        for k in 1:length(slopes)
            kmeans_res = kmeans(X, k, init=1:k, distance=distance)
            cost = kmeans_res.totalcost
            rel_cost = isinf(prev_cost) ? Inf : (prev_cost - cost) / prev_cost
            debug && println("error for $k clusters: $cost (relative: $rel_cost)")
            if rel_cost < threshold
                break
            end
            prev_cost = kmeans_res.totalcost
        end
    else
        kmeans_res = kmeans(X, k, init=1:k, distance=distance)
    end

    clusters = columns(kmeans_res.centers)
    mode_sequence = kmeans_res.assignments

    if !isempty(delays)
        # project out delay dimension again
        clusters = [cl[1:end-1] for cl in clusters]
    end

    return clusters, mode_sequence
end

# helper function to enumerate all column vectors of a matrix
columns(M) = [M[:, j] for j in 1:size(M, 2)]

#####################
# manual clustering #
#####################

function cluster(elements::AbstractVector{T},
                 issimilar,
                 clusters=Vector{Vector{T}}()) where {T}
    mapping = Vector{Int}(undef, length(elements))
    for (i, E) in enumerate(elements)
        found_cluster = false
        for (j, cluster) in enumerate(clusters)
            F = first(cluster)
            if issimilar(E, F)
                # add to existing cluster
                push!(cluster, E)
                @inbounds mapping[i] = j
                found_cluster = true
                break
            end
        end

        if !found_cluster
            # create new cluster
            m = length(clusters)
            cluster = Vector{T}()
            push!(cluster, E)
            push!(clusters, cluster)
            @inbounds mapping[i] = m + 1
        end
    end
    return clusters, mapping
end

# two arrays are always similar
function issimilar_always(u::AbstractArray, v::AbstractArray)
    return true
end

# two arrays are never similar
function issimilar_never(u::AbstractArray, v::AbstractArray)
    return false
end

# two arrays are similar if the entries have the same sign
function issimilar_sign(u::AbstractArray, v::AbstractArray)
    return sign.(u) == sign.(v)
end

# two arrays are similar if the absolute difference of the entries is below `d`
function issimilar_abs(d::Number)
    function _issimilar_abs(u::AbstractArray, v::AbstractArray)
        return @inbounds all(abs(u[i] - v[i] <= d) for i in length(u))
    end
    return _issimilar_abs
end
