__precompile__(true)

module HySynthParametric

using HybridSystems, LazySets, LinearAlgebra, Plots, ReachabilityAnalysis
import CDDLib, Clustering, DelimitedFiles, DifferentialEquations, Distances,
       Polyhedra, Random

export reset_seed, simulate, analyze, simulate_imitate, cpds_arg_all_pieces,
       SqEuclidean, CLUSTERING_THRESHOLD, TUBE_RADIUS, SingleEpsilonMode,
       cpds_arg_dist,
       plot_ts!, plot_tube!, plot_results_individual

include("constants.jl")
include("util.jl")

include("simplify_time_series.jl")
include("clustering.jl")
include("merge_same_slopes.jl")
include("parameter_range.jl")
include("intersect_parameters.jl")
include("minimize.jl")
include("construct_automaton.jl")
include("analyze.jl")

include("simulate.jl")

end  # module
