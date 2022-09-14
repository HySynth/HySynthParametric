using HySynthParametric

# Running example (takes a few seconds)
println("------------------------------------------")
println("running 'Thermostat' (the running example)")
println("------------------------------------------")
include("running_example_thermostat.jl")
run_running_example_thermostat()

# Case study "Cell cycle" (takes ~5 minutes the first time)
println("------------------------------")
println("running 'Cell cycle' benchmark")
println("------------------------------")
include("cell_cycle.jl")
benchmark_cell_cycle(seed=10)

# Case study "Diauxic shift" (takes ~5 minutes the first time)
println("---------------------------------")
println("running 'Diauxic shift' benchmark")
println("---------------------------------")
include("diauxic_shift.jl")
benchmark_diauxic_shift(seed=10)

# Case study "Scalability" (takes several hours)
println("----------------------------------------")
println("running 'Multiple thermostats' benchmark")
println("----------------------------------------")
include("multi_thermostat.jl")
results = Dict()
benchmark_multi_thermostat!(results; warmup=true)
benchmark_multi_thermostat!(results)
plot_scalability_results(results)