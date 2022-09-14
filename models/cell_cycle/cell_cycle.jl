# ==================================================================================
# Cell-cycle model
# See https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1001077
# ==================================================================================
using LazySets, HybridSystems, MathematicalSystems, LinearAlgebra, SparseArrays
using LazySets: HalfSpace  # resolve name-space conflicts with Polyhedra

include("equations.jl")

## random sample from a normal (Gaussian) distribution
function rand_normal(mean, stdev)
    if stdev <= 0.0
        error("standard deviation must be positive")
    end
    u1 = rand()
    u2 = rand()
    r = sqrt( -2.0*log(u1) )
    theta = 2.0*pi*u2
    mean + stdev*r*sin(theta)
end

## delay of an exponential distribution given the mean
# (returns a negative number)
function delay_exponential(mean)
    return -mean
end

function wrap_mode(dyn, invariant)
    A, b = dyn
    return ConstrainedAffineContinuousSystem(A, b, invariant)
end

function cell()
    # variables
    CycA = 1	# cyclin A concentration
    CycB = 2	# cyclin B concentration
    CycE = 3	# cyclin E concentration
    M = 4		# cell mass
    t = 5		# time
    # number of variables
    n = 5


    # constants
    θA = 12.5
    θB1 = 21.25
    θB2 = 3.
    θE = 80.

    cM = 4.

    # random means (some of them are ignored)
    λ1 = 2
    λ2 = 0
    λ3 = 0.01
    λ4 = 1
    λ5 = 0.5
    λ6 = 0.75
    λ7 = 1.5
    λ8 = 0.5
    λ9 = 0.025

    # uniform random number r from [0, 1], the mean would be r = 0.5 and then ln(r) ~ -0.7
    deviation = -0.7

    # time delays
    Tr1 = delay_exponential(λ1) * deviation
    Tr4 = delay_exponential(λ4) * deviation
    Tr6 = delay_exponential(λ6) * deviation
    Tr7 = delay_exponential(λ7) * deviation
    Tr8 = delay_exponential(λ8) * deviation

    # bloating of equality constraints for floating-point issues
    ε = 1.0
    εt = 0.01	# bloating for time

    G = rand_normal(1.0, 3.33)
    δ = 0.5  # cell growth parameter (we use a deterministic value)

    # discrete structure (graph)
    automaton = GraphAutomaton(9)

    # mode 1 (G1a)
    F = G1a!
    invariant = HPolyhedron([
        HalfSpace(sparsevec([CycA], [-1.], n), 0.),  # CycA >= 0
        HalfSpace(sparsevec([CycB], [-1.], n), 0.),  # CycB >= 0
        HalfSpace(sparsevec([CycE], [-1.], n), 0.),  # CycE >= 0
        HalfSpace(sparsevec([M], [-1.], n), 0.),  # M >= 0
        HalfSpace(sparsevec([t], [1.], n), Tr1 + εt)	# t <= Tr1 + εt 
       ])
    m1 = wrap_mode(F, invariant)


    # mode 2 (early_G1b)
    F = early_G1b!
    invariant = HPolyhedron([
        HalfSpace(sparsevec([CycA], [-1.], n), 0.),  # CycA >= 0
        HalfSpace(sparsevec([CycB], [-1.], n), 0.),  # CycB >= 0
        HalfSpace(sparsevec([CycE], [-1.], n), 0.),  # CycE >= 0
        HalfSpace(sparsevec([M], [-1.], n), 0.),  # M >= 0
        HalfSpace(sparsevec([CycE], [cM], n), θE + ε)  # CycE * M <= θE + ε
       ])
    m2 = wrap_mode(F, invariant)


    # mode 3 (late_G1b)
    F = late_G1b!
    invariant = HPolyhedron([
        HalfSpace(sparsevec([CycA], [-1.], n), 0.),  # CycA >= 0
        HalfSpace(sparsevec([CycB], [-1.], n), 0.),  # CycB >= 0
        HalfSpace(sparsevec([CycE], [-1.], n), 0.),  # CycE >= 0
        HalfSpace(sparsevec([M], [-1.], n), 0.),  # M >= 0
        HalfSpace(sparsevec([CycA], [1.], n), θA + ε)  # CycA <= θA + ε
       ])
    m3 = wrap_mode(F, invariant)


    # mode 4 (S)
    F = S!
    invariant = HPolyhedron([
        HalfSpace(sparsevec([CycA], [-1.], n), 0.),  # CycA >= 0
        HalfSpace(sparsevec([CycB], [-1.], n), 0.),  # CycB >= 0
        HalfSpace(sparsevec([CycE], [-1.], n), 0.),  # CycE >= 0
        HalfSpace(sparsevec([M], [-1.], n), 0.),  # M >= 0
        HalfSpace(sparsevec([t], [1.], n), 7 + Tr4 + εt)  # t <= 7 + Tr4 + εt
       ])
    m4 = wrap_mode(F, invariant)


    # mode 5 (G2)
    F = G2!
    invariant = HPolyhedron([
        HalfSpace(sparsevec([CycA], [-1.], n), 0.),  # CycA >= 0
        HalfSpace(sparsevec([CycB], [-1.], n), 0.),  # CycB >= 0
        HalfSpace(sparsevec([CycE], [-1.], n), 0.),  # CycE >= 0
        HalfSpace(sparsevec([M], [-1.], n), 0.),  # M >= 0
        HalfSpace(sparsevec([CycB], [1.], n), θB1 + ε)  # CycB <= θB1 + ε
       ])
    m5 = wrap_mode(F, invariant)


    # mode 6 (prophase)
    F = prophase!
    invariant = HPolyhedron([
        HalfSpace(sparsevec([CycA], [-1.], n), 0.),  # CycA >= 0
        HalfSpace(sparsevec([CycB], [-1.], n), 0.),  # CycB >= 0
        HalfSpace(sparsevec([CycE], [-1.], n), 0.),  # CycE >= 0
        HalfSpace(sparsevec([M], [-1.], n), 0.),  # M >= 0
        HalfSpace(sparsevec([t], [1.], n), Tr6 + εt)  # t <= Tr6 + εt
       ])
    m6 = wrap_mode(F, invariant)


    # mode 7 (metaphase)
    F = metaphase!
    invariant = HPolyhedron([
        HalfSpace(sparsevec([CycA], [-1.], n), 0.),  # CycA >= 0
        HalfSpace(sparsevec([CycB], [-1.], n), 0.),  # CycB >= 0
        HalfSpace(sparsevec([CycE], [-1.], n), 0.),  # CycE >= 0
        HalfSpace(sparsevec([M], [-1.], n), 0.),  # M >= 0
        HalfSpace(sparsevec([t], [1.], n), Tr7 + εt)  # t <= Tr7 + εt
       ])
    m7 = wrap_mode(F, invariant)


    # mode 8 (anaphase)
    F = anaphase!
    invariant = HPolyhedron([
        HalfSpace(sparsevec([CycA], [-1.], n), 0.),  # CycA >= 0
        HalfSpace(sparsevec([CycB], [-1.], n), 0.),  # CycB >= 0
        HalfSpace(sparsevec([CycE], [-1.], n), 0.),  # CycE >= 0
        HalfSpace(sparsevec([M], [-1.], n), 0.),  # M >= 0
        HalfSpace(sparsevec([t], [1.], n), Tr8 + εt)  # t <= Tr8 + εt
       ])
    m8 = wrap_mode(F, invariant)


    # mode 9 (telophase)
    F = telophase!
    invariant = HPolyhedron([
        HalfSpace(sparsevec([CycA], [-1.], n), 0.),  # CycA >= 0
        HalfSpace(sparsevec([CycB], [-1.], n), 0.),  # CycB >= 0
        HalfSpace(sparsevec([CycE], [-1.], n), 0.),  # CycE >= 0
        HalfSpace(sparsevec([M], [-1.], n), 0.),  # M >= 0
        HalfSpace(sparsevec([CycB], [-1.], n), -θB2 + ε)  # CycB >= θB2 - ε
       ])
    m9 = wrap_mode(F, invariant)

    # modes
    modes = [m1, m2, m3, m4, m5, m6, m7, m8, m9]


    # common reset: time
    reset = Dict(t => 0.0)

    # transition 1 -> 2
    add_transition!(automaton, 1, 2, 1)
    guard = HPolyhedron([
        HalfSpace(sparsevec([t], [-1.], n), -(Tr1 - εt)),  # t >= Tr1 - εt
        HalfSpace(sparsevec([t], [1.], n), Tr1 + εt)  # t <= Tr1 + εt
       ])
    t1 = ConstrainedResetMap(n, guard, reset)

    # transition 2 -> 3
    add_transition!(automaton, 2, 3, 2)
    guard = HalfSpace(sparsevec([CycE], [-cM], n), -θE)  # CycE * M >= θE  !! M ~ cM (solve later)
    t2 = ConstrainedResetMap(n, guard, reset)


    # transition 3 -> 4
    add_transition!(automaton, 3, 4, 3)
    guard = HalfSpace(sparsevec([CycA], [-1.], n), θA)  # CycA >= θA
    t3 = ConstrainedResetMap(n, guard, reset)

    # transition 4 -> 5
    # DNA synthesis at least 7 hours
    add_transition!(automaton, 4, 5, 4)
    guard = HPolyhedron([
        HalfSpace(sparsevec([t], [-1.], n), -(7 + Tr4))  # t >= 7 + Tr4
       ])
    t4 = ConstrainedIdentityMap(n, guard)

    # transition 5 -> 6
    add_transition!(automaton, 5, 6, 5)
    guard = HalfSpace(sparsevec([CycB], [-1.], n), -θB1)  # CycB >= θB1
    t5 = ConstrainedResetMap(n, guard, reset)

    # transition 6 -> 7
    add_transition!(automaton, 6, 7, 6)
    guard = HalfSpace(sparsevec([t], [-1.], n), -Tr6)  # t >= Tr6
    t6 = ConstrainedResetMap(n, guard, reset)

    # transition 7 -> 8
    add_transition!(automaton, 7, 8, 7)
    guard = HalfSpace(sparsevec([t], [-1.], n), -Tr7)  # t >= Tr7
    t7 = ConstrainedResetMap(n, guard, reset)

    # transition 8 -> 9
    add_transition!(automaton, 8, 9, 8)
    guard = HalfSpace(sparsevec([t], [-1.], n), -Tr8)  # t >= Tr8
    t8 = ConstrainedResetMap(n, guard, reset)

    # transition 9 -> 1
    add_transition!(automaton, 9, 1, 9)
    guard = HalfSpace(sparsevec([CycB], [1.], n), θB2)  # CycB <= θB2
    # mass of daughter cells at division: M' = δ * M
    a = ones(n)
    a[M] = δ
    a[t] = 0
    A = Diagonal(a)
    c = zeros(n)
    t9 = ConstrainedAffineMap(A, c, guard)

    # transition annotations
    resetmaps = [t1, t2, t3, t4, t5, t6, t7, t8, t9]

    # switching
    switchings = [AutonomousSwitching()]

    H = HybridSystem(automaton, modes, resetmaps, switchings)

    return H
end