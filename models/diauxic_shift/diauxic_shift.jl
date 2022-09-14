# =============================================
# Diauxic shift model
# See https://pubmed.ncbi.nlm.nih.gov/31342219/
# =============================================
using HybridSystems, MathematicalSystems, SparseArrays, TaylorIntegration
using LazySets: HalfSpace  # resolve name-space conflicts with Polyhedra

include("equations.jl")

const CBBCS = ConstrainedBlackBoxContinuousSystem

"""
    diauxic_shift()
Construct the bacterial chemotaxis model.
"""
function diauxic()
    # variables
    C1 = 1  #
    C2 = 2  #
    M = 3  #
    Q = 4  #
    R = 5  #
    T1 = 6  #
    T2 = 7  #
    RP = 8  #

    # number of variables
    n = 8

    # constants
    γ = 20.
    α = 0.03

    # discrete structure (graph)
    automaton = GraphAutomaton(4)

    # mode 1 ((RP,T2) = (off, on))
    F = diauxic_off_on!
    invariant = HPolyhedron([
        HalfSpace(sparsevec([C1], [-1.], n), 0.),  # C1 >= 0
        HalfSpace(sparsevec([C1], [1.], n), γ),  # C1 < γ
        HalfSpace(sparsevec([RP], [-1.], n), 0.),  # RP >= 0
        HalfSpace(sparsevec([RP], [1.], n), α),  # RP < α
        HalfSpace(sparsevec([C2], [-1.], n), 0.),  # C2 >= 0
        HalfSpace(sparsevec([M], [-1.], n), 0.),  # M >= 0
        HalfSpace(sparsevec([T1], [-1.], n), 0.),  # T1 >= 0
        HalfSpace(sparsevec([T2], [-1.], n), 0.),  # T2 >= 0
        HalfSpace(sparsevec([R], [-1.], n), 0.),  # R >= 0
        HalfSpace(sparsevec([Q], [-1.], n), 0.)  # Q >= 0
       ])
    m1 = CBBCS(F, n, invariant)

    # mode 2 ((RP,T2) = (on, on))
    F = diauxic_on_on!
    invariant = HPolyhedron([
        HalfSpace(sparsevec([C1], [-1.], n), -γ),  # C1 >= γ
        HalfSpace(sparsevec([RP], [-1.], n), 0.),  # RP >= 0
        HalfSpace(sparsevec([RP], [1.], n), α),  # RP < α
        HalfSpace(sparsevec([C2], [-1.], n), 0.),  # C2 >= 0
        HalfSpace(sparsevec([M], [-1.], n), 0.),  # M >= 0
        HalfSpace(sparsevec([T1], [-1.], n), 0.),  # T1 >= 0
        HalfSpace(sparsevec([T2], [-1.], n), 0.),  # T2 >= 0
        HalfSpace(sparsevec([R], [-1.], n), 0.),  # R >= 0
        HalfSpace(sparsevec([Q], [-1.], n), 0.)  # Q >= 0
       ])
    m2 = CBBCS(F, n, invariant)

    # mode 3 ((RP,T2) = (on, off))
    F = diauxic_on_off!
    invariant = HPolyhedron([
        HalfSpace(sparsevec([C1], [-1.], n), -γ),  # C1 >= γ
        HalfSpace(sparsevec([RP], [-1.], n), -α),  # RP >= α
        HalfSpace(sparsevec([C2], [-1.], n), 0.),  # C2 >= 0
        HalfSpace(sparsevec([M], [-1.], n), 0.),  # M >= 0
        HalfSpace(sparsevec([T1], [-1.], n), 0.),  # T1 >= 0
        HalfSpace(sparsevec([T2], [-1.], n), 0.),  # T2 >= 0
        HalfSpace(sparsevec([R], [-1.], n), 0.),  # R >= 0
        HalfSpace(sparsevec([Q], [-1.], n), 0.)  # Q >= 0
       ])
    m3 = CBBCS(F, n, invariant)

    # mode 4 ((RP,T2) = (off, off))
    F = diauxic_off_off!
    invariant = HPolyhedron([
        HalfSpace(sparsevec([C1], [-1.], n), 0.),  # C1 >= 0
        HalfSpace(sparsevec([C1], [1.], n), γ),  # C1 < γ
        HalfSpace(sparsevec([RP], [-1.], n), -α),  # RP >= α
        HalfSpace(sparsevec([C2], [-1.], n), 0.),  # C2 >= 0
        HalfSpace(sparsevec([M], [-1.], n), 0.),  # M >= 0
        HalfSpace(sparsevec([T1], [-1.], n), 0.),  # T1 >= 0
        HalfSpace(sparsevec([T2], [-1.], n), 0.),  # T2 >= 0
        HalfSpace(sparsevec([R], [-1.], n), 0.),  # R >= 0
        HalfSpace(sparsevec([Q], [-1.], n), 0.)  # Q >= 0
       ])
    m4 = CBBCS(F, n, invariant)

    # modes
    modes = [m1, m2, m3, m4]

    # transition 1 -> 2
    add_transition!(automaton, 1, 2, 1)
    guard = HPolyhedron([
        HalfSpace(sparsevec([C1], [-1.], n), -γ),  # C1 >= γ
        HalfSpace(sparsevec([RP], [-1.], n), 0.),  # RP >= 0
        HalfSpace(sparsevec([RP], [1.], n), α)  # RP < α
       ])
    t1 = ConstrainedIdentityMap(n, guard)

    # transition 2 -> 3
    add_transition!(automaton, 2, 3, 2)
    guard = HPolyhedron([
        HalfSpace(sparsevec([C1], [-1.], n), -γ),  # C1 >= γ
        HalfSpace(sparsevec([RP], [-1.], n), -α)  # RP >= α
       ])
    t2 = ConstrainedIdentityMap(n, guard)

    # transition 3 -> 4
    add_transition!(automaton, 3, 4, 3)
    guard = HPolyhedron([
        HalfSpace(sparsevec([C1], [-1.], n), 0.),  # C1 >= 0
        HalfSpace(sparsevec([C1], [1.], n), γ),  # C1 < γ
        HalfSpace(sparsevec([RP], [-1.], n), -α)  # RP >= α
       ])
    t3 = ConstrainedIdentityMap(n, guard)

    # transition 4 -> 1
    add_transition!(automaton, 4, 1, 4)
    guard = HPolyhedron([
        HalfSpace(sparsevec([C1], [-1.], n), 0.),  # C1 >= 0
        HalfSpace(sparsevec([C1], [1.], n), γ),  # C1 < γ
        HalfSpace(sparsevec([RP], [-1.], n), 0.),  # RP >= 0
        HalfSpace(sparsevec([RP], [1.], n), α)  # RP < α
       ])
    t4 = ConstrainedIdentityMap(n, guard)

    # transition 1 -> 4
    add_transition!(automaton, 1, 4, 5)
    guard = HPolyhedron([
        HalfSpace(sparsevec([C1], [-1.], n), 0.),  # C1 >= 0
        HalfSpace(sparsevec([C1], [1.], n), γ),  # C1 < γ
        HalfSpace(sparsevec([RP], [-1.], n), -α)  # RP >= α
       ])
    t5 = ConstrainedIdentityMap(n, guard)

    # transition 4 -> 3
    add_transition!(automaton, 4, 3, 6)
    guard = HPolyhedron([
        HalfSpace(sparsevec([C1], [-1.], n), -γ),  # C1 >= γ
        HalfSpace(sparsevec([RP], [-1.], n), -α)  # RP >= α
       ])
    t6 = ConstrainedIdentityMap(n, guard)

    # transition 3 -> 2
    add_transition!(automaton, 3, 2, 7)
    guard = HPolyhedron([
        HalfSpace(sparsevec([C1], [-1.], n), -γ),  # C1 >= γ
        HalfSpace(sparsevec([RP], [-1.], n), 0.),  # RP >= 0
        HalfSpace(sparsevec([RP], [1.], n), α)  # RP < α
       ])
    t7 = ConstrainedIdentityMap(n, guard)

    # transition 2 -> 1
    add_transition!(automaton, 2, 1, 8)
    guard = HPolyhedron([
        HalfSpace(sparsevec([C1], [-1.], n), 0.),  # C1 >= 0
        HalfSpace(sparsevec([C1], [1.], n), γ),  # C1 < γ
        HalfSpace(sparsevec([RP], [-1.], n), 0.),  # RP >= 0
        HalfSpace(sparsevec([RP], [1.], n), α)  # RP < α
       ])
    t8 = ConstrainedIdentityMap(n, guard)

    # transition 1 -> 3
    add_transition!(automaton, 1, 3, 9)
    guard = HPolyhedron([
        HalfSpace(sparsevec([C1], [-1.], n), -γ),  # C1 >= γ
        HalfSpace(sparsevec([RP], [-1.], n), -α)  # RP >= α
       ])
    t9 = ConstrainedIdentityMap(n, guard)

    # transition 3 -> 1
    add_transition!(automaton, 3, 1, 10)
    guard = HPolyhedron([
        HalfSpace(sparsevec([C1], [-1.], n), 0.),  # C1 >= 0
        HalfSpace(sparsevec([C1], [1.], n), γ),  # C1 < γ
        HalfSpace(sparsevec([RP], [-1.], n), 0.),  # RP >= 0
        HalfSpace(sparsevec([RP], [1.], n), α)  # RP < α
       ])
    t10 = ConstrainedIdentityMap(n, guard)

    # transition 2 -> 4
    add_transition!(automaton, 2, 4, 11)
    guard = HPolyhedron([
        HalfSpace(sparsevec([C1], [-1.], n), 0.),  # C1 >= 0
        HalfSpace(sparsevec([C1], [1.], n), γ),  # C1 < γ
        HalfSpace(sparsevec([RP], [-1.], n), -α)  # RP >= α
       ])
    t11 = ConstrainedIdentityMap(n, guard)

    # transition 4 -> 2
    add_transition!(automaton, 4, 2, 12)
    guard = HPolyhedron([
        HalfSpace(sparsevec([C1], [-1.], n), -γ),  # C1 >= γ
        HalfSpace(sparsevec([RP], [-1.], n), 0.),  # RP >= 0
        HalfSpace(sparsevec([RP], [1.], n), α)  # RP < α
       ])
    t12 = ConstrainedIdentityMap(n, guard)

    # transition annotations
    resetmaps = [t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12]

    # switching
    switchings = [AutonomousSwitching()]

    H = HybridSystem(automaton, modes, resetmaps, switchings)

    # mode 2 is the INITIAL MODE

    return H
end
