using LinearAlgebra

function multi_thermostat(n=2)
    b1 = fill(40.0, n)
    b2 = fill(30.0, n)
    A = Diagonal(fill(-0.5, n))
    ubI = 75
    lbI = 65
    ubG = 74.5
    lbG = 65.5

    @variables x[1:n]

    function thermostat_on()
        invariant = HPolyhedron([HalfSpace(x[i] <= ubI, x) for i in 1:n])
        @system(x' = A*x + b1, x ∈ invariant)
    end

    function thermostat_off()
        invariant = HPolyhedron([HalfSpace(x[i] >= lbI, x) for i in 1:n])
        @system(x' = A*x + b2, x ∈ invariant)
    end

    function thermostat_hybrid()
        automaton = GraphAutomaton(2)
        add_transition!(automaton, 1, 2, 1)
        add_transition!(automaton, 2, 1, 2)

        mode1 = thermostat_on()
        mode2 = thermostat_off()
        modes = [mode1, mode2]

        ## transition on -> off
        guard = HalfSpace(x[1] >= ubG, x)
        trans1 = ConstrainedIdentityMap(1, guard)
        ## transition off -> on
        guard = HalfSpace(x[1] <= lbG, x)
        trans2 = ConstrainedIdentityMap(1, guard)
        resetmaps = [trans1, trans2]

        return HybridSystem(automaton, modes, resetmaps, [AutonomousSwitching()])
    end

    H = thermostat_hybrid()

    return H
end
