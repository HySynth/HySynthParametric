function single_thermostat()
    a1 = 40.0
    a2 = 30.0
    b = 0.5
    ubI = 75
    lbI = 65
    ubG = 74.5
    lbG = 65.5

    @variables x

    function thermostat_on()
        invariant = HalfSpace(x <= ubI)
        @system(x' = -b*x + a1, x ∈ invariant)
    end

    function thermostat_off()
        invariant = HalfSpace(x >= lbI)
        @system(x' = -b*x + a2, x ∈ invariant)
    end

    function thermostat_hybrid()
        automaton = GraphAutomaton(2)
        add_transition!(automaton, 1, 2, 1)
        add_transition!(automaton, 2, 1, 2)

        mode1 = thermostat_on()
        mode2 = thermostat_off()
        modes = [mode1, mode2]

        ## transition on -> off
        guard = HalfSpace(x >= ubG)
        trans1 = ConstrainedIdentityMap(1, guard)
        ## transition off -> on
        guard = HalfSpace(x <= lbG)
        trans2 = ConstrainedIdentityMap(1, guard)
        resetmaps = [trans1, trans2]

        return HybridSystem(automaton, modes, resetmaps, [AutonomousSwitching()])
    end

    H = thermostat_hybrid()

    return H
end
