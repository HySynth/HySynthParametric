# modifying version
function merge_same_slopes!(switches, mode_sequence)
    @assert length(switches) == length(mode_sequence) - 1 "incompatible inputs: " *
        "$(length(switches)) switches and $(length(mode_sequence)) modes"
    todelete = Int[]
    for i in eachindex(switches)
        if mode_sequence[i] == mode_sequence[i+1]
            push!(todelete, i)
        end
    end
    deleteat!(switches, todelete)
    deleteat!(mode_sequence, todelete)
    return switches, mode_sequence
end

# non-modifying version (old, not used)
function merge_same_slopes(switches, mode_sequence)
    @assert length(switches) == length(mode_sequence) - 1 "incompatible inputs: " *
        "$(length(switches)) switches and $(length(mode_sequence)) modes"
    switches_new = Int[]
    mode_sequence_new = Int[]
    for i in eachindex(switches)
        if mode_sequence[i] != mode_sequence[i+1]
            push!(switches_new, switches[i])
            push!(mode_sequence_new, mode_sequence[i])
        end
    end
    push!(mode_sequence_new, mode_sequence[end])
    return switches_new, mode_sequence_new
end
