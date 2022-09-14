"""
compute the polyhedron satisfying the tube constraints (constant dynamics)

# Inputs

- ts: the time series
- switches: the indices when a mode switch happens
- mode_sequence: the mode sequence
- p: the total number of modes (optional)
- idx: the index in a list of time series (determines the offset for x0)
- total: the total length of the list of time series (determines the offset for x0)
- ε_mode: mode for ε parameter(s)

# Output

A polyhedron with n dimensions and m constraints such that:

- n = p * dim + total * dim + #ε
- The first p * dim dimensions are the mode slopes.
- The next dim dimensions are the coordinates of x0.
- The last dimensions are reserved for ε.

- m = 2 * dim * k + 2 * dim
- There are 2 * dim * k constraints to bound the slope parameters.
- There are 2 * dim constraints to bound x0.
- The number of dimensions for ε depends on the ε_mode.

Here k is the number of pieces.

# Features

- variable x0
- variable ε
- n dimensions
- p modes
"""
function parameter_range(ts::TimeSeries, switches, mode_sequence, n_without_ε,
                         n_total; p=maximum(mode_sequence), idx::Int=1,
                         ε_mode=SingleEpsilonMode(), ε_offset::Int=0)
    @assert length(switches) == length(mode_sequence) - 1 "incompatible inputs: " *
        "$(length(switches)) switches and $(length(mode_sequence)) modes"

    d = ts.data
    t = ts.time
    k = length(d)-1  # number of pieces
    dim = length(d[1])  # problem dimension
    m = 2 * dim * k + 2 * dim  # number of constraints
    x0_offset = (idx - 1) * dim
    sw_idx = 1  # index of the next switch
    row_jump = 2 * dim  # amount by which one iteration jumps
    col_ε = _get_col_ε(ε_mode, n_without_ε, ε_offset)
    col_slope(mode, j) = dim * (mode - 1) + j  # column of slope constraint in dimension j
    col_x0(j) = p * dim + x0_offset + j  # column of x0 constraint in dimension j
    column_copy, copy_ε_entries = _get_column_copy(ε_mode, n_without_ε)

    Δt = [t[i+1] - t[i] for i in 1:k]  # durations of each piece

    mode = mode_sequence[sw_idx]  # starting mode

    # define system of linear constraints
    A = zeros(m, n_total)
    b = zeros(m)

    # constraints for slope parameters
    row = 0
    for i in 1:k  # i-th piece
        next_d = d[i+1]
        row_tmp = row
        if i > 1
            # copy previous lhs
            for j in 1:row_jump
                row = row_tmp + j
                A[row, column_copy] = A[row - row_jump, column_copy]
            end
        end
        row = row_tmp
        for j in 1:dim
            # x0_j + m_ij * Δt_i - ε_i <= d_ij
            row += 1
            if i == 1
                # only set the first time (copied later)
                A[row, col_x0(j)] = 1.0  # x0_j
            end
            if !copy_ε_entries || i == 1
                A[row, col_ε(i)] = -1.0  # -ε_i
            end
            A[row, col_slope(mode, j)] += Δt[i]
            b[row] = next_d[j]

            # -x0_j - m_ij * Δt_i - ε_i <= -d_ij
            row += 1
            if i == 1
                # only set the first time (copied later)
                A[row, col_x0(j)] = -1.0  # -x0_j
            end
            if !copy_ε_entries || i == 1
                A[row, col_ε(i)] = -1.0  # -ε_i
            end
            A[row, col_slope(mode, j)] -= Δt[i]
            b[row] = -next_d[j]
        end

        if sw_idx <= length(switches) && switches[sw_idx] == i + 1
            sw_idx += 1
            mode = mode_sequence[sw_idx]
        end
    end
    # constraints for x0
    d0 = d[1]
    for j in 1:dim
        # x0_j - ε_0 <= d_0j
        row += 1
        A[row, col_x0(j)] = 1.0  # x0_j
        A[row, col_ε(0)] = -1.0  # -ε_0
        b[row] = d0[j]  # d_ij

        # -x0_j - ε_0 <= -d_0j
        row += 1
        A[row, col_x0(j)] = -1.0  # -x0_j
        A[row, col_ε(0)] = -1.0  # -ε_0
        b[row] = -d0[j]  # -d_ij
    end
    return HPolyhedron(A, b)  # polyhedron A*x <= b
end

function _get_col_ε(ε_mode::SingleEpsilonMode, n_without_ε::Int, ε_offset::Int)
    return i -> n_without_ε + 1
end

function _get_col_ε(ε_mode::MultiEpsilonMode, n_without_ε::Int, ε_offset::Int)
    return i -> n_without_ε + ε_offset + i + 1
end

function _get_column_copy(ε_mode::SingleEpsilonMode, n_without_ε::Int)
    column_copy = 1:(n_without_ε + 1)
    copy_ε_entries = true
    return (column_copy, copy_ε_entries)
end

function _get_column_copy(ε_mode::MultiEpsilonMode, n_without_ε::Int)
    column_copy = 1:n_without_ε
    copy_ε_entries = false
    return (column_copy, copy_ε_entries)
end
