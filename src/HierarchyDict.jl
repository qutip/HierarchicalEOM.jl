# generate index to ado vector
function _IDX2ADO(dims::Vector{Int}, N_exc::Int)
    len = length(dims)
    state = zeros(Int, len)
    result = [copy(state)]
    nexc = 0

    while true
        idx = len
        state[end] += 1
        nexc += 1
        if state[idx] < dims[idx]
            push!(result, copy(state))
        end
        while (nexc == N_exc) || (state[idx] == dims[idx])
            #state[idx] = 0
            idx -= 1
            if idx < 1
                return result
            end

            nexc -= state[idx + 1] - 1
            state[idx + 1] = 0
            state[idx] += 1
            if state[idx] < dims[idx]
                push!(result, copy(state))
            end
        end
    end
end

function ADOs_dictionary(dims::Vector{Int}, N_exc::Int)
    ado2idx = OrderedDict{Vector{Int}, Int}()
    idx2ado = _IDX2ADO(dims, N_exc)
    for (idx, ado) in enumerate(idx2ado)
        ado2idx[ado] = idx
    end

    return length(idx2ado), ado2idx, idx2ado
end