fb(m) = 2 * m + 1

δ(j, k) = j == k ? 1 : 0

_fermi(x) = 1.0 / (exp(x) + 1.0)  # Fermi-Dirac distribution

function i_eigval(N, M, p)
    pf = []
    A_matrix = zeros(M, M)
    for j in 1:M
        for k in 1:M
            A_matrix[j, k] = (δ(j, k+1) + δ(j, k-1)) / √((fb(j-1) + p) * (fb(k-1) + p))
        end
    end

    for val in eigvals(A_matrix)[1:N]
        push!(pf, -2 / val)
    end

    return pf
end

function matsubara(N; fermion::Bool)
    ϵ = [0.0]

    if fermion
        for k in 1:N
            push!(ϵ, (2 * k - 1) * π)
        end
    else
        for k in 1:N
            push!(ϵ, 2 * π * k)
        end
    end
    return ϵ
end

# thoss's spectral (N-1/N) pade
function pade_NmN(N; fermion::Bool)

    if fermion
        local ζ = i_eigval(N    , 2 * N    , 0)
        local χ = i_eigval(N - 1, 2 * N - 1, 2)
        prefactor = 0.5 * N * fb(N)
    else
        local ζ = i_eigval(N    , 2 * N    , 2)
        local χ = i_eigval(N - 1, 2 * N - 1, 4)
        prefactor = 0.5 * N * (fb(N) + 2)
    end

    κ = []
    for j in 1:N
        term = prefactor
        for k1 in 1:(N - 1)
            term *= (χ[k1] ^ 2 - ζ[j] ^ 2) / (ζ[k1] ^ 2 - ζ[j] ^ 2 + δ(j, k1))
        end

        term /= (ζ[N] ^ 2 - ζ[j] ^ 2 + δ(j, N))

        push!(κ, term)
    end

    return append!([0.0], κ), append!([0.0], ζ)
end

function _fermi_pade(x, κ, ζ, N)
    f = 0.5
    for l in 2:(N + 1)
        f -= 2 * κ[l] * x / (x ^ 2 + ζ[l] ^ 2)
    end
    return f
end