fb(m) = 2 * m + 1

function δ(j, k)
    if j == k
        return 1
    else
        return 0
    end
end

function i_eigval(r, p, s)
    pf = []
    A_matrix = zeros(r, r)
    for j in 1:r
        for k in 1:r
            A_matrix[j, k] = (δ(j, k+1) + δ(j, k-1)) / √((fb(j-1) + p) * (fb(k-1) + p))
        end
    end

    for val in eigvals(A_matrix)[1:s]
        push!(pf, -2 / val)
    end

    return pf
end

#thoss's spectral (N-1/N) pade
function pade_NmN(lmax,fermion::Bool=true)

    if fermion
        local ϵ = i_eigval(2 * lmax    , 0, lmax)
        local χ = i_eigval(2 * lmax - 1, 2, lmax - 1)
        prefactor = 0.5 * lmax * fb(lmax)
    else
        local ϵ = i_eigval(2 * lmax    , 2, lmax)
        local χ = i_eigval(2 * lmax - 1, 4, lmax - 1)
        prefactor = 0.5 * lmax * (fb(lmax) + 2)
    end

    η_list = []
    for j in 1:lmax
        term = prefactor
        for k1 in 1:(lmax - 1)
            term *= (χ[k1] ^ 2 - ϵ[j] ^ 2) / (ϵ[k1] ^ 2 - ϵ[j] ^ 2 + δ(j, k1))
        end

        term /= (ϵ[lmax] ^ 2 - ϵ[j] ^ 2 + δ(j, lmax))

        push!(η_list, term)
    end

    return append!([0.0], η_list), append!([0.0], ϵ)
end

function f_approx(x, η, ϵ, lmax)
    f = 0.5
    for l in 2:(lmax + 1)
        f -= 2 * η[l] * x / (x ^ 2 + ϵ[l] ^ 2)
    end
    return f
end 

# Correlation function
function Correlation(σ, μ, η, ϵ, β, W, Γ, lmax)
    local η_list = [0.5 * Γ * W * f_approx(1.0im * β * W, η, ϵ, lmax)]
    local γ_list = [W - σ * 1.0im * μ]

    if lmax > 0
        for l in 2:(lmax + 1)
            append!(η_list, -1.0im * (η[l] / β) * Γ * W ^ 2 / (-(ϵ[l] / β) ^ 2 + W ^ 2))
            append!(γ_list, ϵ[l] / β - σ * 1.0im * μ)
        end
    end
    return η_list, γ_list
end