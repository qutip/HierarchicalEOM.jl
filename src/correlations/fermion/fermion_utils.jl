function f_approx(x, κ, ϵ, N)
    f = 0.5
    for l in 2:(N + 1)
        f -= 2 * κ[l] * x / (x ^ 2 + ϵ[l] ^ 2)
    end
    return f
end