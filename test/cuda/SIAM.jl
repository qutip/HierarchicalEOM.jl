using Test
using HierarchicalEOM
using CUDA
import LinearSolve: KrylovJL_BICGSTAB

@testset "CUDA (Single impurity Anderson model)" begin
    ϵ = -5
    U = 10
    σm = sigmam() ## σ-
    σz = sigmaz() ## σz
    II = qeye(2)  ## identity matrix
    d_up = tensor(σm, II)
    d_dn = tensor(-1 * σz, σm)
    ψ0 = tensor(basis(2, 0), basis(2, 0))
    Hsys = ϵ * (d_up' * d_up + d_dn' * d_dn) + U * (d_up' * d_up * d_dn' * d_dn)
    Γ = 2
    μ = 0
    W = 10
    kT = 0.5
    N = 5
    tier = 3
    bath_up = Fermion_Lorentz_Pade(d_up, Γ, μ, W, kT, N)
    bath_dn = Fermion_Lorentz_Pade(d_dn, Γ, μ, W, kT, N)
    bath_list = [bath_up, bath_dn]

    ## solve stationary state
    L_even_cpu = M_Fermion(Hsys, tier, bath_list; verbose = false)
    L_even_cpu_lazy = M_Fermion(Hsys, tier, bath_list; verbose = false, assemble = Val(:combine))
    L_even_gpu = cu(L_even_cpu)
    L_even_gpu_lazy = cu(L_even_cpu_lazy)
    ados_cpu = steadystate(L_even_cpu; verbose = false)
    ados_gpu1 = steadystate(L_even_gpu; verbose = false)
    ados_gpu2 = steadystate(CUDA.cuSPARSE.CuSparseMatrixCSR(L_even_cpu); verbose = false)
    ados_gpu3 = steadystate(L_even_gpu, ψ0, 10; verbose = false)
    ados_gpu_lazy = steadystate(L_even_gpu_lazy; verbose = false)
    @test L_even_gpu.data.A isa CUDA.cuSPARSE.CuSparseMatrixCSC{ComplexF64, Int32}
    @test all(isapprox.(ados_cpu.data, ados_gpu1.data; atol = 1.0e-6))
    @test all(isapprox.(ados_cpu.data, ados_gpu2.data; atol = 1.0e-6))
    @test all(isapprox.(ados_cpu.data, ados_gpu3.data; atol = 1.0e-6))
    @test all(isapprox.(ados_cpu.data, ados_gpu_lazy.data; atol = 1.0e-6))

    ## solve density of states
    ωlist = -5:0.5:5
    L_odd_cpu = M_Fermion(Hsys, tier, bath_list, ODD; verbose = false)
    L_odd_cpu_lazy = M_Fermion(Hsys, tier, bath_list, ODD; verbose = false, assemble = Val(:combine))
    L_odd_gpu_32 = cu(L_odd_cpu, word_size = Val(32))
    L_odd_gpu_64 = cu(L_odd_cpu, word_size = Val(64))
    L_odd_gpu_32_lazy = cu(L_odd_cpu_lazy, word_size = Val(32))
    L_odd_gpu_64_lazy = cu(L_odd_cpu_lazy, word_size = Val(64))
    dos_cpu = DensityOfStates(L_odd_cpu, ados_cpu, d_up, ωlist; progress_bar = Val(false))
    dos_gpu_32 = DensityOfStates(
        L_odd_gpu_32,
        ados_cpu,
        d_up,
        ωlist;
        progress_bar = Val(false),
        alg = KrylovJL_BICGSTAB(rtol = 1.0f-12, atol = 1.0f-14), # somehow KrylovJL_GMRES doesn't work for Float32 (it takes forever to solve)
    )
    dos_gpu_64 = DensityOfStates(L_odd_gpu_64, ados_cpu, d_up, ωlist; progress_bar = Val(false))
    dos_gpu_32_lazy = DensityOfStates(
        L_odd_gpu_32_lazy,
        ados_cpu,
        d_up,
        ωlist;
        progress_bar = Val(false),
        alg = KrylovJL_BICGSTAB(rtol = 1.0f-12, atol = 1.0f-14), # somehow KrylovJL_GMRES doesn't work for Float32 (it takes forever to solve)
    )
    dos_gpu_64_lazy = DensityOfStates(L_odd_gpu_64_lazy, ados_cpu, d_up, ωlist; progress_bar = Val(false))
    @test L_odd_gpu_32.data.A isa CUDA.cuSPARSE.CuSparseMatrixCSC{ComplexF32, Int32}
    @test L_odd_gpu_64.data.A isa CUDA.cuSPARSE.CuSparseMatrixCSC{ComplexF64, Int32}
    for (i, ω) in enumerate(ωlist)
        @test dos_cpu[i] ≈ dos_gpu_32[i] atol = 1.0e-6
        @test dos_cpu[i] ≈ dos_gpu_64[i] atol = 1.0e-6
        @test dos_cpu[i] ≈ dos_gpu_32_lazy[i] atol = 1.0e-6
        @test dos_cpu[i] ≈ dos_gpu_64_lazy[i] atol = 1.0e-6
    end
end
