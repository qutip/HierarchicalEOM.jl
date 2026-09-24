using Test
using HierarchicalEOM
using CUDA
import SparseArrays: spzeros

@testset "CUDA (Basic)" begin
    λ = 0.01
    W = 0.5
    kT = 0.5
    μ = 0
    N = 3
    tier = 3

    # System Hamiltonian
    Hsys = Qobj(spzeros(ComplexF64, 2, 2))

    # system-bath coupling operator
    Qb = sigmax()
    Qf = sigmam()

    e_ops = [sigmax()]

    # initial state
    ψ0 = basis(2, 1)

    Bbath = Boson_DrudeLorentz_Pade(Qb, λ, W, kT, N)
    Fbath = Fermion_Lorentz_Pade(Qf, λ, μ, W, kT, N)

    # Solving time Evolution
    ## Schrodinger HEOMLS
    L_cpu = M_S(Hsys; verbose = false)
    L_gpu = cu(L_cpu)
    L_gpu_csc = CUDA.cuSPARSE.CuSparseMatrixCSC(L_cpu)
    L_gpu_csr = CUDA.cuSPARSE.CuSparseMatrixCSR(L_cpu)
    sol_cpu = heomsolve(L_cpu, ψ0, [0, 10]; e_ops = e_ops, progress_bar = Val(false))
    sol_gpu = heomsolve(L_gpu, ψ0, [0, 10]; e_ops = e_ops, progress_bar = Val(false))
    @test L_gpu.data.A isa CUDA.cuSPARSE.CuSparseMatrixCSC
    @test L_gpu_csc.data.A isa CUDA.cuSPARSE.CuSparseMatrixCSC
    @test L_gpu_csr.data.A isa CUDA.cuSPARSE.CuSparseMatrixCSR
    @test all(isapprox.(sol_cpu.expect[1, :], sol_gpu.expect[1, :], atol = 1.0e-4))
    @test isapprox(getRho(sol_cpu.ados[end]), getRho(sol_gpu.ados[end]), atol = 1.0e-4)

    ## Boson HEOMLS
    L_cpu = M_Boson(Hsys, tier, Bbath; verbose = false)
    L_gpu = cu(L_cpu)
    L_cpu_lazy = M_Boson(Hsys, tier, Bbath; verbose = false, assemble = Val(:combine))
    L_gpu_lazy = cu(L_cpu_lazy)
    sol_cpu = heomsolve(L_cpu, ψ0, [0, 10]; e_ops = e_ops, progress_bar = Val(false))
    sol_gpu = heomsolve(L_gpu, ψ0, [0, 10]; e_ops = e_ops, progress_bar = Val(false))
    sol_gpu_lazy = heomsolve(L_gpu_lazy, ψ0, [0, 10]; e_ops = e_ops, progress_bar = Val(false))
    @test L_gpu.data.A isa CUDA.cuSPARSE.CuSparseMatrixCSC{ComplexF64, Int32}
    @test all(isapprox.(sol_cpu.expect[1, :], sol_gpu.expect[1, :], atol = 1.0e-4))
    @test all(isapprox.(sol_cpu.expect[1, :], sol_gpu_lazy.expect[1, :], atol = 1.0e-4))
    @test isapprox(getRho(sol_cpu.ados[end]), getRho(sol_gpu.ados[end]), atol = 1.0e-4)
    @test isapprox(getRho(sol_cpu.ados[end]), getRho(sol_gpu_lazy.ados[end]), atol = 1.0e-4)

    ## Fermion HEOMLS
    L_cpu = M_Fermion(Hsys, tier, Fbath; verbose = false)
    L_cpu_lazy = M_Fermion(Hsys, tier, Fbath; verbose = false, assemble = Val(:combine))
    L_gpu = cu(L_cpu)
    L_gpu_lazy = cu(L_cpu_lazy)
    sol_cpu = heomsolve(L_cpu, ψ0, [0, 10]; e_ops = e_ops, progress_bar = Val(false))
    sol_gpu = heomsolve(L_gpu, ψ0, [0, 10]; e_ops = e_ops, progress_bar = Val(false))
    sol_gpu_lazy = heomsolve(L_gpu_lazy, ψ0, [0, 10]; e_ops = e_ops, progress_bar = Val(false))
    @test L_gpu.data.A isa CUDA.cuSPARSE.CuSparseMatrixCSC{ComplexF64, Int32}
    @test all(isapprox.(sol_cpu.expect[1, :], sol_gpu.expect[1, :], atol = 1.0e-4))
    @test all(isapprox.(sol_cpu.expect[1, :], sol_gpu_lazy.expect[1, :], atol = 1.0e-4))
    @test isapprox(getRho(sol_cpu.ados[end]), getRho(sol_gpu.ados[end]), atol = 1.0e-4)
    @test isapprox(getRho(sol_cpu.ados[end]), getRho(sol_gpu_lazy.ados[end]), atol = 1.0e-4)

    ## Boson Fermion HEOMLS
    L_cpu = M_Boson_Fermion(Hsys, tier, tier, Bbath, Fbath; verbose = false)
    L_cpu_lazy = M_Boson_Fermion(Hsys, tier, tier, Bbath, Fbath; verbose = false, assemble = Val(:combine))
    L_gpu = cu(L_cpu)
    L_gpu_lazy = cu(L_cpu_lazy)
    tlist = 0:1:10
    sol_cpu = heomsolve(L_cpu, ψ0, tlist; e_ops = e_ops, saveat = tlist, progress_bar = Val(false))
    sol_gpu = heomsolve(L_gpu, ψ0, tlist; e_ops = e_ops, saveat = tlist, progress_bar = Val(false))
    sol_gpu_lazy = heomsolve(L_gpu_lazy, ψ0, tlist; e_ops = e_ops, saveat = tlist, progress_bar = Val(false))
    @test L_gpu.data.A isa CUDA.cuSPARSE.CuSparseMatrixCSC{ComplexF64, Int32}
    @test all(isapprox.(sol_cpu.expect[1, :], sol_gpu.expect[1, :], atol = 1.0e-4))
    @test all(isapprox.(sol_cpu.expect[1, :], sol_gpu_lazy.expect[1, :], atol = 1.0e-4))
    for i in 1:length(tlist)
        @test isapprox(getRho(sol_cpu.ados[i]), getRho(sol_gpu.ados[i]), atol = 1.0e-4)
        @test isapprox(getRho(sol_cpu.ados[i]), getRho(sol_gpu_lazy.ados[i]), atol = 1.0e-4)
    end
end
