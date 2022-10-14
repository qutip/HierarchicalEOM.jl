using Heom
using Test

@testset "Print version information" begin
    @test Heom.versioninfo() == nothing
end

@testset "Bath and Exponent" begin
    include("bath.jl")
end

@testset "Correlation functions" begin
    include("corr_func.jl")
end

@testset "Heom API" begin 
    include("heom_api.jl")
end

@testset "Time evolution solvers" begin
    include("evolution.jl")
end

@testset "Steady state solvers" begin
    include("steady.jl")
end

@testset "Spectrum" begin
    include("spectrum.jl")
end