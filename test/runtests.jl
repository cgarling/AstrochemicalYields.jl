using AstrochemicalYields
using Test
using SafeTestsets: @safetestset

@testset "AstrochemicalYields.jl" begin
    # Write your tests here.
    @safetestset "Vincenzo2016 Stellar Lifetimes" include("Vincenzo2016_test.jl")
end
