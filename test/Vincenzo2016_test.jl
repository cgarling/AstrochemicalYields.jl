using AstrochemicalYields
using Test

v = @test_nowarn Vincenzo2016()
@test v(0.02, 5) ≈ 0.18342380700629182
@test AstrochemicalYields.Lifetimes.derivative(v, 0.02)(5) ≈ -0.06565310290439594
@test inverse(v)(0.02, 0.18342380700629182) ≈ 5
@test AstrochemicalYields.Lifetimes.derivative(inverse, v, 0.02)(5) ≈ -0.07962087433398927