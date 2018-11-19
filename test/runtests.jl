using NewtonsMethod
using Test

@testset "NewtonsMethod.jl" begin

#Test 1.0

f1(x) = (x-1.0)^3
 f1′(x) = 3.0*(x-1.0)^2
 x1₀ = 0.1
 @test newtonroot(f1, f1′, x₀ = x1₀).root ≈ 0.9999998168636032
 @test newtonroot(f1,      x₀ = x1₀).root ≈ 0.9999998168636032


end
