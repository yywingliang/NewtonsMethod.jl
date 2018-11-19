using NewtonsMethod, LinearAlgebra
using Test

@testset "NewtonsMethod.jl" begin

#Test with analytical derivates
f(x)=x^2
f′(x)=2x
@test newtonroot(f,f′,x₀=0.5).root ≈ 0.0 atol=0.0000001

f(x)=x^3-3x^2+2
f′(x)=3x^2-6x
@test newtonroot(f,f′,x₀=0.5).root ≈ 1.0 atol=0.0000001

f(x)=log(x)-1
f′(x)=1/x
@test newtonroot(f,f′,x₀=0.5).root ≈ 2.71828182845904523536


#Test with automatic derivatives
f(x)=x^2
@test newtonroot(f,x₀=0.5).root ≈ 0.0 atol=0.0000001

f(x)=x^3-3x^2+2
@test newtonroot(f,x₀=0.5).root ≈ 1.0 atol=0.0000001

f(x)=log(x)-1
@test newtonroot(f,x₀=0.5).root ≈ 2.71828182845904523536


#Test Big Float
f(x)=sin(x)-1
f′(x)=cos(x)
@test newtonroot(f, f′, x₀ = BigFloat(1.0), tolerance = 1E-32).root ≈ BigFloat(pi/2) atol = 1E-16
@test newtonroot(f, x₀ = BigFloat(1.0), tolerance = 1E-32).root ≈ BigFloat(pi/2) atol = 1E-16


#Test non-convergence
#f(x)=x^2+2
#@test newtonroot(f,x₀=0.5).root == nothing


#Test maxiter
#f(x)=log(x)-1
#@test newtonroot(f,x₀=0.5,maxiter=5).root == nothing


#Test tolerance
f(x)=log(x)-1
@test newtonroot(f,x₀=0.5).root ≈ 2.71828182845904523536
r1=newtonroot(f,x₀=1.0).root
r2=newtonroot(f,x₀=1.0, tolerance=0.00001).root
r3=newtonroot(f,x₀=1.0, tolerance=0.001).root
@test norm(f(r1)) <= norm(f(r2)) <= norm(f(r3)) # accuracy decreases as tolerance increases


end
