module NewtonsMethod

using LinearAlgebra, Statistics, Compat, ForwardDiff

function newtonroot(f, f′; x₀, tolerance = 1E-7, maxiter = 1000)
    # setup the algorithm
    x_old = x₀
    normdiff = Inf
    iter = 1
    while normdiff > tolerance && iter <= maxiter
        x_new = x_old - f(x_old)/f′(x_old) # use the passed in map
        normdiff = norm(x_new - x_old)
        x_old = x_new
        iter = iter + 1
    end

    if iter < maxiter+1
        return (root = x_old, normdiff = normdiff, iter = iter) # A named tuple
    else
        return nothing
    end
end

# for auto-differentiation, can just forward to the existing one
D(f) = x -> ForwardDiff.derivative(f, x)

# same name! different parameters
newtonroot(f; x₀, tolerance = 1E-7, maxiter = 1000) = newtonroot(f, D(f), x₀=x₀, tolerance = tolerance, maxiter = maxiter)

export newtonroot
end # module
