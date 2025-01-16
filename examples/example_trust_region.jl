# Include the Trust Region algorithm

include("../src/algorithms/trust_region.jl")
include("../src/algorithms/cauchy_point.jl")

using BenchmarkTools
# Define the Rosenbrock function
function rosenbrock(x::Vector{T}) where T
    return (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
end

# Define the gradient of the Rosenbrock function
function gradient_rosenbrock(x::Vector{T}) where T
	grad = similar(x)
	grad[1] = -2 * (1 - x[1]) - 400 * x[1] * (x[2] - x[1]^2)
	grad[2] = 200 * (x[2] - x[1]^2)
	return grad
end

# Define the Hessian of the Rosenbrock function
function hessian_rosenbrock(x::Vector{T}) where T
	hess = Matrix{T}(undef, 2, 2)
	hess[1, 1] = 1200 * x[1]^2 - 400 * x[2] + 2
	hess[1, 2] = -400 * x[1]
	hess[2, 1] = -400 * x[1]
	hess[2, 2] = 200
	return hess
end

# Initial guess
x0 = [3.2, 2.2]  # Close to the global minimum but not at it

# Trust region radius
delta = .5

x0 = cauchy_point(delta, gradient_rosenbrock(x0), hessian_rosenbrock(x0))

gradient_rosenbrock(x0)

# Run the trust-region algorithm
solution = trust_region(rosenbrock, gradient_rosenbrock, hessian_rosenbrock, delta, x0, verbose = 1)

# Print the solution
println("Optimized solution: $solution")
println("Objective value at solution: ", rosenbrock(solution))

# 44.944 μs (1171 allocations: 55.66 KiB)
