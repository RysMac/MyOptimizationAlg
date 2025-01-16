using MyOptimizationAlg
using LinearAlgebra
using BenchmarkTools

# Positive definite matrix example
function test_positive_definite()
	println("=== Positive Definite Matrix ===")
	hess = [4.0 1.0; 1.0 3.0]    # Positive definite matrix
	grad = [-1.0, -1.0]          # Gradient vector
	delta = 0.01                  # Trust region radius
	sol = subproblem(hess, grad, delta, verbose = 1)
	println("Solution:", sol)
	println("Norm of solution:", norm(sol))
end

# Non-positive definite matrix example
function test_non_positive_definite()
    println("\n=== Non-Positive Definite Matrix ===")
    hess = [0.0 1.0; 1.0 0.0]    # Not positive definite
    grad = [-1.0, -1.0]          # Gradient vector
    delta = 0.01                  # Trust region radius
    @time sol = subproblem(hess, grad, delta, verbose = 2)
    println("Solution:", sol)
    println("Norm of solution:", norm(sol))
end

# Hard case example
function test_hard_case()
    println("\n=== Hard Case Example ===")
    hess = [2.0 0.0; 0.0 -1.0]   # Indefinite matrix (one negative eigenvalue)
    grad = [-2.0, 1.0]           # Gradient vector
    delta = 1.5                  # Trust region radius (solution must be on the boundary)
    @time sol = subproblem(hess, grad, delta, verbose = 2)
    println("Solution:", sol)
    println("Norm of solution:", norm(sol))
end

# Run the tests
function run_tests()
    test_positive_definite()
    test_non_positive_definite()
    test_hard_case()
end
run_tests()
