using Test
using MyOptimizationAlg
import LinearAlgebra

f(x) = x[1]^2 + x[2]^2
grad_f(x) = [2x[1], 2x[2]]

@testset "Gradient Descent Tests" begin
	x0 = [10.0, 10.0]
	result = gradient_descent(f, grad_f, x0, alpha=0.1)
	@test LinearAlgebra.norm(result) ≈ 0.0 atol=1e-4
end
