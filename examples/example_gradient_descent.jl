using MyOptimizationAlg

# Define a simple quadratic function
f(x) = x[1]^2 + x[2]^2
grad_f(x) = [2x[1], 2x[2]]

# Test gradient descent
x0 = [10.0, 8.0]
result = gradient_descent(f, grad_f, x0, max_iters=1000, alpha=0.01)
println("Result using gradient descent: ", result)
