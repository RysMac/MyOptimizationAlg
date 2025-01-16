module MyOptimizationAlg

# Include algorithm files
include("algorithms/gradient_descent.jl")
include("algorithms/subproblem.jl")

# Export functions for external use
export gradient_descent, subproblem

end # module MyOptimizationAlg
