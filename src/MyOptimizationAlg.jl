module MyOptimizationAlg

# Include algorithm files
include("algorithms/gradient_descent.jl")
include("algorithms/cauchy_point.jl")
include("algorithms/subproblem.jl")
include("algorithms/trust_region.jl")
include("algorithms/trust_region_LA.jl")

# Export functions for external use
export gradient_descent, subproblem

end # module MyOptimizationAlg
