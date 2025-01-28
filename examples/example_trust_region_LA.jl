
using ForwardDiff
using LinearAlgebra
using StaticArrays
include("../src/algorithms/cauchy_point.jl")
include("../src/algorithms/trust_region_LA.jl")


function inc_energy(γ::Vector{Float64})
	return γ' * γ
end

function inc_energy_grad(γ::Vector{Float64})
    return 2 * γ #ForwardDiff.gradient(inc_energy, γ)
end

function inc_energy_hess(γ::Vector{Float64})
    return [2 0; 0 2] #ForwardDiff.hessian(inc_energy, γ)
end



dγ = [5., 3.] .* 4

delta = 1.;
penalty = 100.;
lower_bounds = [5., 2.]

# inc_energy(dγ)
# inc_energy_grad(dγ)
# inc_energy_hess(dγ)

dγ += cauchy_point(delta, inc_energy_grad(dγ), inc_energy_hess(dγ) )

@time for i in 1:1
	γsol = trust_region_LA(inc_energy, inc_energy_grad, inc_energy_hess, delta, dγ, penalty, lower_bounds; verbose = 1)
end


println(inc_energy_grad(γsol))

println(γsol)

count_positives_negatives(γsol)

