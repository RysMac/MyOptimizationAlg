# for now Trust Region uses only subproblem alg
include("trust_region.jl")
using LinearAlgebra

# Lagranfian for penaltyzing out of bounds values
function la_lagrangian(x::AbstractVector{T}, λ::AbstractVector{T}, LowerBounds::AbstractVector{T}, penalty::T) where T <: Float64
	total = zero(T)  # Initialize total to zero
	@inbounds for i in eachindex(x)
		diff = x[i] - LowerBounds[i]
		term = λ[i] + penalty * diff
		if term < 0
			total += λ[i] * diff + 0.5 * penalty * diff^2
		end
	end
	return total
end

function la_grad(x::AbstractVector{T}, λ::AbstractVector{T}, LowerBounds::AbstractVector{T}, penalty::T) where T <: Float64
	# grad should be in-place modified?
	grad = similar(x)  # Preallocate gradient array with the same type as x
	@inbounds for i in eachindex(x)
		diff = x[i] - LowerBounds[i]
		term = λ[i] + penalty * diff
		grad[i] = term < 0 ? term : 0
	end
	return grad
end

function la_hess(x::AbstractVector{T}, λ::AbstractVector{T}, LowerBounds::AbstractVector{T}, penalty::T) where T <: Float64
	hess = zeros(T, length(x))  # Preallocate diagonal Hessian (vector form)
	@inbounds for i in eachindex(x)
		diff = x[i] - LowerBounds[i]
		term = λ[i] + penalty * diff
		hess[i] = term < 0 ? penalty : 0
	end
	return Diagonal(hess)  # Return as a Diagonal matrix
end

# Update Lagrangea multipliers
function la_update!(λ, penalty, variables, LowerBounds)
	for i in eachindex(λ, variables, LowerBounds)
		if λ[i] + penalty * (variables[i] - LowerBounds[i]) < 0
			λ[i] += penalty * (variables[i] - LowerBounds[i])
		else
			λ[i] = 0
		end
	end
	return λ
end

# Measure of constraints violation
function constraints_violation(x::AbstractVector{T}, LowerBounds::AbstractVector{T}) where T <: Float64
	result = similar(x)  # Preallocate output array with the same type as `x`
	@inbounds for i in eachindex(x, LowerBounds)
		result[i] = x[i] <= LowerBounds[i] ? x[i] - LowerBounds[i] : 0
	end
	return result
end

function trust_region_LA(	fun,
							grad,
							hess,
							delta::T,
							x::Vector{T},
							la::Vector{T},
							penalty::T,
							lower_bounds::Vector{T};
							verbose = 0) where T <: Float64




	sol = copy(x)
	λ 	= copy(la) #la_update!(la, penalty, x, lower_bounds)
	w1 	= 1 / (penalty^0.1)
	w2 	= 1 / penalty;

	for i in 1:100

		fun_LA(z)	= fun(z) + la_lagrangian(z, λ, lower_bounds, penalty)
		grad_LA(z) 	= grad(z) + la_grad(z, λ, lower_bounds, penalty)
		hess_LA(z) 	= hess(z) + la_hess(z, λ, lower_bounds, penalty)

		sol .= trust_region(fun_LA, grad_LA, hess_LA, delta, x, verbose = 1)
		# if verbose > 0
		# 	println("trust region sol = ", sol)
		# end

		if norm(constraints_violation(sol, lower_bounds)) ≤ w1
			# penalty = penalty
			λ .= la_update!(λ, penalty, sol, lower_bounds)
			# if verbose > 0
			# 	println("multipliers updated = ", λ)
			# end
			w1 = w1 / (penalty^0.8)
			w2 = w2 / penalty
		else
			penalty = penalty * 2
			# if verbose > 0
			# 	println("penalty updated = ", penalty)
			# end
			w1 = 1 / (penalty^0.1)
			w2 = 1 / penalty

		end
	end
	return sol
end
