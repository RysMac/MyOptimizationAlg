# for now Trust Region uses only subproblem alg
include("subproblem.jl")

function trust_region(	fun,
						grad,
						hess,
						delta::T,
						x0::Vector{T};
						verbose = 0) where T <: Float64

	sol 		= copy(x0)
	sol_inner 	= similar(x0)
	sol_try		= similar(x0)
	hess_val	= hess(x0)
	grad_val	= grad(x0)
	func_val	= fun(x0)
	println("  delta = ", delta)
	for i in 1:1000
		if verbose > 1
			println("iteration = ", i , "  delta = ", delta)
		end
		@inbounds sol_inner .= subproblem(hess_val, grad_val, delta)

		@inbounds sol_try 	.= sol + sol_inner
		f_try 		= fun(sol_try)
		@views numerator 	= f_try - func_val
		@views denominator	= grad_val' * sol_inner + 0.5 * sol_inner' * hess_val * sol_inner

		if abs(numerator) < 10^-16 && abs(numerator - denominator) < 10^-16
			rho = 1
		else
			rho = numerator/denominator
		end
		if verbose > 1
			println("rho = ", rho, " numerator = ", numerator, " denominator = ", denominator)
		end

		if rho > 0.01 # in HSL lib it is 0.01
			@inbounds sol .= sol_try
			#println("accepted solution = ", sol)
			func_val = fun(sol)
			grad_val .= grad(sol)
			hess_val .= hess(sol)
		end

		if rho ≤ 0.5
			delta = delta/4.
		end

		if rho > 0.8 && abs(delta - norm(sol_inner)) ≤ 10^-8 && delta < 10.
			delta *= 2.0
		end
		println("delata = ", delta)
		#println("solution = ", sol)
		if norm(grad_val) < 10^-12 || i == 1000
			if verbose > 0
				println("gradient = ", norm(grad_val), "  outer iterations = ", i)
				# println("solution = ", sol, "delta = ", delta)
			end
			break
		end
	end
	return sol, delta
end
