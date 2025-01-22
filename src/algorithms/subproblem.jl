using LinearAlgebra
using Roots

function subproblem(hess::AbstractMatrix{T}, grad::Vector{T}, delta::T; tol::T = 1e-4, max_iters::Int = 100, verbose::Int = 2) where T <: Float64

	n 					= length(grad)
	sol 				= zeros(T, n)
	regularized_hess 	= similar(hess)
	lambda 				= zero(T)

	if !LinearAlgebra.isposdef(hess)
		eigen_sys = eigen(hess)
		eigen_val_min = minimum(real(eigen_sys.values))
		lambda = abs(eigen_val_min) + 0.00001
		if verbose > 1
			println("not positive definite: eigenvalue min = ", eigen_val_min)
		end
	end

	for i in 1:max_iters

		# Regularize Hessian if needed
		regularized_hess .= 0.5(hess' + hess) + lambda * I

		# Newton step
		#RTR = regularized_hess_ch' * regularized_hess_ch
		sol .= regularized_hess \ -grad
		println("newton step = ", sol);
		println("Lambda = ", lambda)
		sol_norm = norm(sol)
		println("sol_norm = ", sol_norm)

		if sol_norm ≤ delta
			if lambda == 0 || abs(sol_norm - delta ) < tol * delta
				println("Newton step finished")
				if verbose > 0
					println("Newton step ")
				end
				return sol
			else
				eigen_vec_min = (eigen(hess, sortby=minimum).vectors)[:,1]
				fn = x -> sqrt((sol + x * eigen_vec_min)' * (sol + x * eigen_vec_min)) - delta
				root = find_zero(fn, delta/100.)
				if verbose > 1
					prinln("hard case solution")
				end
				return sol + root * eigen_vec_min
			end
		else
			ql_norm = norm((LinearAlgebra.cholesky(regularized_hess).U)' \ sol)
			lambda += (sol_norm/ql_norm)^2 * ((sol_norm - delta)/delta);
			#println("inner error = ", sol_norm - delta )
		end

		# if abs(sol_norm - delta ) < tol * delta
		# 	if (verbose > 0)
		# 		println("sub problem solution was found after i = ", i , " iterations", "  sol_nomr = ", sol_norm)
		# 	end
		# 	break
		# end

		if i == max_iters
			println("no solution in inner loop !!!!")
		end

	end
	return sol
end





# function subproblem(hess::AbstractMatrix{T}, grad::Vector{T}, delta::T; tol::T = 1e-4, max_iters::Int = 100) where T <: Real
#     # Initialization
#     n = length(grad)
#     sol = zeros(T, n)
#     λ = zero(T)

#     # Ensure Hessian is positive definite
#     if !isposdef(hess)
#         eigen_sys = eigen(hess)
#         eigen_val_min = minimum(real(eigen_sys.values))
#         λ = abs(eigen_val_min) + 1e-4  # Regularization factor
#     end

#     # Iterative refinement
#     for iter in 1:max_iters
#         # Regularize Hessian
#         regularized_hess = hess + λ * I

#         # Solve Newton step
#         sol_newton = regularized_hess \ -grad
#         sol_norm = norm(sol_newton)

#         if sol_norm ≤ delta
#             # Check if the Newton step is valid or needs refinement
#             if λ == 0 || abs(sol_norm - delta) < tol * delta
#                 return sol_newton
#             else
#                 # Hard case: Find root to project onto the boundary
#                 eigen_sys = eigen(hess, sortby=minimum)
#                 eigen_vec_min = eigen_sys.vectors[:, 1]
#                 fn = x -> norm(sol_newton + x * eigen_vec_min) - delta
#                 root = find_zero(fn, delta / 100)
#                 return sol_newton + root * eigen_vec_min
#             end
#         else
#             # Update λ to refine the regularization
#             ql_norm = norm(cholesky(regularized_hess).U' \ sol_newton)
#             λ += (sol_norm / ql_norm)^2 * ((sol_norm - delta) / delta)
#         end

#         # Convergence check
#         if abs(sol_norm - delta) < tol * delta
#             break
#         end

#         if iter == max_iters
#             println("Warning: Subproblem did not converge within max iterations!")
#         end
#     end

#     return sol
# end
