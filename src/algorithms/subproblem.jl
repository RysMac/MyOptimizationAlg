using LinearAlgebra
using Roots

function subproblem(hess::AbstractMatrix{T}, grad::Vector{T}, delta::T; tol::T = 1e-4, max_iters::Int = 100, verbose::Int = 0) where T <: Float64

	n = length(grad)
	sol = zeros(T, n)
	lambda = zero(T)

	if !LinearAlgebra.isposdef(hess)
		eigen_sys = eigen(hess)
		eigen_val_min = minimum(real(eigen_sys.values))
		if verbose > 1
			println("not positive definite: eigenvalue min = ", eigen_val_min)
		end
		lambda = abs(eigen_val_min) + 0.0001
	end

	for i in 1:max_iters

		# Regularize Hessian if needed
		regularized_hess = 0.5(hess' + hess) + lambda * I
		regularized_hess_ch = LinearAlgebra.cholesky(regularized_hess).U

		# Newton step
		RTR = regularized_hess_ch' * regularized_hess_ch
		sol_newton = RTR \ -grad
		sol_newton_norm = norm(sol_newton)
		#println("sol Newton = ", mod_hess_ch)

		if sol_newton_norm ≤ delta

			if lambda == 0 || abs(sol_newton_norm - delta ) < tol * delta

				# since it is the solution just return ??
				sol .= sol_newton
				sol_norm = sol_newton_norm
				if verbose > 0
					println("Newton step ")
				end
				break

			else
				eigen_vec_min = (eigen(hess, sortby=minimum).vectors)[:,1]
				fn = x -> sqrt((sol_newton + x * eigen_vec_min)' * (sol_newton + x * eigen_vec_min)) - delta
				root = find_zero(fn, delta/100.)
				sol .= sol_newton + root * eigen_vec_min
				sol_norm = norm(sol)
				if verbose > 1
					prinln("hard case solution")
				end
				# also use return ?
				break

			end

		else
			sol .= sol_newton
			sol_norm = sol_newton_norm
			ql_norm = norm(regularized_hess_ch' \ sol)
			lambda = lambda + (sol_norm/ql_norm)^2 * ((sol_norm - delta)/delta);
			#println("inner error = ", sol_norm - delta )
		end

		if abs(sol_norm - delta ) < tol #* delta
			if (verbose > 0)
				println("sub problem solution was found after i = ", i , " iterations", "  sol_nomr = ", sol_norm)
			end
			break
		end

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
