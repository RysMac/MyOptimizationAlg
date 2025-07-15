# for now Trust Region uses only subproblem alg
include("subproblem.jl")
include("subproblem_gqtpar.jl")
import .MoreSorensen: gqtpar

"""
    projection!(pg, x, g, lb)

Compute the projected gradient vector for lower bounds.
`pg[i] = P[x[i] - g[i]] - x[i]`,
where P projects onto x_i >= lb[i].
"""
function projection!(pg::AbstractVector,
					x::AbstractVector,
					g::AbstractVector,
					lb::AbstractVector)
    @assert length(pg) == length(x) == length(g) == length(lb)
    @inbounds for i in eachindex(x)
        y = x[i] - g[i]
        yproj = y < lb[i] ? lb[i] : y
        pg[i] = yproj - x[i]
    end
    return nothing
end

"""
    generalized_cauchy(H, g, x, lb, delta; eps=1e-12)

Compute the generalized Cauchy point for the trust-region subproblem
with simple lower bounds. Returns a step vector p_c.
"""
function generalized_cauchy(H::AbstractMatrix, g::AbstractVector,
                            x::AbstractVector, lb::AbstractVector,
                            delta::Float64; eps::Float64=1e-12)
    n = length(x)
    # compute projected gradient
    pg = similar(g)
    projection!(pg, x, g, lb)
    # norm of pg
    pg_norm = norm(pg)
    if pg_norm < eps
        return zeros(n)
    end
    # Hessian-vector product
    Hpg = H * pg
    # optimal step length along pg direction
    alpha = dot(pg, pg) / dot(pg, Hpg)
    # trust-region truncation
    tau = min(alpha, delta / pg_norm)
    # generalized Cauchy point
    return -tau * pg
end

function build_free(x, g, lb; ϵx=1e-10, ϵd=1e-8)
    n = length(x)
    free = Int[]
    for i in 1:n
        if x[i] > lb[i] + ϵx
            # safely off the lower‐bound
            push!(free, i)

        else
            # at the bound: only free if descent wants to move x_i *up*
            # i.e. d_i = -g_i > ϵd
            if -g[i] > ϵd
                push!(free, i)
            end
        end
    end
    return free
end

function build_free_from_p(x, p_full, lb; ϵx=1e-10)
  n = length(x)
  free = Int[]
  for i in 1:n
    # (1) off the bound?  ↦ free
    # (2) or the model step is positive (increase x_i)? ↦ free
    if x[i] > lb[i] + ϵx || p_full[i] > 0
      push!(free, i)
    end
  end
  return free
end

function trust_region_ftb(fun, grad, hess,
                      δ::T, x0::Vector{T};
                      tol = 1e-12,
                      maxit = 100,
                      verbose = 0) where T<:Float64

    #— state variables —#
    x    = copy(x0)                    # current iterate (never mutated until accept)
    f    = fun(x)
    g    = grad(x)
    H    = hess(x)
	lb   = zeros(24)
    # lb   = [-1.0, 1.5]                  # your simple lower bounds
    n    = length(x)

    #— working storage —#
    p    = zeros(T,n)                  # candidate step
    pg   = similar(x)                  # projected gradient

    #— TR parameters —#
    η₁, η₂ = 0.25, 0.75                # acceptance thresholds
    ϵx, ϵg = tol, tol                  # active‐set tolerances
	# free = [1,2]
    for k in 1:maxit
        # 0) check projected‐gradient convergence right away
        projection!(pg, x, g, lb)
        if norm(pg) < tol
			println("Final delta: ", δ)
            verbose>0 && println("→ converged (‖pg‖=$(norm(pg))) in $k iters")
            return x
        end

        # 1) build free set by complementarity
		p_full = generalized_cauchy(H, g, x, lb, delta)
		# println("p_full = ", p_full)

        free = build_free_from_p(x, p_full, lb; ϵx=1e-10)
		println("free set = ", free)

        # 2) if nothing is free, fall back to a simple projected‐gradient step
        if isempty(free)
            p .= -δ * pg / norm(pg)
        else
            # 2a) inner “freeze‐bound‐violators” loop
            fill!(p, 0)                # wipe out last iteration’s step
            while true
                # solve TR subproblem on the free coords
                p_free = subproblem(H[free,free], g[free], δ; tol=1e-4)
				println("subproblem")
                p[free] .= p_free

                # detect which free indices would cross the bound
                viol = [i for i in free if x[i] + p[i] < lb[i] - ϵx]
                if isempty(viol)
                    break                # p is now a feasible step
                end

                # shorten p so it exactly hits the first bound
                α = minimum(
                    (x[i] - lb[i]) / (-p[i])
                    for i in viol if p[i] < 0;
                    init = 1.0
                )
                p .*= α*0.98

                # freeze those violators at the bound, remove from free set
                for i in viol
                    p[i] = lb[i] - x[i]
                end
                free = setdiff(free, viol)
            end
        end

        # 3) compute predicted vs actual reduction
        x_trial = x .+ p
        f_trial = fun(x_trial)
        m_decr  = dot(g,p) + 0.5 * dot(p, H*p)

        # guard against zero model‐decr
        ρ = if abs(m_decr) < eps()
            f_trial < f ? 1.0 : -Inf
        else
            (f_trial - f) / m_decr
        end

        # 4) accept or reject
        if ρ > η₁
            x .= x_trial
			# println("accepted x = ", x)
            f  = f_trial
            g .= grad(x)
            H .= hess(x)

			# println("gradient at sol = ", g)
        end

        # 5) update δ
        if ρ < η₁
            δ *= 0.5
        elseif ρ > η₂ && abs(norm(p) - δ) < 1e-8
            δ = min(2δ, 10.0)
        end

        verbose>1 && println(
            "iter $k │ f=$(round(f,8)) │ ‖pg‖=$(round(norm(pg),3)) │ δ=$(round(δ,3))"
        )
    end
	println("Final delta: ", δ)
	return x
    # error("trust-region failed to converge in $maxit iterations")
end



function trust_region_old(	fun,
						grad,
						hess,
						delta::T,
						x0::Vector{T},
						active::Vector{Int64};
						verbose = 0) where T <: Float64

	sol 		= copy(x0)
	sol_inner 	= similar(x0)
	sol_try		= similar(x0)
	hess_val	= hess(x0)
	grad_val	= grad(x0)
	# println("grad_val in TR = ", grad_val)
	# println("hess_val in TR = ", hess_val)
	func_val	= fun(x0)
	active = collect(eachindex(sol))

	theta = 0.99
	tau = 1.0
	# lower bounds
	lb = [0., 1.5]
	gprj = [0., 0.]
	for i in 1:100
		if verbose > 1
			println("iteration = ", i , "  delta = ", delta)
		end

		# Build free set   (x_i ≤ εb  &&  grad_i > 0)  ⇒  inactive
		ϵx = 1e-10
		ϵg = 1e-10
		active = findall(i ->
			!( sol[i] ≤ lb[i] + ϵx  &&  grad_val[i] ≥ ϵg ),
			1:length(sol)
			)
		println("active set", active)

		hess_free = @view hess_val[active, active]
		grad_free = grad_val[active]
		@inbounds sol_inner[active] .= subproblem(hess_free, grad_free, delta, tol=10^-10)

		# info, sol_inner, iter = gqtpar(hess_val, 'U', grad_val, delta, 10^-8, 10^-8, 100, 0.)[[1,3,5]]
		println("inner sol: ", sol_inner)

		# once subproblem is solved check for crossing boundaries
		tau = 1.0
		for i in eachindex(sol_inner)
			if (sol_inner[i] < 0 && abs(sol_inner[i]) > 1e-16 )
				tau = min(tau, ((sol[i]) - lb[i]) / (-sol_inner[i]))
			end
		end
		sol_inner .*= tau*theta
		println("tau = ", tau)
		println("sol inner = ", sol_inner)
		println("sol inner cut = ", tau .* sol_inner)

		@inbounds sol_try .= sol + sol_inner
		f_try 		= fun(sol_try)
		numerator 	= f_try - func_val
		denominator	= grad_val' * sol_inner + 0.5 * sol_inner' * hess_val * sol_inner

		if abs(numerator) < 10^-16 && abs(numerator - denominator) < 10^-16
			rho = 1
		else
			rho = numerator/denominator
		end
		if verbose > 1
			println("rho = ", rho, " numerator = ", numerator, " denominator = ", denominator)
		end

		if rho > 0.25 # in HSL lib it is 0.01
			@inbounds sol .= sol_try
			println("accepted solution = ", sol)
			func_val = fun(sol)
			grad_val .= grad(sol)
			hess_val .= hess(sol)
			println("grad_val in TR = ", grad_val)
			# println("hess_val in TR = ", hess_val)
		end
		println("TR gradient norm = ", norm(grad_val), "  outer iterations = ", i)
		if rho ≤ 0.5
			delta = delta/4.
		end

		if rho > 0.8 && abs(delta - norm(sol_inner)) ≤ 10^-8 && delta < 10.
			delta *= 2.0
		end



		projection!(gprj, sol, grad(sol), lb)
		println("projected = ", gprj)
		if norm(gprj) < 10^-12 || i == 1000
			println("TR gradient = ", norm(grad_val), "  outer iterations = ", i)
			if verbose > 0
				println("gradient = ", norm(grad_val), "  outer iterations = ", i)
				# println("solution = ", sol, "delta = ", delta)
			end
			break
		end
	end
	return sol, delta
end
