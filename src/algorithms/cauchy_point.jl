function cauchy_point(delta, gk, Bk)
	"""
	delta: positive scalar, trust-region radius at current iteration
	gk: gradient at current iterate
	Bk: symmetric matrix, e.g., Hessian

	Returns the cauchy point, a vector.
	"""
	tau = 1
	if (gk' * Bk * gk)[1] > 0
		tau = min(norm(gk,2)^3 ./ (delta * gk'*Bk*gk), 1)[1]
	end

	# Cauchy point equation
	return -tau * delta / norm(gk, 2) * gk
end
