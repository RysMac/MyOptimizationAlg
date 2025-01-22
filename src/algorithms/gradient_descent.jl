function	gradient_descent(f, grad_f, x0; max_iters=100, alpha=0.01)
	x = x0
	for i in 1:max_iters
		x -= alpha * grad_f(x)
	end
	return x
end
