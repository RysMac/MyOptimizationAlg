using LinearAlgebra
using BenchmarkTools

# Define the function
function projection!(x::Vector{Float64}, lb::Vector{Float64}, ub::Vector{Float64})
	@inbounds @simd for i in eachindex(x, lb, ub)
		x[i] = min(max(x[i], lb[i]), ub[i])
	end
end

# Create test data
n = 10^6
x = randn(n)
lb = fill(-1.0, n)
ub = fill(1.0, n)

# Benchmark the function
@btime projection!($x, $lb, $ub)


# example
function inc_energy(γ::Vector{Float64})
	return γ' * γ
end

function inc_energy_grad(γ::Vector{Float64})
    return 2 * γ #ForwardDiff.gradient(inc_energy, γ)
end

function inc_energy_hess(γ::Vector{Float64})
    return 2*I #ForwardDiff.hessian(inc_energy, γ)
end

x0 = [1.,2.,3.,4.,5.].^2
inc_energy(x0)
inc_energy_grad(x0)
inc_energy_hess(x0)*x0

lb = -1*[10000.,-2.,-3.,-4.,-5.]/100
ub = [1.,2.,3.,4.,5.]*10

# model
@inline m(x) = inc_energy_grad(x)' * x + 0.5 * x' * (inc_energy_hess(x) * x);
@inline m_grad(x) = inc_energy_grad(x) + inc_energy_hess(x) * x;
@inline m_hess(x) = inc_energy_hess(x)

x = x0

g = m_grad(x0)
d = -m_grad(x0)  # Steepest descent direction

n = 5;
t_breaks = zeros(5);
n_breaks = 0
for i in 1:n
	if d[i] > 0 && ub[i] < Inf
		n_breaks += 1
		t_breaks[n_breaks] = (ub[i] - x[i]) / d[i]
	elseif d[i] < 0 && lb[i] > -Inf
		n_breaks += 1
		t_breaks[n_breaks] = (lb[i] - x[i]) / d[i]
	end
end

t_breaks
t_breaks_sorted = sort(t_breaks)
unique!(t_breaks_sorted)

p = similar(x)
t_prev = 0;
t_sum = 0
j = 5;
t = t_breaks_sorted[j]
Δt = t - t_prev

# g = m_grad(x)
# d = -m_grad(x)

@inbounds for i in 1:n
	p[i] = (t <= t_breaks[i]) ? d[i] : 0.0
	println("t = ")
end

d = p

q_prime = dot(g, d) + dot(x, m_hess(x0) * d)
q_double_prime = dot(d, m_hess(x0) * d)
Δt_star = - q_prime / q_double_prime

if q_prime > 0 || ( q_prime == 0 && q_double_prime > 0 )
	println("minimizer lies at t_i i.e. ", t)
else
	println("minimizer is not in t_i = ", t)
end

if (q_prime <= 0 && q_double_prime <= 0) || ( q_prime < 0 && q_double_prime > 0 && Δt_star >= Δt)
	println("minimizer lies at or beyond t_i+1 ")
	x = x + Δt * p
	t_sum += Δt
else
	println("minimizer is not in t_i+1 ?? ")
end


if q_prime < 0 && q_double_prime > 0 && Δt_star < Δt
	println("minimizer is at t_i - q_prime/q_double_prime ")
else
	println("3rd condition not true go further ")
end

t_sum

t

(-q_prime / q_double_prime)


t_prev = t

x

x = x - q_prime / q_double_prime * d


prepend!(t_breaks_sorted, [0])

t_prev 		= 0.0
x_prev 		= similar(x)
p 			= similar(x)
x_cauchy 	= similar(x)
x_segment 	= similar(x)

copyto!(x_prev, x)
q_prev = inc_energy(x_prev)
found_cauchy = false


for j in 1:1
	t = t_breaks_sorted[j]
	Δt = t - t_prev
	@inbounds for i in 1:n
		p[i] = (t <= t_breaks[i]) ? -d[i] : 0.0
		println("t = ")
	end
	t_prev = t
end

p

t = t_breaks[1]
Δt = t - t_prev

x_segment .= x_prev .+ Δt .* p

projection!(x_segment, lb, ub)

for j in 1:n_breaks
	t = t_breaks[j]
	Δt = t - t_prev
	@inbounds for i in 1:n
		p[i] = (x_prev[i] == x[i]) ? d[i] : 0.0
	end
	println("ppp = ", p , "  j = $j")
	x_segment .= x_prev .+ Δt .* p
	projection!(x_segment, lb, ub)

	# Quadratic function on this segment
	q_segment = m(x_segment)
	q_diff = q_segment - q_prev

	if q_diff >= 0  # Local minimum found
		found_cauchy = true
		break
	end

	t_prev = t
	copyto!(x_prev, x_segment)
	q_prev = q_segment
end

x_prev
