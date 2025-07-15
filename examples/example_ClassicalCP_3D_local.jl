
using ForwardDiff
using LinearAlgebra
using StaticArrays
include("../src/algorithms/cauchy_point.jl")
include("../src/algorithms/trust_region_LA.jl")
include("../src/algorithms/trust_region_frac_to_boundary.jl")

function count_positives_negatives(vec::AbstractVector{T}) where T <: Real
	positives = sum(x -> x > 1e-14, vec)  # Count positives
	negatives = sum(x -> x < 1e-14, vec)  # Count negatives

	println("Number of positive elements: $positives")
	println("Number of negative elements: $negatives")
end

# elastic matrix
struct mat_data{T, M <: Matrix{T}}
	c11::T
	c12::T
	c44::T
	Dᵉ::M
end

function elastic_data(c11 = 170.0, c12 = 124.0, c44 = 75.0)
    Cv6o = [
        c11 c12 c12 0.0 0.0 0.0
        c12 c11 c12 0.0 0.0 0.0
        c12 c12 c11 0.0 0.0 0.0
        0.0 0.0 0.0 2.0*c44 0.0 0.0
        0.0 0.0 0.0 0.0 2.0*c44 0.0
        0.0 0.0 0.0 0.0 0.0 2.0*c44
        ]
    return mat_data(c11, c12, c44, Cv6o)
end;

struct history_data{T, S <: Matrix}
    # Store "converged" values
    ϵᵖ	::S # plastic strain
    σ 	::S # stress
    λ 	::T # deforamation state
    γ_tot::T # total slip
end;

function history_data()
	return history_data(
				zeros(3, 3),
				zeros(3, 3),
				0.0,
				0.0)
end;

# vector to matrix transformation
function V6toM(vec)
	sq2 = sqrt(2.0)
	m = [
		vec[1]			vec[6] / sq2	vec[5] / sq2;
		vec[6] / sq2 	vec[2]			vec[4] / sq2;
		vec[5] / sq2 	vec[4] / sq2	vec[3]
    	]
    return m
end;

# matrix to vector transformation
function MtoV6(mat) #symmetric ??
	sq2 = sqrt(2.0)
	return [mat[1, 1], mat[2, 2], mat[3, 3], sq2 * mat[2, 3], sq2 * mat[1, 3], sq2 * mat[1, 2]]
end;

# mα ⊗ nβ
function Pα()

	# Define mFCC
	mFCC = [
	[1, 1, 1], [1, 1, 1], [1, 1, 1],
	[1, 1, -1], [1, 1, -1], [1, 1, -1],
	[1, -1, 1], [1, -1, 1], [1, -1, 1],
	[-1, 1, 1], [-1, 1, 1], [-1, 1, 1],
	[1, 1, 1], [1, 1, 1], [1, 1, 1],
	[1, 1, -1], [1, 1, -1], [1, 1, -1],
	[1, -1, 1], [1, -1, 1], [1, -1, 1],
	[-1, 1, 1], [-1, 1, 1], [-1, 1, 1]
	] / sqrt(3.0)

	# Define sFCC
	sFCC = [
	[1, -1, 0], [1, 0, -1], [0, 1, -1],
	[1, -1, 0], [1, 0, 1], [0, 1, 1],
	[1, 1, 0], [1, 0, -1], [0, 1, 1],
	[1, 1, 0], [1, 0, 1], [0, 1, -1],
	-[1, -1, 0], -[1, 0, -1], -[0, 1, -1],
	-[1, -1, 0], -[1, 0, 1], -[0, 1, 1],
	-[1, 1, 0], -[1, 0, -1], -[0, 1, 1],
	-[1, 1, 0], -[1, 0, 1], -[0, 1, -1]
	] / sqrt(2.0)

	return [0.5 * ((mFCC[i] * sFCC[i]')' + (mFCC[i] * sFCC[i]')) for i in 1:24]
end;

# interaction matrix
function hαβ()
    dim = 24;
    MxH = 1.4 * ones(Float64, dim, dim);
    MxH[1:dim+1:end] .= 1.0;
    return MxH
end;

function ϵt(λ)
	Id3 = [1 0 0; 0 1 0; 0 0 1]
	F = [1. λ 0.; 0. 1. 0.; 0. 0. 1.]
	ϵt = 0.5 * (F' * F - Id3)
	return ϵt
end;

material_data = elastic_data();
history_state = history_data();
Dᵉ = material_data.Dᵉ;
P = Pα();
hab = hαβ();
mCm = [MtoV6(P[i])' * Dᵉ * MtoV6(P[j]) for i in 1:24, j in 1:24];

function yield_f(γ::AbstractVector{Float64}, state::history_data)

	λn = state.λ
	λ = λn + 0.001
	ϵpn = state.ϵᵖ

	# Replace temporary allocation with StaticArrays
	ϵp = ϵpn + sum(P[i] * γ[i] for i in 1:24)

	ϵ = ϵt(λ)
	ϵe = ϵ - ϵp
	σ = V6toM(Dᵉ * MtoV6(ϵe))

	τ = [sum(σ .* P[i]) for i in 1:24]
	τc = hab * γ
	f = τ - τc .- 0.001
	return f
end;

yield_TR(x::AbstractVector{Float64}) = yield_f(x, history_state)

yield_TR(zeros(24))
ftrial = yield_f(zeros(24), history_state)

active = [ftrial[i] > 0. ? 1 : 0 for i in 1:24]

count_positives_negatives(ftrial)
println(ftrial)

function inc_energy(γ::Vector{Float64})
	return - ftrial' * γ + 0.5 * γ' * ( ( hab + mCm ) * γ)
end

function inc_energy_grad(γ::Vector{Float64})
	yield = - ftrial + ( ( hab + mCm ) * γ) #ForwardDiff.gradient(inc_energy, γ)
	return yield
end

function inc_energy_hess(γ::Vector{Float64})
	return (hab + mCm) #ForwardDiff.hessian(inc_energy, γ)
end

ftrial .+ inc_energy_grad(dγ)

hess = inc_energy_hess(dγ)

println(eigen(hess).values)

delta = 1.;
lb = zeros(24)
ub = [Inf for i in 1:24]
dγ = zeros(24)

dγ = [dγ[i] > 0. ? dγ[i] : 0. for i in 1:24]

active
c = -inc_energy_grad(dγ)
println(c)

penalty = 100.;

# inc_energy(dγ)
# inc_energy_grad(dγ)
# inc_energy_hess(dγ)

# dγ = cauchy_point(delta, inc_energy_grad(dγ), inc_energy_hess(dγ) )


active_f = [i for i in 1:24] .* active

active_final = [active_f[i] for i in 1:24 if active[i] > 0]

count_positives_negatives(dγ)
sol1 = trust_region_LA(inc_energy, inc_energy_grad, inc_energy_hess, yield_TR, delta, dγ, penalty, lb; verbose = 0)

start = zeros(24)
start[active_final] .+= 0.0001 .+randn(8)/1000000
start
sol2 = trust_region_ftb(inc_energy, inc_energy_grad, inc_energy_hess,  delta, start)
println(ftrial)
println(sol2)
println(sol1)

println(yield_f(sol1, history_state))
println(yield_f(sol2, history_state))

inc_energy(sol1)
inc_energy(sol2)

println(yield_f(dγ, history_state))
yield_f(dγ, history_state)' * dγ
yield_f(dγ, history_state) + inc_energy_grad(dγ, active)

println(yield_f(dγ, history_state))

-inc_energy_grad(dγ, active)' * dγ

inc_energy(sol)

dγ
count_positives_negatives(dγ)

ftrial' * dγ

- inc_energy_grad(dγ) .* dγ


ftrial = yield_f(dγ, history_state)
println(ftrial)
