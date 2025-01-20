
using ForwardDiff
using LinearAlgebra
using StaticArrays
include("../src/algorithms/cauchy_point.jl")
include("../src/algorithms/trust_region_LA.jl")
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
	F = [1. 0. 0.; 0. 1. 0.; 0. 0. -λ*1.]
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
	λ = λn + 0.0001
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

ftrial = yield_f(zeros(24), history_state)

function inc_energy(γ::Vector{Float64})
	return - ftrial' * γ + 0.5 * γ' * ( ( hab + mCm ) * γ)
end

function inc_energy_grad(γ::Vector{Float64})
    return - ftrial + ( ( hab + mCm ) * γ) #ForwardDiff.gradient(inc_energy, γ)
end

function inc_energy_hess(γ::Vector{Float64})
    return hab + mCm  #ForwardDiff.hessian(inc_energy, γ)
end


dγ = zeros(24)
γsol = zeros(24)
count_positives_negatives(ftrial)

delta = 1.1;
penalty = 10000.;
lower_bounds = zeros(24)

# inc_energy(dγ)
# inc_energy_grad(dγ)
# inc_energy_hess(dγ)

dγ = cauchy_point(delta, inc_energy_grad(dγ), inc_energy_hess(dγ) )
160*160
count_positives_negatives(dγ)
@time for i in 1:10000
	γsol .= trust_region_LA(inc_energy, inc_energy_grad, inc_energy_hess, delta, dγ, penalty, lower_bounds; verbose = 0)
end
count_positives_negatives(γsol)

ftrial = yield_f(γsol, history_state)
println(ftrial)
println(inc_energy_grad(γsol))
println(γsol)
count_positives_negatives(γsol)

norm(ftrial .* γsol)
