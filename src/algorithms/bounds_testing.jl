include("bounds.jl")

using .SimpleBounds
using LinearAlgebra

f(x::Vector{Float64}) = x'* x

grad(x::Vector{Float64}) = 2x

hess(x::Vector{Float64}) = I * 2


dst = [0.,0.,0., 0.]
x0 = [60., 70., 65., 25.]
xnew = zeros(4)

lb = [35., 60., -10., 10.]
ub = [800., 800., 800., 800.]

g = grad(x0)

xnew = x0 - hess(x0) \ g

project_variables!(dst, x0, lb, ub)

project_direction!(dst, x0, lb, ub, -, g)


project_gradient!(dst, x0, lb, ub, g)

project_gradient!(dst, xnew, lb, ub, g)
