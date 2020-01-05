# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.0
#   kernelspec:
#     display_name: Julia (2 threads) 1.3.0
#     language: julia
#     name: julia-(2-threads)-1.3
# ---

# # Testing

# ### Random Matrices

using LinearAlgebra, Statistics, PyPlot, Revise
import mUtils

# +
N = 1000; M = 100000; c = N/M

println(c)

testMat = randn(N,M);
W = testMat * testMat' ./ M

λval = [ λ for λ in mUtils.MPpdf(0,c)[2]:0.01:mUtils.MPpdf(0,c)[3] ]
MPval = [ mUtils.MPpdf(λ,c)[1] for λ in mUtils.MPpdf(0,c)[2]:0.01:mUtils.MPpdf(0,c)[3] ];

hist(eigvals(Symmetric(W)),40,density=true,log=false,rwidth=.9); plot(λval,MPval)

# +
N = 336; M = 10000; c = N/M

testMat = rand((-.8,.8),N,M);
W = testMat * testMat' ./ M

λval = [ λ for λ in mUtils.MPpdf(0,c)[2]:0.01:mUtils.MPpdf(0,c)[3] ]
MPval = [ mUtils.MPpdf(λ,c)[1] for λ in mUtils.MPpdf(0,c)[2]:0.01:mUtils.MPpdf(0,c)[3] ];

hist(eigvals(Symmetric(W)),40,density=true,log=false,rwidth=.9); plot(λval,MPval)

# +
N = 360; M = 1000; c = N/M

testMat = -1.5 .+ rand(N,M) .* 3.0;
W = testMat * testMat' ./ M

λval = [ λ for λ in mUtils.MPpdf(0,c)[2]:0.01:mUtils.MPpdf(0,c)[3] ]
MPval = [ mUtils.MPpdf(λ,c)[1] for λ in mUtils.MPpdf(0,c)[2]:0.01:mUtils.MPpdf(0,c)[3] ];

hist(eigvals(Symmetric(W)),40,density=true,log=false,rwidth=.9); plot(λval,MPval)
# -

N = 1000; M = 10000; c = 3
testMat = rand((0.0,1.0),N,N*c);
W = testMat * testMat'
hist(eigvals(Symmetric(W)),40,log=false,rwidth=.9);

N = 1000; M = 10000; c = 3
testMat = rand((-0.7,0.0,1.0),N,N*c);
W = testMat * testMat'
hist(eigvals(Symmetric(W)),40,log=false,rwidth=.9);

# ### Complex Stuff

G = [ exp(im*rand((0,pi/2,pi,3pi/2))) for i in 1:10 ]
mUtils.chopArray!(G)
mean(G)

G[1]*conj(G[2])

B = [ 1.0, im, -1.0, -im ]
real([ B[i] * conj(B)[j] for i in 1:4, j in 1:4 ])

# ## Plotting Grids

using Revise, mUtils, PyPlot

# +
dict = Dict([1,0] => [0,-1], [0,1] => [1,0], [-1,0] => [0,1], [0,-1] => [-1,0])
dx = .05
Δx = .70

V = [ [[1,0],[0,0]], [[0,1],[0,0]], [[0,0],[1,0]], [[1,1],[1,0]], [[1,1],[1,0]], [[1,0],[1,1]] ]

subplots(1,1,figsize=(5,5))
xticks(collect(-1:2))
yticks(collect(-1:2))

for v in V
    arrow( (v[2] + 2(dx*(v[1] - v[2])) + dx*dict[v[1] - v[2]])..., (Δx*(v[1] - v[2]))..., head_length=2dx, width= .05*rand(), ec="tab:blue",fc="tab:blue")
end
# -

[ range(N+1:N+3) for N in 0:3:10 ]


