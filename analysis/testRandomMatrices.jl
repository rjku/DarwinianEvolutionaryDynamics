# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.2.4
#   kernelspec:
#     display_name: Julia 1.2.0
#     language: julia
#     name: julia-1.2
# ---

# # Test Games: Random Matrices

using LinearAlgebra, Statistics, PyPlot

N = 1000; M = 10000; c = 3
testMat = randn(N,N*c);
W = testMat * testMat'
hist(eigvals(Symmetric(W)),40,log=false,rwidth=.9);

N = 1000; M = 10000; c = 3
testMat = rand((-1.0,1.0),N,N*c);
W = testMat * testMat'
hist(eigvals(Symmetric(W)),40,log=false,rwidth=.9);

N = 1000; M = 10000; c = 3
testMat = rand((0.0,1.0),N,N*c);
W = testMat * testMat'
hist(eigvals(Symmetric(W)),40,log=false,rwidth=.9);

N = 1000; M = 10000; c = 3
testMat = rand((-0.7,0.0,1.0),N,N*c);
W = testMat * testMat'
hist(eigvals(Symmetric(W)),40,log=false,rwidth=.9);

N = 1000; M = 10000; c = 3
testMat = rand((-1.0,0.0,1.0),N,N*c);
W = testMat * testMat'
hist(eigvals(Symmetric(W)),40,log=false,rwidth=.9);

W
