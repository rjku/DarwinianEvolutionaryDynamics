# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.1
#   kernelspec:
#     display_name: Julia (4 threads) 1.4.0
#     language: julia
#     name: julia-(4-threads)-1.4
# ---

# # Testing

using LinearAlgebra, Statistics, PyPlot, Revise, Distances
import mUtils

# ### Fictitious Fitness Landscape

logistic(x,x0,k,l) = l/( 1 + exp(k*(x-x0)) )

β = 0.005
F = 1
l = 10
k = 1
L = 40
s1 = [ exp( -β*((i-L)^2 + (j-L)^2) ) + exp( -β*(i^2 + j^2) ) for i in 0:L, j in 0:L ];
s2 = [ logistic(i,l,k,F)*logistic(j,l,k,F) for i in 0:L, j in 0:L ];

imshow(s2)

# ### Random Matrices

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

# ## Plotting Box Plots

using PyPlot
data = rand(100);
subplots(1,1,figsize=(5,5))
boxplot([data,data], positions=[.5, 2]);

# ## Plotting Grids

using Revise, mUtils, PyPlot, mGraphs

# +
X(i,L) = (i-1)%L + 1
Y(i,L) = (i-1)÷L + 1

L = 3
G = FiniteSquareLattice(L)

neighbors(G, 1)

# +
Vo = map( e -> [X(e,G.L),Y(e,G.L)], G.V[1] )
Vt = map( e -> [X(e,G.L),Y(e,G.L)], G.V[2] )

fig = figure(figsize=(5,5))

for v in V
    arrow( (v[2] + 2(dx*(v[1] - v[2])) + dx*dict[v[1] - v[2]])..., (Δx*(v[1] - v[2]))...,
        head_length=2dx, width= 0.05*rand(), ec="tab:blue",fc="tab:blue")   
end

# +
L = 5
N = 10

aPG1 = [ rand() for i in 1:L^2 ]
aPG2 = [ rand(4) for i in 1:L^2 ]

dx = 0.1
ℓ = 1.0 - 3dx

d = Dict( 0 => [0,0], 1 => [1,0], 2 => [0,1], 3 => [-1,0], 4 => [0,-1] )
o = Dict( 0 => [0,0], 1 => [0,-1], 2 => [1,0], 3 => [0,1], 4 => [-1,0] )

V = [ [X(i,L),Y(i,L)] for i in 1:L^2 ]

fig = figure(figsize=(7,7))
for (i,v) in enumerate(V)
    for (j,p) in enumerate(aPG2[i])
        arrow( (v + dx*d[ j ] + dx/2*o[ j ])..., ℓ*d[ j ]...,
            head_length=dx, width= 0.05*p, ec="tab:blue",fc="tab:blue")
    end
    
    scatter([v[1]], [v[2]], s=[(30 * aPG1[i])^2], cmap="cividis", c=10*[aPG1[i]], alpha=0.5)
end

# N = 2
# x = [1,2]
# y = [1,2]
# colors = rand(N)
# area = (30 * rand(N)).^2  # 0 to 15 point radii

# scatter(x, y, s=area, cmap="cividis", c=colors, alpha=0.5)

xticks(collect(0:L+1))
yticks(collect(0:L+1))

# +
dict = Dict([1,0] => [0,-1], [0,1] => [1,0], [-1,0] => [0,1], [0,-1] => [-1,0])
dx = 0.05
Δx = 0.70

V = [ [[1,0],[0,0]], [[0,1],[0,0]], [[0,0],[1,0]], [[1,1],[1,0]], [[1,1],[1,0]], [[1,0],[1,1]] ]

fig = figure(figsize=(5,5))
ax = fig.add_subplot(111)
xticks(collect(-1:2))
yticks(collect(-1:2))

for v in V
    arrow( (v[2] + 2(dx*(v[1] - v[2])) + dx*dict[v[1] - v[2]])..., (Δx*(v[1] - v[2]))..., head_length=2dx, width= 0.05*rand(), ec="tab:blue",fc="tab:blue")   
end

N = 2
x = [0,1]
y = [0,1]
colors = rand(N)
area = (30 * rand(N)).^2  # 0 to 15 point radii

scatter(x, y, s=area, c=colors, alpha=0.5)
show()
# -

[ range(N+1:N+3) for N in 0:3:10 ]
