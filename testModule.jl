
using mEvolution, mIsingEvolution

SYSTEMSIZE = Int32(10)
NGEN, NPOP, DGTY = 1, Int32(5), 2SYSTEMSIZE^2

# the (initial) state vector
# n = rand([-1,1], 	L, L)
# n = ones(Int8, SYSTEMSIZE, SYSTEMSIZE)

isingEty = tIsingSigTransEty(SYSTEMSIZE,0.7,0.5)

isingEnv = tCompEnv{Float64}([ [-5.0,-1.0], [5.0,1.0] ])

aIsingGty = [ tTupGty{Float64}(rand(Float64,2SYSTEMSIZE^2)) for i in 1:3NPOP ]

isingPop = tLivingPop{Float64}(NPOP,isingEty,isingEnv,aIsingGty)

# the (initial) interaction matrix
# interactionmatrix = ones(Float64,2*SYSTEMSIZE^2)
# interactionmatrix = tVecGty{Float64}(2SYSTEMSIZE^2,rand(Float64,2SYSTEMSIZE^2))
# @time println(fitness(isingST, interactionmatrix, isingSTenv))

# aves = Array{Float64}(undef, SYSTEMSIZE, SYSTEMSIZE)
#
# @time println(metropolis(isingST,n,interactionmatrix,-10.0))
# @time metropolis(isingST,n,interactionmatrix,-10.0,aves)
# Juno.@run metropolis(SYSTEMSIZE,interactionmatrix)
