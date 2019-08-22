
using mEvoTypes, mEvoFunc, PyPlot

const NGEN, NPOP, DGTY = Int32(1000), Int32(1000), 2
const REPRATE, MUTRATE = 10., 0.5
const DELTAX, FITNESSOFFSET, DELTATOFFSET = 0.01, 1.0, 1.0

# TRIVIAL TEST
trivialEty = tTrivialEty()
trivialEnv = tTrivialEnv()
aTrivialGty = [ tVecGty{Float64}(zeros(Float64,DGTY)) for i in 1:3NPOP ]
trivialPop = tLivingPop{Float64}( Int32(NPOP),trivialEty,trivialEnv,aTrivialGty,REPRATE,MUTRATE,DELTAX,DELTATOFFSET )

trivialData = tEvoData(NGEN)
@time evolution!(trivialPop, trivialData)

plot(collect(1:NGEN),trivialData.aveFitness); gcf()

# return plot(collect(1:NGEN),μ)
# ISING TEST
# SYSTEMSIZE = Int32(4)
# NGEN, NPOP, DGTY = 1, Int32(4), 2 #2SYSTEMSIZE^2

# isingEty = tIsingSigTransEty(SYSTEMSIZE,0.7,0.5)
# isingEnv = tCompEnv{Float64}([ [-5.0,-1.0], [5.0,1.0] ])
# aIsingGty = [ tVecGty{Float64}(ones(Float64,2SYSTEMSIZE^2)) for i in 1:3NPOP ]
# isingPop = tLivingPop{Float64}( NPOP,isingEty,isingEnv,aIsingGty,REPRATE,MUTRATE )

# trivialData = tEvoData(NGEN)
# evolution!(trivialPop, trivialData)
#
# plot(collect(1:NGEN),trivialData.aveFitness); gcf()

# the (initial) interaction matrix
# interactionmatrix = ones(Float64,2*SYSTEMSIZE^2)
# interactionmatrix = tVecGty{Float64}(2SYSTEMSIZE^2,rand(Float64,2SYSTEMSIZE^2))
# @time println(fitness(isingST, interactionmatrix, isingSTenv))

# aves = Array{Float64}(undef, SYSTEMSIZE, SYSTEMSIZE)
#
# @time println(metropolis(isingST,n,interactionmatrix,-10.0))
# @time metropolis(isingST,n,interactionmatrix,-10.0,aves)
# Juno.@run metropolis(SYSTEMSIZE,interactionmatrix)
