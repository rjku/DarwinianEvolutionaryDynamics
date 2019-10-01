
using Revise, mEvoTypes, Statistics
import mEvoFunc

# *******************
# *   ISING TEST	*
# *******************

const NBATCHES, NGEN, NPOP, SYSTEMSIZE = 1, Int32(20), Int32(100), Int32(6)			# 100, 30, 30 -> 15'
const INVTEMPERATURE, EXTERNALFIELD = 0.8, 0.1
const REPRATE, MUTRATE, SELSTRENGTH = 10.0^7, 10.0^3, 1.
const DELTAG, DELTATOFFSET = 1/3, 0.01
const NITERATION, NSAMPLINGS, NTRIALS = Int32(10^3), Int32(10^2), Int32(2)

# myFancyG = Float64[ i%(20*SYSTEMSIZE) <= 3*SYSTEMSIZE ? 1.5 : -3. for i in 1:2SYSTEMSIZE^2 ]
# myFancyIsing = tVecGty{Vector{Float64}}(myFancyG)

isingDTMCprm = tDTMCprm( NITERATION,NSAMPLINGS,NTRIALS )
aIsingMGty = [ tIsingSigTransMGty(SYSTEMSIZE,INVTEMPERATURE,EXTERNALFIELD,isingDTMCprm) ]
# aIsingGty = [ tVecGty( [aIsingMGty[1]], rand(-2:DELTAG:2,2SYSTEMSIZE^2) ) for i in 1:3NPOP ]
aIsingGty = [ tAlphaGty( [aIsingMGty[1]], rand(-2:DELTAG:2,2SYSTEMSIZE^2), collect(-2:DELTAG:2) ) for i in 1:3NPOP ]

isingEnv = tCompEnv([ [-10.0^10,-1.0], [10.0^10,1.0] ],SELSTRENGTH)
isingEty = mEvoFunc.tEty(REPRATE,MUTRATE,DELTATOFFSET,aIsingGty[1])

isingPop = mEvoFunc.initLivingPop( NPOP,isingEty,isingEnv,aIsingMGty,aIsingGty )
# isingPop = tLivingPop( Int32[NPOP,NPOP,length(aIsingGty)],isingEty,isingEnv,aIsingMGty,aIsingGty )

aIsingData = tEvoData[]
for i in 1:NBATCHES
	Fave = mean( [isingPop.aGty[i].pF[1] for i in 1:isingPop.pN[2]] )
	push!( aIsingData, tEvoData(NGEN, Fave + ( 1. - ( Fave % 1 ) ) % (1/3) + .2, 1/3) )
	@time mEvoFunc.evolution!(isingPop,aIsingData[end],ubermode=true)
end

push!( aIsingData[end].aLivingPop,deepcopy(isingPop) )

# mEvoFunc.write_tLivingPop(isingPop, "test_#2") # * "#$(length(aIsingData))")

# check upgrade
# mEvoFunc.upgradeGtyG!(isingPop.aGty[1],1/3)
# push!( isingPop.aMetaGty,tIsingSigTransMGty(
# 	isingPop.aGty[1].pMetaGty[1].L+Int32(2), isingPop.aGty[1].pMetaGty[1].β, isingPop.aGty[1].pMetaGty[1].he, isingPop.aGty[1].pMetaGty[1].prms
# 	))
# isingPop.aGty[1].pMetaGty[1] = isingPop.aMetaGty[end]


# *******************
# TRIVIAL TEST
# *******************

# const NGEN, NPOP, DGTY = Int32(2000), Int32(50), 2(30^2)
# const REPRATE, MUTRATE = 10000., .1
# const DELTAG, DELTATOFFSET = 0.01, 1.0
#
# trivialEty = tTrivialEty()
# trivialEnv = tTrivialEnv()
# aTrivialGty = [ tVecGty{Array{Float64,1}}(zeros(Float64,DGTY)) for i in 1:3NPOP ]
# trivialPop = tLivingPop{Float64,tTrivialEty,tTrivialEnv,Array{tVecGty,1}}(
# 	Int32(NPOP),trivialEty,trivialEnv,aTrivialGty,REPRATE,MUTRATE,DELTAG,DELTATOFFSET )
#
# trivialData = tEvoData(NGEN)
# @time evolution!(trivialPop, trivialData, ubermode=true)
#
# plot(collect(1:NGEN),trivialData.aveFitness); gcf()
#
# clf()
