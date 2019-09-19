
using Revise, mEvoTypes
import mEvoFunc

# *******************
# *   ISING TEST	*
# *******************

const NGEN, NPOP, SYSTEMSIZE = Int32(10), Int32(50), Int32(6)			# 100, 30, 30 -> 15'
const INVTEMPERATURE, EXTERNALFIELD = 0.8, 0.1
const REPRATE, MUTRATE, SELSTRENGTH = 10.0^5, 10.0^4, 1.
const DELTAX, DELTATOFFSET = .5, 0.1
const NITERATION, NSAMPLINGS, NTRIALS = Int32(10^3), Int32(50), Int32(2)

# myFancyX = Float64[ i%(20*SYSTEMSIZE) <= 3*SYSTEMSIZE ? 1.5 : -3. for i in 1:2SYSTEMSIZE^2 ]
# myFancyIsing = tVecGty{Vector{Float64}}(myFancyX)

isingDTMCprm = tDTMCprm( NITERATION,NSAMPLINGS,NTRIALS )
aIsingMGty = [ tIsingSigTransMGty(SYSTEMSIZE,INVTEMPERATURE,EXTERNALFIELD,isingDTMCprm) ]
aIsingGty = [ tVecGty( [aIsingMGty[1]], rand(-4:.5:5,2SYSTEMSIZE^2) ) for i in 1:3NPOP ]

isingEnv = tCompEnv([ [-10.0^10,-1.0], [10.0^10,1.0] ],SELSTRENGTH)
isingEty = mEvoFunc.tEty{Float64}(REPRATE,MUTRATE,DELTATOFFSET,aIsingMGty[end].dX,DELTAX)

isingPop = mEvoFunc.initLivingPop( NPOP,isingEty,isingEnv,aIsingMGty,aIsingGty )
isingPop = tLivingPop( Int32[NPOP,NPOP,length(aIsingGty)],isingEty,isingEnv,aIsingMGty,aIsingGty )

# aIsingMGtyClone, aIsingGtyClone = read_aIsingSigTransGty( "",INVTEMPERATURE,EXTERNALFIELD )
# isingPopClone = tLivingPop(
# 	Int32[length(aIsingGtyClone),length(aIsingGtyClone),length(aIsingGtyClone)],
# 	tEty{Float64}( REPRATE,MUTRATE,DELTATOFFSET,Int32(2SYSTEMSIZE^2),DELTAX ),isingEnv,aIsingMGtyClone,aIsingGtyClone )

# for iterations in 1:10

aIsingData = tEvoData[]

push!(aIsingData,tEvoData(NGEN))
@time mEvoFunc.evolution!(isingPop,aIsingData[end],ubermode=true)

mEvoFunc.write_aGty(isingPop, "test" * "#$(length(aIsingData))")

aAves = [ zeros(Float64,isingPop.aGty[1].pMetaGty[1].L,isingPop.aGty[1].pMetaGty[1].L) for i in isingEnv.idealInputOutput ]
@time mEvoFunc.showPhenotype!(isingEnv,isingPop.aGty[1],aAves)
# @time showPhenotype!(isingEnv,isingPopClone.aGty[2],aAves)
# @time showPhenotype!(isingEnv,myFancyIsing,aAves)

JijMat = Array{Float64}(undef,2isingPop.aGty[1].pMetaGty[1].L,2isingPop.aGty[1].pMetaGty[1].L)
@time mEvoFunc.showGenotype!(isingPop.aGty[1],JijMat)
# @time showGenotype!(isingPopClone.aGty[1],JijMat)
# @time showGenotype!(myFancyIsing,JijMat)

# check upgrade
# mEvoFunc.upgradeGtyX!(isingPop.aGty[1])
# push!( isingPop.aMetaGty,tIsingSigTransMGty(
# 	isingPop.aGty[1].pMetaGty[1].L+Int32(2), isingPop.aGty[1].pMetaGty[1].β, isingPop.aGty[1].pMetaGty[1].he, isingPop.aGty[1].pMetaGty[1].prms
# 	))
# isingPop.aGty[1].pMetaGty[1] = isingPop.aMetaGty[end]


# *******************
# TRIVIAL TEST
# *******************

# const NGEN, NPOP, DGTY = Int32(2000), Int32(50), 2(30^2)
# const REPRATE, MUTRATE = 10000., .1
# const DELTAX, DELTATOFFSET = 0.01, 1.0
#
# trivialEty = tTrivialEty()
# trivialEnv = tTrivialEnv()
# aTrivialGty = [ tVecGty{Array{Float64,1}}(zeros(Float64,DGTY)) for i in 1:3NPOP ]
# trivialPop = tLivingPop{Float64,tTrivialEty,tTrivialEnv,Array{tVecGty,1}}(
# 	Int32(NPOP),trivialEty,trivialEnv,aTrivialGty,REPRATE,MUTRATE,DELTAX,DELTATOFFSET )
#
# trivialData = tEvoData(NGEN)
# @time evolution!(trivialPop, trivialData, ubermode=true)
#
# plot(collect(1:NGEN),trivialData.aveFitness); gcf()
#
# clf()
