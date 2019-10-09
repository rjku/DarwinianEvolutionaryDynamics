
using Revise, mEvoTypes, Statistics
import mEvoFunc

# *******************
# *   ISING TEST	*
# *******************

const NBATCHES, NGEN, NPOP, SYSTEMSIZE = 1, Int32(50), Int32(100), Int32(6)			# 100, 30, 30 -> 15'
const INVTEMPERATURE, EXTERNALFIELD = 0.8, 0.1
const REPRATE, MUTRATE, SELSTRENGTH = 10.0^5, 10.0^4, 1.
const DELTAG, DELTATOFFSET = 1/3, 0.1
const NITERATION, NSAMPLINGS, NTRIALS = Int32(5*10^3), Int32(50), Int32(2)

# myFancyG = Float64[ i%(20*SYSTEMSIZE) <= 3*SYSTEMSIZE ? 1.5 : -3. for i in 1:2SYSTEMSIZE^2 ]
# myFancyIsing = tVecGty{Vector{Float64}}(myFancyG)

isingDTMCprm = tDTMCprm( NITERATION,NSAMPLINGS,NTRIALS )
aIsingMGty, aIsingGty, Npop = mEvoFunc.read_aIsingSigTransGty(isingDTMCprm,"test_#1")

isingEnv = tCompEnv([ [-10.0^10,-1.0], [10.0^10,1.0] ],SELSTRENGTH)
isingEty = mEvoFunc.tDscdtEty{Float64}(REPRATE,MUTRATE,DELTATOFFSET,aIsingMGty[end].dG,DELTAG)

isingPop = tLivingPop( Int32[Npop,Npop,Npop],isingEty,isingEnv,aIsingMGty,aIsingGty )

aIsingData = tEvoData[]
# for i in 1:NBATCHES
# 	Fave = mean( [isingPop.aGty[i].pF[1] for i in 1:isingPop.pN[2]] )
# 	push!( aIsingData, tEvoData(NGEN, Fave + ( 1. - ( Fave % 1 ) ) % (1/3) + .2, 1/3) )
# 	@time mEvoFunc.evolution!(isingPop,aIsingData[end],ubermode=true)
# end

Fave = mean( [isingPop.aGty[i].pF[1] for i in 1:isingPop.pN[2]] )
push!(aIsingData,tEvoData(NGEN,Fave + ( 1. - ( Fave % 1 ) ) % (1/3) + .2, 1/3) )
push!(aIsingData[end].aLivingPop,deepcopy(isingPop))
