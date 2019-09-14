
using mEvoTypes, mEvoFunc

# *******************
# *   ISING TEST	*
# *******************

const NGEN, NPOP, SYSTEMSIZE = Int32(3), Int32(3), Int32(6)			# 100, 30, 30 -> 15'
const INVTEMPERATURE, EXTERNALFIELD = 1.0, 0.2
const REPRATE, MUTRATE, SELSTRENGTH = 10^5., 50., 1.
const DELTAX, DELTATOFFSET = .5, 0.1
const NITERATION, NSAMPLINGS = Int32(10^4)*SYSTEMSIZE^2, Int32(100)

# myFancyX = Float64[ i%(20*SYSTEMSIZE) <= 3*SYSTEMSIZE ? 1.5 : -3. for i in 1:2SYSTEMSIZE^2 ]
# myFancyIsing = tVecGty{Vector{Float64}}(myFancyX)

# isingEty = tEty(REPRATE,MUTRATE,DELTATOFFSET,Int32(2SYSTEMSIZE^2),DELTAX)
isingDTMCprm = tDTMCprm( NITERATION,NSAMPLINGS )
isingEnv = tCompEnv([ [-10.0^10,-1.0], [10.0^10,1.0] ],SELSTRENGTH)
aIsingMGty = [ tIsingSigTransMGty(SYSTEMSIZE,INVTEMPERATURE,EXTERNALFIELD,isingDTMCprm) ]
aIsingGty = [ tVecGty( [aIsingMGty[1]], rand(-4:.5:5,2SYSTEMSIZE^2) ) for i in 1:3NPOP ]

isingPop = initLivingPop( NPOP,tEty{Float64}( REPRATE,MUTRATE,DELTATOFFSET,Int32(2SYSTEMSIZE^2),DELTAX ),isingEnv,aIsingMGty,aIsingGty )

# write_aGty(isingPopClone)

# aIsingMGtyClone, aIsingGtyClone = read_aIsingSigTransGty( "",INVTEMPERATURE,EXTERNALFIELD )
# isingPopClone = tLivingPop(
# 	Int32[length(aIsingGtyClone),length(aIsingGtyClone),length(aIsingGtyClone)],
# 	tEty{Float64}( REPRATE,MUTRATE,DELTATOFFSET,Int32(2SYSTEMSIZE^2),DELTAX ),isingEnv,aIsingMGtyClone,aIsingGtyClone )

isingData = tEvoData(NGEN)
@time evolution!(isingPop,isingData,ubermode=true)

aAves = [ zeros(Float64,SYSTEMSIZE,SYSTEMSIZE) for i in isingEnv.idealInputOutput ]
@time showPhenotype!(isingEnv,isingPop.aGty[1],aAves)
# @time showPhenotype!(isingEnv,isingPopClone.aGty[2],aAves)
# @time showPhenotype!(isingEnv,myFancyIsing,aAves)

JijMat = Array{Float64}(undef,2isingPop.aGty[1].pMetaGty[1].L,2isingPop.aGty[1].pMetaGty[1].L)
@time showGenotype!(isingPop.aGty[1],JijMat)
# @time showGenotype!(isingPopClone.aGty[1],JijMat)
# @time showGenotype!(myFancyIsing,JijMat)

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
