
using mEvoTypes, mEvoFunc

# *******************
# *   ISING TEST	*
# *******************

const NGEN, NPOP, SYSTEMSIZE = Int32(10), Int32(2), Int32(6)			# 100, 30, 30 -> 15'
const INVTEMPERATURE, EXTERNALFIELD = 1.0, 0.2
const REPRATE, MUTRATE, SELSTRENGTH = 10^5., 2., 1.
const DELTAX, DELTATOFFSET = .5, 0.1
const NITERATION, NSAMPLINGS = Int32(10^4)*SYSTEMSIZE^2, Int32(100)

# myFancyX = Float64[ i%(20*SYSTEMSIZE) <= 3*SYSTEMSIZE ? 1.5 : -3. for i in 1:2SYSTEMSIZE^2 ]
# myFancyIsing = tVecGty{Vector{Float64}}(myFancyX)

# isingEty = tEty(REPRATE,MUTRATE,DELTATOFFSET,Int32(2SYSTEMSIZE^2),DELTAX)
isingEnv = tCompEnv([ [-10.0^10,-1.0], [10.0^10,1.0] ],SELSTRENGTH)
aIsingMGty = [ tIsingSigTransMGty(SYSTEMSIZE,INVTEMPERATURE,EXTERNALFIELD) ]
aIsingGty = [ tVecGty( [aIsingMGty[1]], rand(-4:.5:5,2SYSTEMSIZE^2) ) for i in 1:3NPOP ]
isingDTMCprm = tDTMCprm(NITERATION,NSAMPLINGS)

isingPop = initLivingPop( NPOP,tEty{Float64}( REPRATE,MUTRATE,DELTATOFFSET,Int32(2SYSTEMSIZE^2),DELTAX ),isingEnv,aIsingMGty,aIsingGty,isingDTMCprm )

write_aGty(isingPop)

aIsingMGtyClone, aIsingGtyClone = read_aIsingSigTransGty( "population_2019-09-12T17:25:39.612.dat",INVTEMPERATURE,EXTERNALFIELD )
isingPopClone = initLivingPop(
	NPOP,tEty{Float64}( REPRATE,MUTRATE,DELTATOFFSET,Int32(2SYSTEMSIZE^2),DELTAX ),isingEnv,aIsingMGtyClone,aIsingGtyClone,isingDTMCprm )

isingData = tEvoData(NGEN)
@time evolution!(isingPop,isingData,isingDTMCprm,ubermode=true)

aAves = [ zeros(Float64,SYSTEMSIZE,SYSTEMSIZE) for i in isingEnv.idealInputOutput ]
@time showPhenotype!(isingEnv,isingPop.aGty[1],aAves,isingDTMCprm)
# @time showPhenotype!(isingEnv,isingPopClone.aGty[2],aAves,isingDTMCprm)
# @time showPhenotype!(isingEnv,myFancyIsing,aAves,isingDTMCprm)

JijMat = Array{Float64}(undef,2SYSTEMSIZE,2SYSTEMSIZE)
@time showGenotype!(isingPop.aGty[1],JijMat)
# @time showGenotype!(isingPopClone.aGty[2],JijMat)
# @time showGenotype!(myFancyIsing,JijMat)

# figure(1); clf()
# matshow(aAves[1],cmap="Greys_r",vmin=-1,vmax=1); gcf()
# matshow(aAves[2],cmap="Greys_r",vmin=-1,vmax=1); gcf()

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
