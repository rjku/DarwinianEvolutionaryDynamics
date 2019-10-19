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
#     display_name: Julia (4 threads) 1.2.0
#     language: julia
#     name: julia-(4-threads)-1.2
# ---

using Revise, BenchmarkTools, Statistics, LinearAlgebra, Statistics, PyPlot, mEvoTypes
import mEvoFunc

# %matplotlib notebook

const NBATCHES, NGEN, NPOP, SYSTEMSIZE, MAXSIZE = 8, Int32(25), Int32(25), Int32(6), Int32(12)
const INVTEMPERATURE, EXTERNALFIELD = 0.8, 0.1
const REPRATE, MUTRATE, SELSTRENGTH = 10.0^7, 10.0^3, 1.0
const REPFACTOR, PNTMUTFACTOR, SELSTRENGTH = 10.0, 0.05, 1.0
const DELTAG, DELTATOFFSET = 1/3, 0.01
const GMIN, GMAX = -2.0, 2.0
const NSAMPLINGS, NMCSPS, NTRIALS = Int32(20), Int32(10^5), Int32(1);

isingEnv = tCompEnv([ [-1.0,-1.0], [1,1.0] ],SELSTRENGTH)
isingEty = tPntMutEty(REPFACTOR, PNTMUTFACTOR, [ 2L^2 for L in SYSTEMSIZE:2:MAXSIZE ], 12);
isingDTMCprm = tDTMCprm( NSAMPLINGS, NMCSPS, NTRIALS );

aIsingMGty = [ tIsingSigTransMGty(SYSTEMSIZE,INVTEMPERATURE,EXTERNALFIELD,isingDTMCprm) ];
aIsingGty = [ tCntGty( [aIsingMGty[1]], GMIN .+ rand(2SYSTEMSIZE^2) .* (GMAX - GMIN), [GMIN,GMAX] ) for i in 1:3NPOP ];
aIsingData = tEvoData[];

@time isingPop = mEvoFunc.initLivingPop( NPOP,isingEty,isingEnv,aIsingMGty,aIsingGty );

for i in 1:NBATCHES
	Fave = mean( [isingPop.aGty[i].pF[1] for i in 1:isingPop.pN[2]] )
	push!( aIsingData, tEvoData(NGEN, Fave + ( 1. - ( Fave % 1 ) ) % (1/3) + .2, 1/3) )
    @time mEvoFunc.evolutionGKPup!(isingPop,aIsingData[end],ubermode=true)
end
push!( aIsingData[end].aLivingPop,deepcopy(isingPop) );

lastBatch = 2
subplots(3,1,figsize=(5,10))
subplot(311); title("Average Fitness"); plot( vcat([dataBtc.aveFitness for dataBtc in aIsingData[1:lastBatch]]...) );
subplot(312); title("Growth Factor"); plot( vcat([dataBtc.growthFactor for dataBtc in aIsingData[1:lastBatch]]...) );
subplot(313); title("Mutation Factor"); plot( vcat([dataBtc.mutationFactor for dataBtc in aIsingData[1:lastBatch]]...) );
tight_layout();


# +
iBatch, iPop, iIsing = 2, 2, 1

const NSMPLSTAT, NMCSPSSTAT, NTRLSTAT = Int32(2+10^1), Int32(10^5), Int32(2)
isingDTMCprmStat = tDTMCprm( NSMPLSTAT, NMCSPSSTAT, NTRLSTAT )

GMat = Array{Float64}(undef, 2aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L, 2aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L)
mEvoFunc.showG!(aIsingData[iBatch].aLivingPop[iPop].aGty[1], GMat)

aSpinAve = [ zeros(Float64, aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L, aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L) for i in isingEnv.IOidl ]
aSpinCov = [ zeros(Float64, aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L2, aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L2) for i in isingEnv.IOidl ]
sCor = Array{Float64}(undef, aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L2, aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L2)
@time mEvoFunc.getSpinStat!(isingEnv, aIsingData[iBatch].aLivingPop[iPop].aGty[1], aSpinAve, aSpinCov, sCor, isingDTMCprmStat)

my_cmap=get_cmap("YlOrRd"); my_cmap.set_under("Black"); my_cmap.set_over("White")

subplot(131)
imshow(GMat,cmap=my_cmap,vmin=aIsingData[iBatch].aLivingPop[iPop].aGty[1].gbounds[1],vmax=aIsingData[iBatch].aLivingPop[iPop].aGty[1].gbounds[2]);
title("Interaction Strength");
xticks(0:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L:2aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L-1, 1:Int32(aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].halfL):aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L);
yticks(0:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L:2aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L-1, 1:Int32(aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].halfL):aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L);
# colorbar();

subplot(132)
imshow(aSpinAve[1],cmap="Greys_r",vmin=-1,vmax=1);
title("Input Field: High");
xticks(0:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L-1, 1:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L);
yticks(0:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L-1, 1:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L);

subplot(133)
imshow(aSpinAve[2],cmap="Greys_r",vmin=-1,vmax=1);
title("Input Field: Low");
xticks(0:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L-1, 1:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L);
yticks(0:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L-1, 1:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L);

tight_layout();
# -


