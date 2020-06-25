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

# +
# data constants
const NBATCHES, LOGFASMPL = 1, 1/2

# evolution and population constants
const NGEN, NPOP = Int32(5), Int32(5)
const REPFACTOR, PNTMUTFACTOR, NMUTMAX = 10.0, 0.07, Int32(15)

# metagenotypic constants
const SYSTEMSIZE, MAXSIZE = Int32(6), Int32(8)
const INVTEMPERATURE, EXTERNALFIELD = 0.8, 0.1
const GMIN, GMAX = -2.0, 2.0
# const DELTAG, DELTATOFFSET = 1/3, 0.01

# environmental constants
const SELSTRENGTH, NIO = 1.0, Int32(2);
const HIGHFLOW, LOWFLOW = 1.0, 10.0^-3;

# fast-time-scale dynamics constants
const NSAMPLINGS, NMCSPS, NTRIALS = Int32(20), Int32(10^5), Int32(1);

# -

isingEnv = tCompEnv( [[-1.0,-1.0], [1,1.0]], SELSTRENGTH )
isingEty = tPntMutEty(REPFACTOR, PNTMUTFACTOR, [ 2L^2 for L in SYSTEMSIZE:2:MAXSIZE ], NMUTMAX);
isingDTMCprm = tDTMCprm( NSAMPLINGS, NMCSPS, NTRIALS );

aIsingMetaGty = [ tIsingSigTransMetaGty(SYSTEMSIZE,INVTEMPERATURE,EXTERNALFIELD,isingDTMCprm) ];
aIsingGty = [ tCntGty( [aIsingMetaGty[1]], GMIN .+ rand(2SYSTEMSIZE^2) .* (GMAX - GMIN), [GMIN,GMAX] ) for i in 1:REPFACTOR*NPOP ];

@time isingPop = mEvoFunc.initEvoPop( NPOP,isingEty,isingEnv,aIsingMetaGty,aIsingGty );

aIsingData = tEvoData[];

for i in 1:NBATCHES
	push!( aIsingData, tEvoData(NGEN, LOGFASMPL) )
    mEvoFunc.gmsPopEDup!(isingPop,aIsingData[end],MAXSIZE,elite=true)
end

subplots(3,1,figsize=(7.5,7.5))
subplot(311); title("Average Fitness"); plot( vcat([dataBtc.aveFitness for dataBtc in aIsingData]...) );
subplot(312); title("Growth Factor"); plot( vcat([dataBtc.growthFactor for dataBtc in aIsingData]...) );
subplot(313); title("Mutation Factor"); plot( vcat([dataBtc.mutationFactor for dataBtc in aIsingData]...) );
tight_layout();


# +
iBatch, iPop, iIsing = 2, 2, 1

const NSMPLSTAT, NMCSPSSTAT, NTRLSTAT = Int32(2+10^1), Int32(10^5), Int32(2)
isingDTMCprmStat = tDTMCprm( NSMPLSTAT, NMCSPSSTAT, NTRLSTAT )

GMat = Array{Float64}(undef, 2aIsingData[iBatch].aEvoPop[iPop].aGty[1].pMetaGty[1].L, 2aIsingData[iBatch].aEvoPop[iPop].aGty[1].pMetaGty[1].L)
mEvoFunc.showG!(aIsingData[iBatch].aEvoPop[iPop].aGty[1], GMat)

aSpinAve = [ zeros(Float64, aIsingData[iBatch].aEvoPop[iPop].aGty[1].pMetaGty[1].L, aIsingData[iBatch].aEvoPop[iPop].aGty[1].pMetaGty[1].L) for i in isingEnv.IOidl ]
aSpinCov = [ zeros(Float64, aIsingData[iBatch].aEvoPop[iPop].aGty[1].pMetaGty[1].L2, aIsingData[iBatch].aEvoPop[iPop].aGty[1].pMetaGty[1].L2) for i in isingEnv.IOidl ]
sCor = Array{Float64}(undef, aIsingData[iBatch].aEvoPop[iPop].aGty[1].pMetaGty[1].L2, aIsingData[iBatch].aEvoPop[iPop].aGty[1].pMetaGty[1].L2)
@time mEvoFunc.getSpinStat!(isingEnv, aIsingData[iBatch].aEvoPop[iPop].aGty[1], aSpinAve, aSpinCov, sCor, isingDTMCprmStat)

my_cmap=get_cmap("YlOrRd"); my_cmap.set_under("Black"); my_cmap.set_over("White")

subplot(131)
imshow(GMat,cmap=my_cmap,vmin=aIsingData[iBatch].aEvoPop[iPop].aGty[1].gbounds[1],vmax=aIsingData[iBatch].aEvoPop[iPop].aGty[1].gbounds[2]);
title("Interaction Strength");
xticks(0:aIsingData[iBatch].aEvoPop[iPop].aGty[1].pMetaGty[1].L:2aIsingData[iBatch].aEvoPop[iPop].aGty[1].pMetaGty[1].L-1, 1:Int32(aIsingData[iBatch].aEvoPop[iPop].aGty[1].pMetaGty[1].halfL):aIsingData[iBatch].aEvoPop[iPop].aGty[1].pMetaGty[1].L);
yticks(0:aIsingData[iBatch].aEvoPop[iPop].aGty[1].pMetaGty[1].L:2aIsingData[iBatch].aEvoPop[iPop].aGty[1].pMetaGty[1].L-1, 1:Int32(aIsingData[iBatch].aEvoPop[iPop].aGty[1].pMetaGty[1].halfL):aIsingData[iBatch].aEvoPop[iPop].aGty[1].pMetaGty[1].L);
# colorbar();

subplot(132)
imshow(aSpinAve[1],cmap="Greys_r",vmin=-1,vmax=1);
title("Input Field: High");
xticks(0:aIsingData[iBatch].aEvoPop[iPop].aGty[1].pMetaGty[1].L-1, 1:aIsingData[iBatch].aEvoPop[iPop].aGty[1].pMetaGty[1].L);
yticks(0:aIsingData[iBatch].aEvoPop[iPop].aGty[1].pMetaGty[1].L-1, 1:aIsingData[iBatch].aEvoPop[iPop].aGty[1].pMetaGty[1].L);

subplot(133)
imshow(aSpinAve[2],cmap="Greys_r",vmin=-1,vmax=1);
title("Input Field: Low");
xticks(0:aIsingData[iBatch].aEvoPop[iPop].aGty[1].pMetaGty[1].L-1, 1:aIsingData[iBatch].aEvoPop[iPop].aGty[1].pMetaGty[1].L);
yticks(0:aIsingData[iBatch].aEvoPop[iPop].aGty[1].pMetaGty[1].L-1, 1:aIsingData[iBatch].aEvoPop[iPop].aGty[1].pMetaGty[1].L);

tight_layout();
# -


