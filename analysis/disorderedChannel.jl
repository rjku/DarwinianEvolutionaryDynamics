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
const NGEN, NPOP = Int32(100), Int32(10^3)
const REPFACTOR, PNTMUTFACTOR, NMUTMAX = 5.0, 0.1, Int32(15)

# metagenotypic constants
const SYSTEMSIZE, PVIABILITY, KOUT = Int32(12), 1.0, 100.0
const GMIN, GMAX = -1.5, 1.5

# environmental constants
const SELSTRENGTH, NIO = 1.0, Int32(1);
const HIGHFLOW, LOWFLOW = 1.0, 10.0^-3;

# +
# genotypic and metagenotypic types
aDisChnMGty = [ mEvoFunc.tDisChnMGty( SYSTEMSIZE, PVIABILITY, KOUT ) ]

# uniformly distributed initial genotypic variables
aDisChnGty = [ tCntGty( [aDisChnMGty[1]], GMIN .+ rand( aDisChnMGty[1].dG ) .* (GMAX - GMIN), [GMIN,GMAX] ) for i in 1:REPFACTOR*NPOP ];

# normally distributed initial genotypic variables
# aDisChnGty = [ tCntGty( [aDisChnMGty[1]], randn( aDisChnMGty[1].dG ), [GMIN,GMAX] ) for i in 1:REPFACTOR*NPOP ];

# metagenotypes and genotypes from file
# aIsingMGty, aIsingGty, Npop = mEvoFunc.read_aIsingSigTransGty(isingDTMCprm,"test_#1")

# +
# environmental and evotypic types

disChnEnv = tCompEnv([ [ rand((LOWFLOW,HIGHFLOW),SYSTEMSIZE), rand((LOWFLOW,HIGHFLOW),SYSTEMSIZE) ] for i in 1:NIO ],SELSTRENGTH)
for iio in disChnEnv.IOidl iio[2] = iio[2] ./ sum(iio[2]) end
disChnEty = tPntMutEty(REPFACTOR,PNTMUTFACTOR,[ aDisChnMGty[1].dG ],NMUTMAX);

# +
# adding/deleting input--outputs to computational environment
# newOutputIdl = rand((LOWFLOW,HIGHFLOW),SYSTEMSIZE); push!(disChnEnv.IOidl,[ rand((LOWFLOW,HIGHFLOW),SYSTEMSIZE), newOutputIdl ./ sum(newOutputIdl) ]);
# deleteat!(disChnEnv.IOidl,<+iio+>)
# -

@time disChnPop = mEvoFunc.initLivingPop( NPOP,disChnEty,disChnEnv,aDisChnMGty,aDisChnGty );

aDisChnData = tEvoData[];

for i in 1:NBATCHES
	push!( aDisChnData, tEvoData(NGEN, LOGFASMPL) )
    mEvoFunc.evolutionGKP!(disChnPop,aDisChnData[end],ubermode=false)
end

subplots(3,1,figsize=(7.5,7.5))
subplot(311); title("Average Fitness"); plot( vcat([dataBtc.aveFitness for dataBtc in aDisChnData]...) );
subplot(312); title("Growth Factor"); plot( vcat([dataBtc.growthFactor for dataBtc in aDisChnData]...) );
subplot(313); title("Mutation Factor"); plot( vcat([dataBtc.mutationFactor for dataBtc in aDisChnData]...) );
tight_layout();

# ## System Properties

aDisChnMGty[1].dG, 4aDisChnMGty[1].L2mL

disChnEnv.IOidl

[ mEvoFunc.response(aDisChnGty[1], io[1]) for io in disChnEnv.IOidl ]

# ## Genotype Analysis

# +
Nbatch, iPop, iGty = length(aDisChnData), 1, 1

GAve = zeros(Float64, aDisChnData[1].aLivingPop[iPop].aGty[iGty].pMetaGty[1].dG)
GCov = zeros(Float64, aDisChnData[1].aLivingPop[iPop].aGty[iGty].pMetaGty[1].dG, aDisChnData[1].aLivingPop[iPop].aGty[iGty].pMetaGty[1].dG)
GCor = zeros(Float64, aDisChnData[1].aLivingPop[iPop].aGty[iGty].pMetaGty[1].dG, aDisChnData[1].aLivingPop[iPop].aGty[iGty].pMetaGty[1].dG)

for iBatch in eachindex(aDisChnData), iPop in eachindex(aDisChnData[iBatch].aLivingPop)
    gen = aDisChnData[iBatch].aGen[iPop] + ( iBatch > 1 ? sum([ aDisChnData[jBatch].Ngen for jBatch in 1:iBatch-1]) : 0 )
    (iBatch - 1)*( iBatch > 1 ? aDisChnData[iBatch - 1].Ngen : 1 )
    mEvoFunc.getGStat!(aDisChnData[iBatch].aLivingPop[iPop], GAve, GCov, GCor)
    evGCor = replace( e -> e < 10^-12 ? 0.0 : e, eigvals( Symmetric(GCov) ) )

    subplots(2,2,figsize=(6,6))
    subplot(221); title("G values histogram");  hist(vcat( [aDisChnData[iBatch].aLivingPop[iPop].aGty[i].G for i in 1:aDisChnData[iBatch].aLivingPop[iPop].pN[2]]... ), collect(aDisChnData[iBatch].aLivingPop[iPop].aGty[1].gbounds[1]:.1:aDisChnData[iBatch].aLivingPop[iPop].aGty[1].gbounds[2]), density=true,rwidth=.9);
    subplot(222); title("G correlations"); imshow(GCor,cmap="viridis",vmin=-1,vmax=1); colorbar();
    subplot(212); title("G eigenvalue distribution"); hist(evGCor,20,log=false,rwidth=.9); 
    suptitle("Generation $gen", fontsize=16)
    tight_layout(); subplots_adjust(top=0.85);
end

gen = sum([aDisChnData[jBatch].Ngen for jBatch in 1:Nbatch])
mEvoFunc.getGStat!(disChnPop, GAve, GCov, GCor)
evGCor = replace( e -> e < 10^-12 ? 0.0 : e, eigvals( Symmetric(GCov) ) )

subplots(2,2,figsize=(6,6))
subplot(221); title("G values histogram");  hist(vcat( [disChnPop.aGty[i].G for i in 1:disChnPop.pN[2]]... ), collect(disChnPop.aGty[1].gbounds[1]:.1:disChnPop.aGty[1].gbounds[2]), density=true,rwidth=.9);
subplot(222); title("G correlations"); imshow(GCor,cmap="viridis",vmin=-1,vmax=1); colorbar();
subplot(212); title("G eigenvalue distribution"); hist(evGCor,20,log=true,rwidth=.9); 
suptitle("Generation $gen", fontsize=16)
tight_layout(); subplots_adjust(top=0.85);
