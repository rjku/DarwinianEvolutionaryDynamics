# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.0
#   kernelspec:
#     display_name: Julia (2 threads) 1.3.1
#     language: julia
#     name: julia-(2-threads)-1.3
# ---

using Revise, BenchmarkTools, PyPlot, Distances
using Statistics, LinearAlgebra
using mEvoTypes
import mEvoFunc, mUtils

# +
# data constants
const NEPOCHS, LOGFASMPL = 3, 1/5

# evolution and population constants
const NGEN, NPOP, NPOPNICHE = Int32(10), Int32(10^3), Int32(10)
const REPFACTOR, PNTMUTFACTOR, NMUTMAX = 10.0, 0.07, Int32(15)
const SCALINGFACTOR = 1

# metagenotypic constants
const SYSTEMSIZE, PVIABILITY, KOUT = Int32(6), 1.0, 100.0
const GMIN, GMAX = -2.0, 2.0

# environmental constants
const SELSTRENGTH, NIO = 1.0, Int32(4);
const HF, LF = 1.0, 10.0^-4;

# +
# genotypic and metagenotypic types
aDisChnMGty = [ mEvoFunc.tDisChnMGty( SYSTEMSIZE, PVIABILITY, KOUT ) ]

# uniformly distributed initial genotypic variables
aDisChnGty = [ tAlphaGty( [aDisChnMGty[1]], rand( [GMIN,GMAX], aDisChnMGty[1].dG ), [GMIN,GMAX] ) for i in 1:REPFACTOR*NPOP ];

# normally distributed initial genotypic variables
# aDisChnGty = [ tCntGty( [aDisChnMGty[1]], randn( aDisChnMGty[1].dG ), [GMIN,GMAX] ) for i in 1:REPFACTOR*NPOP ];

# niched array of genotypes
aDisChnGtyElite = deepcopy(aDisChnGty[1:NPOP]);
aDisChnGtyNiche = deepcopy(aDisChnGty[1:NPOP]);
aDisChnGtyOne = deepcopy(aDisChnGty[1:NPOP]);

# metagenotypes and genotypes from file
# aIsingMGty, aIsingGty, Npop = mEvoFunc.read_aIsingSigTransGty(isingDTMCprm,"test_#1")
# -

# environmental and evotypic types: RANDOM IO
aIO = [ [ rand((LF,HF),SYSTEMSIZE), rand((LF,HF),SYSTEMSIZE) ] for i in 1:NIO ]
for io in aIO io[2] = io[2] ./ sum(io[2]) end

# adding random input--outputs to computational environment
newIO = [ rand((LOWFLOW,HIGHFLOW),SYSTEMSIZE), rand((LOWFLOW,HIGHFLOW),SYSTEMSIZE) ]
newIO[2] = newIO[2] ./ sum(newIO[2])
push!(disChnEnv.IOidl, newIO);

# environmental and evotypic types: DESIGNED BOOLEAN AND-GATE IO
aIO = [
    [ [HF, HF, HF, HF, HF, HF], [HF, HF, LF, LF, LF, LF] ],
    [ [HF, HF, HF, HF, LF, LF], [LF, LF, LF, LF, HF, HF] ],
    [ [HF, HF, LF, LF, HF, HF], [LF, LF, LF, LF, HF, HF] ],
    [ [HF, HF, LF, LF, LF, LF], [LF, LF, LF, LF, HF, HF] ]
]
for io in aIO io[2] = io[2] ./ sum(io[2]) end

# environmental and evotypic types: DESIGNED BOOLEAN NOT-GATE IO
aIO = [
    [ [LF, LF, HF, HF, LF, LF, LF, LF, LF, LF], [LF, LF,LF, LF, LF, LF, HF, HF, LF, LF] ],
    [ [LF, LF, LF, LF, LF, LF, HF, HF, LF, LF], [LF, LF,HF, HF, LF, LF, LF, LF, LF, LF] ]
]
for io in aIO io[2] = io[2] ./ sum(io[2]) end

# +
disChnEnv = tCompEnv(aIO, SELSTRENGTH)
disChnEty = tPntMutEty(REPFACTOR,PNTMUTFACTOR,[ aDisChnMGty[1].dG ],NMUTMAX);

# print IOideal
for io in disChnEnv.IOidl
    println( replace(e -> e == 1.0 ? 1 : 0, io[1])," → ", replace( x -> x >= 1/SYSTEMSIZE ? 1 : 0, io[2]) )
end

# +
disChnPop = mEvoFunc.initLivingPop( NPOP,disChnEty,disChnEnv,aDisChnMGty,aDisChnGty );
disChnPopElite = mEvoFunc.initLivingPop( NPOP,disChnEty,disChnEnv,aDisChnMGty,aDisChnGtyElite );
disChnPopNiche = mEvoFunc.initLivingPop( NPOP,disChnEty,disChnEnv,aDisChnMGty,aDisChnGtyNiche );
disChnPopOne = mEvoFunc.initLivingPop( NPOP,disChnEty,disChnEnv,aDisChnMGty,aDisChnGtyNiche );

aNichePop = [ mEvoFunc.initLivingPop( NPOPNICHE,disChnEty,disChnEnv,aDisChnMGty,[aDisChnGtyNiche[i] for j in 1:NPOPNICHE] ) for i in 1:NPOP ];
# -

aDisChnData = tEvoData[];
aDisChnDataElite = tEvoData[];
aDisChnDataNiche = tEvoData[];
aDisChnDataOne = tEvoData[];

# +
# simulated annealing with selection strength
for i in 1:NEPOCHS

#     selStrength = 10.0^(-i%3+1)
#     disChnEnv = tCompEnv(aIO, SELSTRENGTH)
#     disChnEty = tPntMutEty(REPFACTOR,PNTMUTFACTOR,[ aDisChnMGty[1].dG ],NMUTMAX);

    disChnPop = mEvoFunc.initLivingPop( NPOP,disChnEty,disChnEnv,aDisChnMGty,aDisChnGty );
    disChnPopElite = mEvoFunc.initLivingPop( NPOP,disChnEty,disChnEnv,aDisChnMGty,aDisChnGtyElite );
    disChnPopNiche = mEvoFunc.initLivingPop( NPOP,disChnEty,disChnEnv,aDisChnMGty,aDisChnGtyNiche );
    disChnPopOne = mEvoFunc.initLivingPop( NPOP,disChnEty,disChnEnv,aDisChnMGty,aDisChnGtyOne );
    
    aNichePop = [ mEvoFunc.initLivingPop( NPOPNICHE,disChnEty,disChnEnv,aDisChnMGty,aNichePop[i].aGty ) for i in 1:NPOP ]

    push!( aDisChnData, tEvoData(NGEN, LOGFASMPL) )
    push!( aDisChnDataElite, tEvoData(NGEN, LOGFASMPL) )
    push!( aDisChnDataNiche, tEvoData(NGEN, LOGFASMPL) )
    push!( aDisChnDataOne, tEvoData(NGEN, LOGFASMPL) )

    mEvoFunc.gmsPopED!(disChnPop,aDisChnData[end],elite=false)
    mEvoFunc.gmsPopED!(disChnPopElite,aDisChnDataElite[end],elite=true)
    mEvoFunc.gmsNicED!(disChnPopNiche,aNichePop,aDisChnDataNiche[end],elite=true)
    mEvoFunc.gmsNicOneED!(disChnPopOne,aDisChnDataOne[end])
end
# -

subplots(4,1,figsize=(7.5,10))
subplot(411); title("Average Fitness");
    plot( vcat([dataBtc.aveFitness for dataBtc in aDisChnData]...) );
    plot( vcat([dataBtc.aveFitness for dataBtc in aDisChnDataElite]...) );
    plot( vcat([dataBtc.aveFitness for dataBtc in aDisChnDataNiche]...) );
    plot( vcat([dataBtc.aveFitness for dataBtc in aDisChnDataOne]...) );
subplot(412); title("Growth Factor");
    plot( vcat([dataBtc.growthFactor for dataBtc in aDisChnData]...) );
    plot( vcat([dataBtc.growthFactor for dataBtc in aDisChnDataElite]...) );
    plot( vcat([dataBtc.growthFactor for dataBtc in aDisChnDataNiche]...) );
    plot( vcat([dataBtc.growthFactor for dataBtc in aDisChnDataOne]...) );
subplot(413); title("Average Teleonomy");
    plot( vcat([dataBtc.aveTeleonomy for dataBtc in aDisChnData]...) );
    plot( vcat([dataBtc.aveTeleonomy for dataBtc in aDisChnDataElite]...) );
    plot( vcat([dataBtc.aveTeleonomy for dataBtc in aDisChnDataNiche]...) );
    plot( vcat([dataBtc.aveTeleonomy for dataBtc in aDisChnDataOne]...) );
subplot(414); title("Mutation Factor");
    plot( vcat([dataBtc.mutationFactor for dataBtc in aDisChnData]...) );
    plot( vcat([dataBtc.mutationFactor for dataBtc in aDisChnDataElite]...) );
    plot( vcat([dataBtc.mutationFactor for dataBtc in aDisChnDataNiche]...) );
    plot( vcat([dataBtc.mutationFactor for dataBtc in aDisChnDataOne]...) );
tight_layout();

# +
Nsmpls = 100

ldf = REPFACTOR*mean(aDisChnData[end].aveFitness[end-Nsmpls:end]) - mean(aDisChnData[end].growthFactor[end-Nsmpls:end])
# -

# ## System Properties

aDisChnMGty[1].dG, 4aDisChnMGty[1].L2mL

[ mEvoFunc.response(aDisChnGty[1], io[1]) for io in disChnEnv.IOidl ];

# ## Genotype and Response Analysis
#
# Observations:
#
# * I/O shape the system. As more and more I/O are imposed, more and more genotypic variables become constrained, and their correlations become more and more polarized. How to quantify the Capacity: i.e. how many I/O can be encoded.
# * For the response variables, it is hard to distinguish correlations induced by the selection from the intrinsic correlations, as the latter are not trivial in the first place ... not to forget phylogeny.
# * Higher selection strength seems to polarize correlations, both genotypic and phenotypic ones. This occured even though the fitness is very low (achieved via selection strength quenching/sudden switch). In this regime, the distribution of response variables seems to be power law, and the eigenvalue distribution sparser.
# * Elite selection speeds up the evolution towards fitter individuals but reduces the "genetic variance" of the population. Indeed, the correlations of both genotypic and response variables become highly polarized, as for high selection strength. It seems like an analogy with capacity can be of help: think of differing genotpyes as liquid molecules in a vessel. When we squeeze the vessel some of the liquid get lost.

# +
Nbatch, iPop, iGty = length(aDisChnData), 1, 1

r = aDisChnMGty[1].dG/NPOP

Gstat = tStat(aDisChnGty[iPop])
Rstat = tStat(aDisChnMGty[1].L2)

λval = [ λ for λ in mUtils.MPpdf(0,r)[2]:0.01:mUtils.MPpdf(0,r)[3] ]
MPval = [ mUtils.MPpdf(λ,r)[1] for λ in mUtils.MPpdf(0,r)[2]:0.01:mUtils.MPpdf(0,r)[3] ];

for iBatch in eachindex(aDisChnData), iPop in eachindex(aDisChnData[iBatch].aLivingPop)
    gen = aDisChnData[iBatch].aGen[iPop] + ( iBatch > 1 ? sum([ aDisChnData[jBatch].Ngen for jBatch in 1:iBatch-1]) : 0 )
    mEvoFunc.getGStat!(aDisChnData[iBatch].aLivingPop[iPop], Gstat)
    aResp = mEvoFunc.getRStat!(aDisChnData[iBatch].aLivingPop[iPop], Rstat)
    
    evGCor = eigvals( Symmetric(Gstat.cov) )
    evRCor = eigvals( Symmetric(Rstat.cov) )

    subplots(4,2,figsize=(7,14))
    subplot(421); title("G values histogram");  hist(vcat( [aDisChnData[iBatch].aLivingPop[iPop].aGty[i].G for i in 1:aDisChnData[iBatch].aLivingPop[iPop].pN[2]]... ), aDisChnDataElite[iBatch].aLivingPop[iPop].aGty[1].pdg[1], density=true,rwidth=.9);
    subplot(422); title("G correlations"); imshow(Gstat.cor,cmap="viridis",vmin=-1,vmax=1); colorbar();
    subplot(412); title("G eigenvalue distribution"); hist(evGCor,30,density=true,log=false,rwidth=.9); plot(λval,MPval)
    
    subplot(425); title("R values histogram"); hist( reshape(aResp,:,1), [10^i for i in -10:.5:0], log=true, density=true, rwidth=.9); xscale("log");
    subplot(426); title("R correlations"); imshow(Rstat.cor,cmap="viridis",vmin=-1,vmax=1); colorbar(); xticks(disChnPop.aMetaGty[1].L-1:disChnPop.aMetaGty[1].L:disChnPop.aMetaGty[1].L2,disChnPop.aMetaGty[1].L:disChnPop.aMetaGty[1].L:disChnPop.aMetaGty[1].L2); yticks(disChnPop.aMetaGty[1].L-1:disChnPop.aMetaGty[1].L:disChnPop.aMetaGty[1].L2,disChnPop.aMetaGty[1].L:disChnPop.aMetaGty[1].L:disChnPop.aMetaGty[1].L2);
    subplot(414); title("R eigenvalue distribution"); hist(evRCor,30,density=true,log=false,rwidth=.9);
    
    suptitle("Generation $gen", fontsize=16)
    tight_layout(); subplots_adjust(top=0.9333);
end

gen = sum([aDisChnData[jBatch].Ngen for jBatch in 1:Nbatch])
mEvoFunc.getGStat!(disChnPop, Gstat)
aResp = mEvoFunc.getRStat!(disChnPop, Rstat)

evGCor = eigvals( Symmetric(Gstat.cov) )
evRCor = eigvals( Symmetric(Rstat.cov) )

subplots(4,2,figsize=(7,14))
subplot(421); title("G values histogram");  hist(vcat( [disChnPop.aGty[i].G for i in 1:disChnPop.pN[2]]... ), aDisChnDataElite[Nbatch].aLivingPop[iPop].aGty[1].pdg[1], density=true,rwidth=.9);
subplot(422); title("G correlations"); imshow(Gstat.cor,cmap="viridis",vmin=-1,vmax=1); colorbar();
subplot(412); title("G eigenvalue distribution"); hist(evGCor,30,density=true,log=false,rwidth=.9); plot(λval,MPval)

subplot(425); title("R values histogram"); hist( reshape(aResp,:,1), [10^i for i in -10:.5:0], log=true, density=true, rwidth=.9); xscale("log");
subplot(426); title("R correlations"); imshow(Rstat.cor,cmap="viridis",vmin=-1,vmax=1); colorbar();
xticks(disChnPop.aMetaGty[1].L-1:disChnPop.aMetaGty[1].L:disChnPop.aMetaGty[1].L2,disChnPop.aMetaGty[1].L:disChnPop.aMetaGty[1].L:disChnPop.aMetaGty[1].L2);
yticks(disChnPop.aMetaGty[1].L-1:disChnPop.aMetaGty[1].L:disChnPop.aMetaGty[1].L2,disChnPop.aMetaGty[1].L:disChnPop.aMetaGty[1].L:disChnPop.aMetaGty[1].L2);
subplot(414); title("R eigenvalue distribution"); hist(evRCor,30,density=true,log=false,rwidth=.9);

suptitle("Generation $gen", fontsize=16)
tight_layout(); subplots_adjust(top=0.9333);

# +
Nbatch, iPop, iGty = length(aDisChnDataElite), 1, 1

r = aDisChnMGty[1].dG/NPOP

Gstat = tStat(aDisChnGtyElite[iPop])
Rstat = tStat(aDisChnMGty[1].L2)

λval = [ λ for λ in mUtils.MPpdf(0,r)[2]:0.01:mUtils.MPpdf(0,r)[3] ]
MPval = [ mUtils.MPpdf(λ,r)[1] for λ in mUtils.MPpdf(0,r)[2]:0.01:mUtils.MPpdf(0,r)[3] ];

for iBatch in eachindex(aDisChnDataElite), iPop in eachindex(aDisChnDataElite[iBatch].aLivingPop)
    gen = aDisChnDataElite[iBatch].aGen[iPop] + ( iBatch > 1 ? sum([ aDisChnData[jBatch].Ngen for jBatch in 1:iBatch-1]) : 0 )
    mEvoFunc.getGStat!(aDisChnDataElite[iBatch].aLivingPop[iPop], Gstat)
    aResp = mEvoFunc.getRStat!(aDisChnDataElite[iBatch].aLivingPop[iPop], Rstat)
    
    evGCor = eigvals( Symmetric(Gstat.cov) )
    evRCor = eigvals( Symmetric(Rstat.cov) )

    subplots(4,2,figsize=(7,14))
    subplot(421); title("G values histogram");  hist(vcat( [aDisChnDataElite[iBatch].aLivingPop[iPop].aGty[i].G for i in 1:aDisChnDataElite[iBatch].aLivingPop[iPop].pN[2]]... ), aDisChnDataElite[iBatch].aLivingPop[iPop].aGty[1].pdg[1], density=true,rwidth=.9);
    subplot(422); title("G correlations"); imshow(Gstat.cor,cmap="viridis",vmin=-1,vmax=1); colorbar();
    subplot(412); title("G eigenvalue distribution"); hist(evGCor,30,density=true,log=false,rwidth=.9); plot(λval,MPval)
    
    subplot(425); title("R values histogram"); hist( reshape(aResp,:,1), [10^i for i in -10:.5:0], log=true, density=true, rwidth=.9); xscale("log");
    subplot(426); title("R correlations"); imshow(Rstat.cor,cmap="viridis",vmin=-1,vmax=1); colorbar(); xticks(disChnPopElite.aMetaGty[1].L-1:disChnPop.aMetaGty[1].L:disChnPop.aMetaGty[1].L2,disChnPop.aMetaGty[1].L:disChnPop.aMetaGty[1].L:disChnPop.aMetaGty[1].L2); yticks(disChnPop.aMetaGty[1].L-1:disChnPop.aMetaGty[1].L:disChnPop.aMetaGty[1].L2,disChnPop.aMetaGty[1].L:disChnPop.aMetaGty[1].L:disChnPop.aMetaGty[1].L2);
    subplot(414); title("R eigenvalue distribution"); hist(evRCor,30,density=true,log=false,rwidth=.9);
    
    suptitle("Generation $gen", fontsize=16)
    tight_layout(); subplots_adjust(top=0.9333);
end

gen = sum([aDisChnDataElite[jBatch].Ngen for jBatch in 1:Nbatch])
mEvoFunc.getGStat!(disChnPopElite, Gstat)
aResp = mEvoFunc.getRStat!(disChnPopElite, Rstat)

evGCor = eigvals( Symmetric(Gstat.cov) )
evRCor = eigvals( Symmetric(Rstat.cov) )

subplots(4,2,figsize=(7,14))
subplot(421); title("G values histogram");  hist(vcat( [disChnPopElite.aGty[i].G for i in 1:disChnPop.pN[2]]... ), aDisChnDataElite[Nbatch].aLivingPop[iPop].aGty[1].pdg[1], density=true,rwidth=.9);
subplot(422); title("G correlations"); imshow(Gstat.cor,cmap="viridis",vmin=-1,vmax=1); colorbar();
subplot(412); title("G eigenvalue distribution"); hist(evGCor,30,density=true,log=false,rwidth=.9); plot(λval,MPval)

subplot(425); title("R values histogram"); hist( reshape(aResp,:,1), [10^i for i in -10:.5:0], log=true, density=true, rwidth=.9); xscale("log");
subplot(426); title("R correlations"); imshow(Rstat.cor,cmap="viridis",vmin=-1,vmax=1); colorbar();
xticks(disChnPopElite.aMetaGty[1].L-1:disChnPop.aMetaGty[1].L:disChnPop.aMetaGty[1].L2,disChnPop.aMetaGty[1].L:disChnPop.aMetaGty[1].L:disChnPop.aMetaGty[1].L2);
yticks(disChnPopElite.aMetaGty[1].L-1:disChnPop.aMetaGty[1].L:disChnPop.aMetaGty[1].L2,disChnPop.aMetaGty[1].L:disChnPop.aMetaGty[1].L:disChnPop.aMetaGty[1].L2);
subplot(414); title("R eigenvalue distribution"); hist(evRCor,30,density=true,log=false,rwidth=.9);

suptitle("Generation $gen", fontsize=16)
tight_layout(); subplots_adjust(top=0.9333);

# +
Nbatch, iPop, iGty = length(aDisChnDataNiche), 1, 1

r = aDisChnMGty[1].dG/NPOP

Gstat = tStat(aDisChnGtyNiche[iPop])
Rstat = tStat(aDisChnMGty[1].L2)

λval = [ λ for λ in mUtils.MPpdf(0,r)[2]:0.01:mUtils.MPpdf(0,r)[3] ]
MPval = [ mUtils.MPpdf(λ,r)[1] for λ in mUtils.MPpdf(0,r)[2]:0.01:mUtils.MPpdf(0,r)[3] ];

for iBatch in eachindex(aDisChnDataNiche), iPop in eachindex(aDisChnDataNiche[iBatch].aLivingPop)
    gen = aDisChnDataNiche[iBatch].aGen[iPop] + ( iBatch > 1 ? sum([ aDisChnDataNiche[jBatch].Ngen for jBatch in 1:iBatch-1]) : 0 )
    mEvoFunc.getGStat!(aDisChnDataNiche[iBatch].aLivingPop[iPop], Gstat)
    aResp = mEvoFunc.getRStat!(aDisChnDataNiche[iBatch].aLivingPop[iPop], Rstat)
    
    evGCor = eigvals( Symmetric(Gstat.cov) )
    evRCor = eigvals( Symmetric(Rstat.cov) )

    subplots(4,2,figsize=(7,14))
    subplot(421); title("G values histogram");  hist(vcat( [aDisChnDataNiche[iBatch].aLivingPop[iPop].aGty[i].G for i in 1:aDisChnDataNiche[iBatch].aLivingPop[iPop].pN[2]]... ), aDisChnDataNiche[iBatch].aLivingPop[iPop].aGty[1].pdg[1], density=true,rwidth=.9);
    subplot(422); title("G correlations"); imshow(Gstat.cor,cmap="viridis",vmin=-1,vmax=1); colorbar();
    subplot(412); title("G eigenvalue distribution"); hist(evGCor,30,density=true,log=false,rwidth=.9); plot(λval,MPval)
    
    subplot(425); title("R values histogram"); hist( reshape(aResp,:,1), [10^i for i in -10:.5:0], log=true, density=true, rwidth=.9); xscale("log");
    subplot(426); title("R correlations"); imshow(Rstat.cor,cmap="viridis",vmin=-1,vmax=1); colorbar(); xticks(disChnPopNiche.aMetaGty[1].L-1:disChnPopNiche.aMetaGty[1].L:disChnPopNiche.aMetaGty[1].L2,disChnPopNiche.aMetaGty[1].L:disChnPopNiche.aMetaGty[1].L:disChnPopNiche.aMetaGty[1].L2); yticks(disChnPopNiche.aMetaGty[1].L-1:disChnPopNiche.aMetaGty[1].L:disChnPopNiche.aMetaGty[1].L2,disChnPopNiche.aMetaGty[1].L:disChnPopNiche.aMetaGty[1].L:disChnPopNiche.aMetaGty[1].L2);
    subplot(414); title("R eigenvalue distribution"); hist(evRCor,30,density=true,log=false,rwidth=.9);
    
    suptitle("Generation $gen", fontsize=16)
    tight_layout(); subplots_adjust(top=0.9333);
end

gen = sum([aDisChnDataNiche[jBatch].Ngen for jBatch in 1:Nbatch])
mEvoFunc.getGStat!(disChnPopNiche, Gstat)
aResp = mEvoFunc.getRStat!(disChnPopNiche, Rstat)

evGCor = eigvals( Symmetric(Gstat.cov) )
evRCor = eigvals( Symmetric(Rstat.cov) )

subplots(4,2,figsize=(7,14))
subplot(421); title("G values histogram");  hist(vcat( [disChnPopNiche.aGty[i].G for i in 1:disChnPopNiche.pN[2]]... ), aDisChnDataNiche[Nbatch].aLivingPop[iPop].aGty[1].pdg[1], density=true,rwidth=.9);
subplot(422); title("G correlations"); imshow(Gstat.cor,cmap="viridis",vmin=-1,vmax=1); colorbar();
subplot(412); title("G eigenvalue distribution"); hist(evGCor,30,density=true,log=false,rwidth=.9); plot(λval,MPval)

subplot(425); title("R values histogram"); hist( reshape(aResp,:,1), [10^i for i in -10:.5:0], log=true, density=true, rwidth=.9); xscale("log");
subplot(426); title("R correlations"); imshow(Rstat.cor,cmap="viridis",vmin=-1,vmax=1); colorbar();
xticks(disChnPopNiche.aMetaGty[1].L-1:disChnPopNiche.aMetaGty[1].L:disChnPopNiche.aMetaGty[1].L2,disChnPopNiche.aMetaGty[1].L:disChnPopNiche.aMetaGty[1].L:disChnPopNiche.aMetaGty[1].L2);
yticks(disChnPopNiche.aMetaGty[1].L-1:disChnPopNiche.aMetaGty[1].L:disChnPopNiche.aMetaGty[1].L2,disChnPopNiche.aMetaGty[1].L:disChnPopNiche.aMetaGty[1].L:disChnPopNiche.aMetaGty[1].L2);
subplot(414); title("R eigenvalue distribution"); hist(evRCor,30,density=true,log=false,rwidth=.9);

suptitle("Generation $gen", fontsize=16)
tight_layout(); subplots_adjust(top=0.9333);

# +
θvec = collect(10^-3:10^-3:maximum(Gstat.cor))
cum = [ mUtils.cumulativeCount(θ, broadcast( e -> abs(e), Gstat.cor )) for θ in θvec]

title("Cumulative Count of High Correlated Genes")
plot(θvec,cum); xscale("log"); xlabel("θ"); yscale("log"); ylabel("count");
# -

# ## Response Analysis

aInputNodes = [ findall(x -> x == 1.0, disChnEnv.IOidl[i][1]) for i in eachindex(disChnEnv.IOidl) ]
aOutputNodes = [ aDisChnMGty[1].L2mL .+ findall(x -> x >= 1/SYSTEMSIZE, disChnEnv.IOidl[i][2]) for i in eachindex(disChnEnv.IOidl) ];;

# +
# GKP Dynamics Statistics: last evolution step
fmax = maximum([ gty.pF[1] for gty in disChnPop.aGty[1:disChnPop.pN[2]] ])
iGty = findall( f -> f == fmax, [ gty.pF[1] for gty in disChnPop.aGty[1:disChnPop.pN[2]] ])[1]

fmin = minimum([ gty.pF[1] for gty in disChnPop.aGty[1:disChnPop.pN[2]] ])
iGtyLow = findall( f -> f == fmin, [ gty.pF[1] for gty in disChnPop.aGty[1:disChnPop.pN[2]] ])[1]

fittestFluxPattern = tFluxPattern(disChnEnv,disChnPop.aGty[iGty])
fittestCrntPtrn = tCrntPattern(disChnEnv,disChnPop.aGty[iGty])

@time mEvoFunc.getFluxStat!(disChnEnv,disChnPop.aGty[iGty],fittestFluxPattern);
@time mEvoFunc.getCrntPtrn!(disChnPop.aGty[iGty],fittestFluxPattern,fittestCrntPtrn);

shittyFluxPattern = tFluxPattern(disChnEnv,disChnPop.aGty[iGtyLow])
shittyCrntPtrn = tCrntPattern(disChnEnv,disChnPop.aGty[iGtyLow])

@time mEvoFunc.getFluxStat!(disChnEnv,disChnPop.aGty[iGtyLow],shittyFluxPattern);
@time mEvoFunc.getCrntPtrn!(disChnPop.aGty[iGtyLow],shittyFluxPattern,shittyCrntPtrn);
# -

mat"""
    figure('Name','fittestEvolved')
"""
for io in eachindex(disChnEnv.IOidl)
    mat"""
        G = digraph($(fittestCrntPtrn.aV[io][2]), $(fittestCrntPtrn.aV[io][1]), $(fittestCrntPtrn.Jave[io]));
        LWidths = ( 10*G.Edges.Weight/max(G.Edges.Weight) );
    
        ioInside = double($io);

        subplot(2,2,ioInside)
        pl = plot(G,'Layout','force','LineWidth',LWidths,'NodeLabel',{});

        highlight(pl,$(aInputNodes[io]),'NodeColor','g');
        highlight(pl,$(aOutputNodes[io]),'NodeColor','r');
    """
end
#plot(G,'Layout','force','EdgeLabel',G.Edges.Weight,'LineWidth',$(LWidths))

mat"""
    figure('Name','shitty')
"""
for io in eachindex(disChnEnv.IOidl)
    mat"""
        G = digraph($(shittyCrntPtrn.aV[io][2]), $(shittyCrntPtrn.aV[io][1]), $(shittyCrntPtrn.Jave[io]));
        LWidths = ( 10*G.Edges.Weight/max(G.Edges.Weight) );
    
        ioInside = double($io);

        subplot(2,2,ioInside)
        pl = plot(G,'Layout','force','LineWidth',LWidths,'NodeLabel',{});

        highlight(pl,$(aInputNodes[io]),'NodeColor','g');
        highlight(pl,$(aOutputNodes[io]),'NodeColor','r');
    """
end
#plot(G,'Layout','force','EdgeLabel',G.Edges.Weight,'LineWidth',$(LWidths))

# +
# within batch
# Nbatch, iPop, Nio = length(aDisChnData), 1, length(disChnEnv.IOidl)

# fittestFluxPattern = tFluxPattern(disChnEnv,aDisChnData[Nbatch].aLivingPop[iPop].aGty[iGty])
# fittestCrntPtrn = tCrntPattern(disChnEnv,aDisChnData[Nbatch].aLivingPop[iPop].aGty[iGty])

# @time mEvoFunc.getFluxStat!(disChnEnv,aDisChnData[Nbatch].aLivingPop[iPop].aGty[1],fittestFluxPattern);
# @time mEvoFunc.getCrntPtrn!(aDisChnData[Nbatch].aLivingPop[iPop].aGty[1],fittestFluxPattern,fittestCrntPtrn);

# inputNodes = findall(x -> x == 1.0, disChnEnv.IOidl[io][1])
# outputNodes = aDisChnMGty[1].L2mL .+ findall(x -> x >= 1/SYSTEMSIZE, disChnEnv.IOidl[io][2]);

# +
# Elite Statistics: last evolution step
fittestFluxPatternElite = tFluxPattern(disChnEnv,disChnPopElite.aGty[1])
fittestCrntPtrnElite = tCrntPattern(disChnEnv,disChnPopElite.aGty[1])

@time mEvoFunc.getFluxStat!(disChnEnv,disChnPopElite.aGty[1],fittestFluxPatternElite);
@time mEvoFunc.getCrntPtrn!(disChnPopElite.aGty[1],fittestFluxPatternElite,fittestCrntPtrnElite);

shittyFluxPatternElite = tFluxPattern(disChnEnv,disChnPop.aGty[disChnPop.pN[2]])
shittyCrntPtrnElite = tCrntPattern(disChnEnv,disChnPop.aGty[disChnPop.pN[2]])

@time mEvoFunc.getFluxStat!(disChnEnv,disChnPop.aGty[disChnPop.pN[2]],shittyFluxPatternElite);
@time mEvoFunc.getCrntPtrn!(disChnPop.aGty[disChnPop.pN[2]],shittyFluxPattern,shittyCrntPtrnElite);
# -

mat"""
    figure('Name','fittestEvolvedElite')
"""
for io in eachindex(disChnEnv.IOidl)
    mat"""
        G = digraph($(fittestCrntPtrnElite.aV[io][2]), $(fittestCrntPtrnElite.aV[io][1]), $(fittestCrntPtrnElite.Jave[io]));
        LWidths = ( 10*G.Edges.Weight/max(G.Edges.Weight) );
    
        ioInside = double($io);

        subplot(2,2,ioInside)
        pl = plot(G,'Layout','force','LineWidth',LWidths,'NodeLabel',{});

        highlight(pl,$(aInputNodes[io]),'NodeColor','g');
        highlight(pl,$(aOutputNodes[io]),'NodeColor','r');
    """
end
#plot(G,'Layout','force','EdgeLabel',G.Edges.Weight,'LineWidth',$(LWidths))

# +
# Niche Statistics: last evolution step
fmax = maximum([ gty.pF[1] for gty in disChnPopNiche.aGty[1:disChnPopNiche.pN[2]] ])
iGty = findall( f -> f == fmax, [ gty.pF[1] for gty in disChnPopNiche.aGty[1:disChnPopNiche.pN[2]] ])[1]

fmin = minimum([ gty.pF[1] for gty in disChnPopNiche.aGty[1:disChnPopNiche.pN[2]] ])
iGtyLow = findall( f -> f == fmin, [ gty.pF[1] for gty in disChnPopNiche.aGty[1:disChnPopNiche.pN[2]] ])[1]

fittestFluxPatternNiche = tFluxPattern(disChnEnv,disChnPopNiche.aGty[iGty])
fittestCrntPtrnNiche = tCrntPattern(disChnEnv,disChnPopNiche.aGty[iGty])

@time mEvoFunc.getFluxStat!(disChnEnv,disChnPopNiche.aGty[iGty],fittestFluxPatternNiche);
@time mEvoFunc.getCrntPtrn!(disChnPopNiche.aGty[iGty],fittestFluxPatternNiche,fittestCrntPtrnNiche);
# -

mat"""
    figure('Name','fittestEvolvedNiche')
"""
for io in eachindex(disChnEnv.IOidl)
    mat"""
        G = digraph($(fittestCrntPtrnNiche.aV[io][2]), $(fittestCrntPtrnNiche.aV[io][1]), $(fittestCrntPtrnNiche.Jave[io]));
        LWidths = ( 10*G.Edges.Weight/max(G.Edges.Weight) );
    
        ioInside = double($io);

        subplot(2,2,ioInside)
        pl = plot(G,'Layout','force','LineWidth',LWidths,'NodeLabel',{});

        highlight(pl,$(aInputNodes[io]),'NodeColor','g');
        highlight(pl,$(aOutputNodes[io]),'NodeColor','r');
    """
end
#plot(G,'Layout','force','EdgeLabel',G.Edges.Weight,'LineWidth',$(LWidths))

mat"""
    figure('Name','shitty')
"""
for io in eachindex(disChnEnv.IOidl)
    mat"""
        G = digraph($(shittyCrntPtrn.aV[io][2]), $(shittyCrntPtrn.aV[io][1]), $(shittyCrntPtrn.Jave[io]));
        LWidths = ( 10*G.Edges.Weight/max(G.Edges.Weight) );
    
        ioInside = double($io);

        subplot(2,2,ioInside)
        pl = plot(G,'Layout','force','LineWidth',LWidths,'NodeLabel',{});

        highlight(pl,$(aInputNodes[io]),'NodeColor','g');
        highlight(pl,$(aOutputNodes[io]),'NodeColor','r');
    """
end
#plot(G,'Layout','force','EdgeLabel',G.Edges.Weight,'LineWidth',$(LWidths))

aInputNodes = [ findall(x -> x == 1.0, disChnEnv.IOidl[i][1]) for i in eachindex(disChnEnv.IOidl) ]
aOutputNodes = [ aDisChnMGty[1].L2mL .+ findall(x -> x >= 1/SYSTEMSIZE, disChnEnv.IOidl[i][2]) for i in eachindex(disChnEnv.IOidl) ];;

mat"""
    figure('Name','shitty')
"""
for io in eachindex(disChnEnv.IOidl)
    mat"""
        G = digraph($(shittyCrntPtrn.aV[io][2]), $(shittyCrntPtrn.aV[io][1]), $(shittyCrntPtrn.Jave[io]));
        LWidths = ( 10*G.Edges.Weight/max(G.Edges.Weight) );
    
        ioInside = double($io);

        subplot(2,2,ioInside)
        pl = plot(G,'Layout','force','LineWidth',LWidths,'NodeLabel',{});

        highlight(pl,$(aInputNodes[io]),'NodeColor','g');
        highlight(pl,$(aOutputNodes[io]),'NodeColor','r');
    """
end
#plot(G,'Layout','force','EdgeLabel',G.Edges.Weight,'LineWidth',$(LWidths))

mat"""
    figure('Name','fittestEvolved')
"""
for io in eachindex(disChnEnv.IOidl)
    mat"""
        G = digraph($(fittestCrntPtrn.aV[io][2]), $(fittestCrntPtrn.aV[io][1]), $(fittestCrntPtrn.Jave[io]));
        LWidths = ( 10*G.Edges.Weight/max(G.Edges.Weight) );
    
        ioInside = double($io);

        subplot(2,2,ioInside)
        pl = plot(G,'Layout','force','LineWidth',LWidths,'NodeLabel',{});

        highlight(pl,$(aInputNodes[io]),'NodeColor','g');
        highlight(pl,$(aOutputNodes[io]),'NodeColor','r');
    """
end
#plot(G,'Layout','force','EdgeLabel',G.Edges.Weight,'LineWidth',$(LWidths))

# ## Principal Component Analysis

Gmat = [ disChnPopNiche.aGty[i].G[g] for i in eachindex(disChnPopNiche.aGty), g in eachindex(disChnPopNiche.aGty[1].G) ];
GmatAve = [ Gmat[i,g] - Gstat.ave[g] for i in eachindex(disChnPopNiche.aGty), g in eachindex(disChnPopNiche.aGty[1].G) ];

F = svd(GmatAve);

aEv = broadcast( e -> e^2/(disChnPopNiche.pN[2]-1), F.S);
cumEnergy( p ) = sum(aEv[1:p])/sum(aEv);

plot([cumEnergy(p) for p in 1:length(aEv)]);

Δi = 0.001
plot( collect(-1:Δi:1), [ mEvoFunc.fitness(disChnEnv,tAlphaGty( [aDisChnMGty[1]], t*F.S[3]*F.Vt[3,:], [GMIN,GMAX] ))[1] for t in -1:Δi:1] );
