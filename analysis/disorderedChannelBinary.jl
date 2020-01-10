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
const NEPOCHS, LOGFASMPL = 2, 1/5

# evolution and population constants
const NGEN, NPOP, NPOPNICHE = Int32(10), Int32(2*10^3), Int32(5)
const REPFACTOR, PNTMUTFACTOR, NMUTMAX = 10.0, 0.07, Int32(15)
const SCALINGFACTOR = 1

# metagenotypic constants
const SYSTEMSIZE, PVIABILITY, KOUT = Int32(6), 1.0, 100.0
const GMIN, GMAX = -2.0, 2.0

# environmental constants
const REPSTRENGTH, SELSTRENGTH, NIO = 1.0, 100.0, Int32(4);
const HF, LF = 1.0, 10.0^-4;

# +
# metagenotypic and genotypic types
aDisChnMetaGty = [ mEvoFunc.tDisChnMetaGty( SYSTEMSIZE, PVIABILITY, KOUT ) ]

aDisChnGty = [ tAlphaGty( aDisChnMetaGty, rand( [GMIN,GMAX], aDisChnMetaGty[1].dG ), [GMIN,GMAX] ) for i in 1:trunc(Int32,REPFACTOR)*NPOP ];
aDisChnGtyElite = [ tAlphaGty( aDisChnMetaGty, deepcopy(aDisChnGty[i].G), [GMIN,GMAX] ) for i in 1:trunc(Int32,REPFACTOR)*NPOP ];
aDisChnGtyNiche = [ tAlphaGty( aDisChnMetaGty, deepcopy(aDisChnGty[i].G), [GMIN,GMAX] ) for i in 1:NPOP ];
aDisChnGtyOne = [ tAlphaGty( aDisChnMetaGty, deepcopy(aDisChnGty[i].G), [GMIN,GMAX] ) for i in 1:NPOP ];
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
disChnEnv = tCompEnv(aIO, [REPSTRENGTH, SELSTRENGTH])
disChnEty = tPntMutEty(REPFACTOR,PNTMUTFACTOR,[ aDisChnMetaGty[1].dG ],NMUTMAX);

# print IOideal
for io in disChnEnv.IOidl
    println( replace(e -> e == 1.0 ? 1 : 0, io[1])," → ", replace( x -> x >= 1/SYSTEMSIZE ? 1 : 0, io[2]) )
end

# +
disChnPop = mEvoFunc.initEvoPop( NPOP,disChnEty,disChnEnv,aDisChnMetaGty,aDisChnGty );
disChnPopElite = mEvoFunc.initEvoPop( NPOP,disChnEty,disChnEnv,aDisChnMetaGty,aDisChnGtyElite );
disChnPopNiche = mEvoFunc.initEvoPop( NPOP,disChnEty,disChnEnv,aDisChnMetaGty,aDisChnGtyNiche );
disChnPopOne = mEvoFunc.initEvoPop( NPOP,disChnEty,disChnEnv,aDisChnMetaGty,aDisChnGtyOne );

aNichePop = [ mEvoFunc.tEvoPop( [NPOPNICHE,NPOPNICHE],disChnEty,disChnEnv,aDisChnMetaGty,[aDisChnGtyNiche[i] for j in 1:NPOPNICHE] ) for i in 1:NPOP ];
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
#     disChnEty = tPntMutEty(REPFACTOR,PNTMUTFACTOR,[ aDisChnMetaGty[1].dG ],NMUTMAX);

    disChnPop = mEvoFunc.initEvoPop( NPOP,disChnEty,disChnEnv,aDisChnMetaGty,aDisChnGty );
    disChnPopElite = mEvoFunc.initEvoPop( NPOP,disChnEty,disChnEnv,aDisChnMetaGty,aDisChnGtyElite );
    disChnPopNiche = mEvoFunc.initEvoPop( NPOP,disChnEty,disChnEnv,aDisChnMetaGty,aDisChnGtyNiche );
    disChnPopOne = mEvoFunc.initEvoPop( NPOP,disChnEty,disChnEnv,aDisChnMetaGty,aDisChnGtyOne );
    
    aNichePop = [ mEvoFunc.initEvoPop( NPOPNICHE,disChnEty,disChnEnv,aDisChnMetaGty,aNichePop[i].aGty ) for i in 1:NPOP ]

    push!( aDisChnData, tEvoData(NGEN, LOGFASMPL, 0.0) )
    push!( aDisChnDataElite, tEvoData(NGEN, LOGFASMPL, 0.0) )
    push!( aDisChnDataNiche, tEvoData(NGEN, LOGFASMPL, 0.0) )
    push!( aDisChnDataOne, tEvoData(NGEN, LOGFASMPL, 0.0) )

    mEvoFunc.gmsPopED!(disChnPop,aDisChnData[end],elite=false)
    mEvoFunc.gmsPopED!(disChnPopElite,aDisChnDataElite[end],elite=true)
    mEvoFunc.gmsNicED!(disChnPopNiche,aNichePop,aDisChnDataNiche[end],elite=false)
    mEvoFunc.gmsNicOneED!(disChnPopOne,aDisChnDataOne[end],elite=false)
end
# -

subplots(3,1,figsize=(7.5,7.5))
subplot(311); title("Average Performance: Loss");
    plot( vcat([dataBtc.avePerformance for dataBtc in aDisChnData]...) );
    plot( vcat([dataBtc.avePerformance for dataBtc in aDisChnDataElite]...) );
    plot( vcat([dataBtc.avePerformance for dataBtc in aDisChnDataNiche]...) );
    plot( vcat([dataBtc.avePerformance for dataBtc in aDisChnDataOne]...) );
subplot(312); title("Growth Factor");
    plot( vcat([dataBtc.growthFactor for dataBtc in aDisChnData]...) );
    plot( vcat([dataBtc.growthFactor for dataBtc in aDisChnDataElite]...) );
    plot( vcat([dataBtc.growthFactor for dataBtc in aDisChnDataNiche]...) );
    plot( vcat([dataBtc.growthFactor for dataBtc in aDisChnDataOne]...) );
subplot(313); title("Mutation Factor");
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

aDisChnMetaGty[1].dG, 4aDisChnMetaGty[1].L2mL

[ mEvoFunc.response(aDisChnGty[1], io[1]) for io in disChnEnv.IOidl ];

# ## Performance Analysis

Nbins=30
subplots(1,3,figsize=(15,5))
subplot(131); title("Population"); hist( [disChnPop.aGty[i].aF[1] for i in 1:disChnPop.pN[2]], Nbins, log=false, density=true, rwidth=.9);
subplot(132); title("Niche"); hist( [disChnPopNiche.aGty[i].aF[1] for i in 1:disChnPopNiche.pN[2]], Nbins, log=false, density=true, rwidth=.9);
subplot(133); title("One"); hist( [disChnPopOne.aGty[i].aF[1] for i in 1:disChnPopOne.pN[2]], Nbins, log=false, density=true, rwidth=.9);
suptitle("Loss distribution", fontsize=16);

Nbins=30
subplots(1,3,figsize=(15,5))
subplot(131); title("Population"); hist( [disChnPop.aGty[i].aF[1] for i in 1:disChnPop.pN[2]], Nbins, log=false, density=true, rwidth=.9);
subplot(132); title("Niche"); hist( [disChnPopNiche.aGty[i].aF[1] for i in 1:disChnPopNiche.pN[2]], Nbins, log=false, density=true, rwidth=.9);
subplot(133); title("One"); hist( [disChnPopOne.aGty[i].aF[1] for i in 1:disChnPopOne.pN[2]], Nbins, log=false, density=true, rwidth=.9);
suptitle("Loss distribution", fontsize=16);

Nbins=30
subplots(1,3,figsize=(15,5))
subplot(131); title("Population"); hist( [disChnPop.aGty[i].aF[1] for i in 1:disChnPop.pN[2]], Nbins, log=false, density=true, rwidth=.9);
subplot(132); title("Niche"); hist( [disChnPopNiche.aGty[i].aF[1] for i in 1:disChnPopNiche.pN[2]], Nbins, log=false, density=true, rwidth=.9);
subplot(133); title("One"); hist( [disChnPopOne.aGty[i].aF[1] for i in 1:disChnPopOne.pN[2]], Nbins, log=false, density=true, rwidth=.9);
suptitle("Loss distribution", fontsize=16);

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

r = aDisChnMetaGty[1].dG/NPOP

Gstat = tStat(aDisChnGty[iPop])
Rstat = tStat(aDisChnMetaGty[1].L2)

λval = [ λ for λ in mUtils.MPpdf(0,r)[2]:0.01:mUtils.MPpdf(0,r)[3] ]
MPval = [ mUtils.MPpdf(λ,r)[1] for λ in mUtils.MPpdf(0,r)[2]:0.01:mUtils.MPpdf(0,r)[3] ];

for iBatch in eachindex(aDisChnData), iPop in eachindex(aDisChnData[iBatch].aEvoPop)
    gen = aDisChnData[iBatch].aGen[iPop] + ( iBatch > 1 ? sum([ aDisChnData[jBatch].Ngen for jBatch in 1:iBatch-1]) : 0 )
    mEvoFunc.getGStat!(aDisChnData[iBatch].aEvoPop[iPop], Gstat)
    aResp = mEvoFunc.getRStat!(aDisChnData[iBatch].aEvoPop[iPop], Rstat)
    
    evGCor = eigvals( Symmetric(Gstat.cov) )
    evRCor = eigvals( Symmetric(Rstat.cov) )

    subplots(4,2,figsize=(7,14))
    subplot(421); title("G values histogram");  hist(vcat( [aDisChnData[iBatch].aEvoPop[iPop].aGty[i].G for i in 1:aDisChnData[iBatch].aEvoPop[iPop].pN[2]]... ), aDisChnData[iBatch].aEvoPop[iPop].aGty[1].pdg[1], density=true,rwidth=.9);
    subplot(422); title("G correlations"); imshow(Gstat.cor,cmap="viridis",vmin=-1,vmax=1); colorbar();
    subplot(412); title("G eigenvalue distribution"); hist(evGCor,30,density=true,log=false,rwidth=.9); plot(λval,MPval)
    
#     subplot(425); title("R values histogram"); hist( reshape(aResp,:,1), [10^i for i in -10:.5:0], log=true, density=true, rwidth=.9); xscale("log");
    subplot(425); title("F values histogram"); hist( [aDisChnData[iBatch].aEvoPop[iPop].aGty[i].aF[1] for i in 1:aDisChnData[iBatch].aEvoPop[iPop].pN[2]], 20, log=false, density=true, rwidth=.9);
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
subplot(421); title("G values histogram");  hist(vcat( [disChnPop.aGty[i].G for i in 1:disChnPop.pN[2]]... ), aDisChnData[Nbatch].aEvoPop[iPop].aGty[1].pdg[1], density=true,rwidth=.9);
subplot(422); title("G correlations"); imshow(Gstat.cor,cmap="viridis",vmin=-1,vmax=1); colorbar();
subplot(412); title("G eigenvalue distribution"); hist(evGCor,30,density=true,log=false,rwidth=.9); plot(λval,MPval)

# subplot(425); title("R values histogram"); hist( reshape(aResp,:,1), [10^i for i in -10:.5:0], log=true, density=true, rwidth=.9); xscale("log");
subplot(425); title("F values histogram"); hist( [disChnPop.aGty[i].aF[1] for i in 1:disChnPop.pN[2]], 60, log=false, density=true, rwidth=.9);
subplot(426); title("R correlations"); imshow(Rstat.cor,cmap="viridis",vmin=-1,vmax=1); colorbar();
xticks(disChnPop.aMetaGty[1].L-1:disChnPop.aMetaGty[1].L:disChnPop.aMetaGty[1].L2,disChnPop.aMetaGty[1].L:disChnPop.aMetaGty[1].L:disChnPop.aMetaGty[1].L2);
yticks(disChnPop.aMetaGty[1].L-1:disChnPop.aMetaGty[1].L:disChnPop.aMetaGty[1].L2,disChnPop.aMetaGty[1].L:disChnPop.aMetaGty[1].L:disChnPop.aMetaGty[1].L2);
subplot(414); title("R eigenvalue distribution"); hist(evRCor,30,density=true,log=false,rwidth=.9);

suptitle("Generation $gen", fontsize=16)
tight_layout(); subplots_adjust(top=0.9333);

# +
Nbatch, iPop, iGty = length(aDisChnDataElite), 1, 1

r = aDisChnMetaGty[1].dG/NPOP

Gstat = tStat(aDisChnGtyElite[iPop])
Rstat = tStat(aDisChnMetaGty[1].L2)

λval = [ λ for λ in mUtils.MPpdf(0,r)[2]:0.01:mUtils.MPpdf(0,r)[3] ]
MPval = [ mUtils.MPpdf(λ,r)[1] for λ in mUtils.MPpdf(0,r)[2]:0.01:mUtils.MPpdf(0,r)[3] ];

for iBatch in eachindex(aDisChnDataElite), iPop in eachindex(aDisChnDataElite[iBatch].aEvoPop)
    gen = aDisChnDataElite[iBatch].aGen[iPop] + ( iBatch > 1 ? sum([ aDisChnData[jBatch].Ngen for jBatch in 1:iBatch-1]) : 0 )
    mEvoFunc.getGStat!(aDisChnDataElite[iBatch].aEvoPop[iPop], Gstat)
    aResp = mEvoFunc.getRStat!(aDisChnDataElite[iBatch].aEvoPop[iPop], Rstat)
    
    evGCor = eigvals( Symmetric(Gstat.cov) )
    evRCor = eigvals( Symmetric(Rstat.cov) )

    subplots(4,2,figsize=(7,14))
    subplot(421); title("G values histogram");  hist(vcat( [aDisChnDataElite[iBatch].aEvoPop[iPop].aGty[i].G for i in 1:aDisChnDataElite[iBatch].aEvoPop[iPop].pN[2]]... ), aDisChnDataElite[iBatch].aEvoPop[iPop].aGty[1].pdg[1], density=true,rwidth=.9);
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
subplot(421); title("G values histogram");  hist(vcat( [disChnPopElite.aGty[i].G for i in 1:disChnPop.pN[2]]... ), aDisChnDataElite[Nbatch].aEvoPop[iPop].aGty[1].pdg[1], density=true,rwidth=.9);
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

r = aDisChnMetaGty[1].dG/NPOP

Gstat = tStat(aDisChnGtyNiche[iPop])
Rstat = tStat(aDisChnMetaGty[1].L2)

λval = [ λ for λ in mUtils.MPpdf(0,r)[2]:0.01:mUtils.MPpdf(0,r)[3] ]
MPval = [ mUtils.MPpdf(λ,r)[1] for λ in mUtils.MPpdf(0,r)[2]:0.01:mUtils.MPpdf(0,r)[3] ];

for iBatch in eachindex(aDisChnDataNiche), iPop in eachindex(aDisChnDataNiche[iBatch].aEvoPop)
    gen = aDisChnDataNiche[iBatch].aGen[iPop] + ( iBatch > 1 ? sum([ aDisChnDataNiche[jBatch].Ngen for jBatch in 1:iBatch-1]) : 0 )
    mEvoFunc.getGStat!(aDisChnDataNiche[iBatch].aEvoPop[iPop], Gstat)
    aResp = mEvoFunc.getRStat!(aDisChnDataNiche[iBatch].aEvoPop[iPop], Rstat)
    
    evGCor = eigvals( Symmetric(Gstat.cov) )
    evRCor = eigvals( Symmetric(Rstat.cov) )

    subplots(4,2,figsize=(7,14))
    subplot(421); title("G values histogram");  hist(vcat( [aDisChnDataNiche[iBatch].aEvoPop[iPop].aGty[i].G for i in 1:aDisChnDataNiche[iBatch].aEvoPop[iPop].pN[2]]... ), aDisChnDataNiche[iBatch].aEvoPop[iPop].aGty[1].pdg[1], density=true,rwidth=.9);
    subplot(422); title("G correlations"); imshow(Gstat.cor,cmap="viridis",vmin=-1,vmax=1); colorbar();
    subplot(412); title("G eigenvalue distribution"); hist(evGCor,30,density=true,log=false,rwidth=.9); plot(λval,MPval)
    
#     subplot(425); title("R values histogram"); hist( reshape(aResp,:,1), [10^i for i in -10:.5:0], log=true, density=true, rwidth=.9); xscale("log");
    subplot(425); title("F values histogram"); hist( [aDisChnDataNiche[iBatch].aEvoPop[iPop].aGty[i].aF[1] for i in 1:aDisChnDataNiche[iBatch].aEvoPop[iPop].pN[2]], 20, log=false, density=true, rwidth=.9);
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
subplot(421); title("G values histogram");  hist(vcat( [disChnPopNiche.aGty[i].G for i in 1:disChnPopNiche.pN[2]]... ), aDisChnDataNiche[Nbatch].aEvoPop[iPop].aGty[1].pdg[1], density=true,rwidth=.9);
subplot(422); title("G correlations"); imshow(Gstat.cor,cmap="viridis",vmin=-1,vmax=1); colorbar();
subplot(412); title("G eigenvalue distribution"); hist(evGCor,30,density=true,log=false,rwidth=.9); plot(λval,MPval)

# subplot(425); title("R values histogram"); hist( reshape(aResp,:,1), [10^i for i in -10:.5:0], log=true, density=true, rwidth=.9); xscale("log");
subplot(425); title("F values histogram"); hist( [disChnPopNiche.aGty[i].aF[1] for i in 1:disChnPopNiche.pN[2]], 20, log=false, density=true, rwidth=.9);
subplot(426); title("R correlations"); imshow(Rstat.cor,cmap="viridis",vmin=-1,vmax=1); colorbar();
xticks(disChnPopNiche.aMetaGty[1].L-1:disChnPopNiche.aMetaGty[1].L:disChnPopNiche.aMetaGty[1].L2,disChnPopNiche.aMetaGty[1].L:disChnPopNiche.aMetaGty[1].L:disChnPopNiche.aMetaGty[1].L2);
yticks(disChnPopNiche.aMetaGty[1].L-1:disChnPopNiche.aMetaGty[1].L:disChnPopNiche.aMetaGty[1].L2,disChnPopNiche.aMetaGty[1].L:disChnPopNiche.aMetaGty[1].L:disChnPopNiche.aMetaGty[1].L2);
subplot(414); title("R eigenvalue distribution"); hist(evRCor,30,density=true,log=false,rwidth=.9);

suptitle("Generation $gen", fontsize=16)
tight_layout(); subplots_adjust(top=0.9333);

# +
Nbatch, iPop, iGty = length(aDisChnDataOne), 1, 1

r = aDisChnMetaGty[1].dG/NPOP

Gstat = tStat(aDisChnGtyOne[iPop])
Rstat = tStat(aDisChnMetaGty[1].L2)

λval = [ λ for λ in mUtils.MPpdf(0,r)[2]:0.01:mUtils.MPpdf(0,r)[3] ]
MPval = [ mUtils.MPpdf(λ,r)[1] for λ in mUtils.MPpdf(0,r)[2]:0.01:mUtils.MPpdf(0,r)[3] ];

for iBatch in eachindex(aDisChnDataOne), iPop in eachindex(aDisChnDataOne[iBatch].aEvoPop)
    gen = aDisChnDataOne[iBatch].aGen[iPop] + ( iBatch > 1 ? sum([ aDisChnDataOne[jBatch].Ngen for jBatch in 1:iBatch-1]) : 0 )
    mEvoFunc.getGStat!(aDisChnDataOne[iBatch].aEvoPop[iPop], Gstat)
    aResp = mEvoFunc.getRStat!(aDisChnDataOne[iBatch].aEvoPop[iPop], Rstat)
    
    evGCor = eigvals( Symmetric(Gstat.cov) )
    evRCor = eigvals( Symmetric(Rstat.cov) )

    subplots(4,2,figsize=(7,14))
    subplot(421); title("G values histogram");  hist(vcat( [aDisChnDataOne[iBatch].aEvoPop[iPop].aGty[i].G for i in 1:aDisChnDataOne[iBatch].aEvoPop[iPop].pN[2]]... ), aDisChnDataOne[iBatch].aEvoPop[iPop].aGty[1].pdg[1], density=true,rwidth=.9);
    subplot(422); title("G correlations"); imshow(Gstat.cor,cmap="viridis",vmin=-1,vmax=1); colorbar();
    subplot(412); title("G eigenvalue distribution"); hist(evGCor,30,density=true,log=false,rwidth=.9); plot(λval,MPval)
    
#     subplot(425); title("R values histogram"); hist( reshape(aResp,:,1), [10^i for i in -10:.5:0], log=true, density=true, rwidth=.9); xscale("log");
    subplot(425); title("F values histogram"); hist( [aDisChnDataOne[iBatch].aEvoPop[iPop].aGty[i].aF[1] for i in 1:aDisChnDataOne[iBatch].aEvoPop[iPop].pN[2]], 20, log=false, density=true, rwidth=.9);
    subplot(426); title("R correlations"); imshow(Rstat.cor,cmap="viridis",vmin=-1,vmax=1); colorbar(); xticks(disChnPopOne.aMetaGty[1].L-1:disChnPopOne.aMetaGty[1].L:disChnPopOne.aMetaGty[1].L2,disChnPopOne.aMetaGty[1].L:disChnPopOne.aMetaGty[1].L:disChnPopOne.aMetaGty[1].L2); yticks(disChnPopOne.aMetaGty[1].L-1:disChnPopOne.aMetaGty[1].L:disChnPopOne.aMetaGty[1].L2,disChnPopOne.aMetaGty[1].L:disChnPopOne.aMetaGty[1].L:disChnPopOne.aMetaGty[1].L2);
    subplot(414); title("R eigenvalue distribution"); hist(evRCor,30,density=true,log=false,rwidth=.9);
    
    suptitle("Generation $gen", fontsize=16)
    tight_layout(); subplots_adjust(top=0.9333);
end

gen = sum([aDisChnDataOne[jBatch].Ngen for jBatch in 1:Nbatch])
mEvoFunc.getGStat!(disChnPopOne, Gstat)
aResp = mEvoFunc.getRStat!(disChnPopOne, Rstat)

evGCor = eigvals( Symmetric(Gstat.cov) )
evRCor = eigvals( Symmetric(Rstat.cov) )

subplots(4,2,figsize=(7,14))
subplot(421); title("G values histogram");  hist(vcat( [disChnPopOne.aGty[i].G for i in 1:disChnPopOne.pN[2]]... ), aDisChnDataOne[Nbatch].aEvoPop[iPop].aGty[1].pdg[1], density=true,rwidth=.9);
subplot(422); title("G correlations"); imshow(Gstat.cor,cmap="viridis",vmin=-1,vmax=1); colorbar();
subplot(412); title("G eigenvalue distribution"); hist(evGCor,30,density=true,log=false,rwidth=.9); plot(λval,MPval)

# subplot(425); title("R values histogram"); hist( reshape(aResp,:,1), [10^i for i in -10:.5:0], log=true, density=true, rwidth=.9); xscale("log");
subplot(425); title("F values histogram"); hist( [disChnPopOne.aGty[i].aF[1] for i in 1:disChnPopOne.pN[2]], 20, log=false, density=true, rwidth=.9);
subplot(426); title("R correlations"); imshow(Rstat.cor,cmap="viridis",vmin=-1,vmax=1); colorbar();
xticks(disChnPopOne.aMetaGty[1].L-1:disChnPopOne.aMetaGty[1].L:disChnPopOne.aMetaGty[1].L2,disChnPopOne.aMetaGty[1].L:disChnPopOne.aMetaGty[1].L:disChnPopOne.aMetaGty[1].L2);
yticks(disChnPopOne.aMetaGty[1].L-1:disChnPopOne.aMetaGty[1].L:disChnPopOne.aMetaGty[1].L2,disChnPopOne.aMetaGty[1].L:disChnPopOne.aMetaGty[1].L:disChnPopOne.aMetaGty[1].L2);
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
aOutputNodes = [ aDisChnMetaGty[1].L2mL .+ findall(x -> x >= 1/SYSTEMSIZE, disChnEnv.IOidl[i][2]) for i in eachindex(disChnEnv.IOidl) ];;

# +
# GKP Dynamics Statistics: last evolution step
fmax = maximum([ gty.aF[2] for gty in disChnPop.aGty[1:disChnPop.pN[2]] ])
iGty = findall( f -> f == fmax, [ gty.aF[2] for gty in disChnPop.aGty[1:disChnPop.pN[2]] ])[1]

fmin = minimum([ gty.aF[2] for gty in disChnPop.aGty[1:disChnPop.pN[2]] ])
iGtyLow = findall( f -> f == fmin, [ gty.aF[2] for gty in disChnPop.aGty[1:disChnPop.pN[2]] ])[1]

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

# fittestFluxPattern = tFluxPattern(disChnEnv,aDisChnData[Nbatch].aEvoPop[iPop].aGty[iGty])
# fittestCrntPtrn = tCrntPattern(disChnEnv,aDisChnData[Nbatch].aEvoPop[iPop].aGty[iGty])

# @time mEvoFunc.getFluxStat!(disChnEnv,aDisChnData[Nbatch].aEvoPop[iPop].aGty[1],fittestFluxPattern);
# @time mEvoFunc.getCrntPtrn!(aDisChnData[Nbatch].aEvoPop[iPop].aGty[1],fittestFluxPattern,fittestCrntPtrn);

# inputNodes = findall(x -> x == 1.0, disChnEnv.IOidl[io][1])
# outputNodes = aDisChnMetaGty[1].L2mL .+ findall(x -> x >= 1/SYSTEMSIZE, disChnEnv.IOidl[io][2]);

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
fmax = maximum([ gty.aF[2] for gty in disChnPopNiche.aGty[1:disChnPopNiche.pN[2]] ])
iGty = findall( f -> f == fmax, [ gty.aF[2] for gty in disChnPopNiche.aGty[1:disChnPopNiche.pN[2]] ])[1]

fmin = minimum([ gty.aF[2] for gty in disChnPopNiche.aGty[1:disChnPopNiche.pN[2]] ])
iGtyLow = findall( f -> f == fmin, [ gty.aF[2] for gty in disChnPopNiche.aGty[1:disChnPopNiche.pN[2]] ])[1]

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
aOutputNodes = [ aDisChnMetaGty[1].L2mL .+ findall(x -> x >= 1/SYSTEMSIZE, disChnEnv.IOidl[i][2]) for i in eachindex(disChnEnv.IOidl) ];;

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
plot( collect(-1:Δi:1), [ mEvoFunc.fitness(disChnEnv,tAlphaGty( [aDisChnMetaGty[1]], t*F.S[3]*F.Vt[3,:], [GMIN,GMAX] ))[1] for t in -1:Δi:1] );
