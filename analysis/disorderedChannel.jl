# -*- coding: utf-8 -*-
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

using Revise, BenchmarkTools, PyPlot, MATLAB
using Statistics, LinearAlgebra
using mEvoTypes
import mEvoFunc, mUtils

# +
# data constants
const NEPOCHS, LOGFASMPL = 1, 1/5

# evolution and population constants
const NGEN, NPOP = Int32(500), Int32(2*10^3)
const REPFACTOR, PNTMUTFACTOR, NMUTMAX = 5.0, 0.07, Int32(15)

# metagenotypic constants
const SYSTEMSIZE, PVIABILITY, KOUT = Int32(6), 1.0, 100.0
const GMIN, GMAX = -1.5, 1.5

# environmental constants
const SELSTRENGTH, NIO = 1.0, Int32(4);
const HIGHFLOW, LOWFLOW = 1.0, 10.0^-4;

# +
# genotypic and metagenotypic types
aDisChnMGty = [ mEvoFunc.tDisChnMGty( SYSTEMSIZE, PVIABILITY, KOUT ) ]

# uniformly distributed initial genotypic variables
# aDisChnGty = [ tCntGty( [aDisChnMGty[1]], GMIN .+ rand( aDisChnMGty[1].dG ) .* (GMAX - GMIN), [GMIN,GMAX] ) for i in 1:REPFACTOR*NPOP ];

# normally distributed initial genotypic variables
aDisChnGty = [ tCntGty( [aDisChnMGty[1]], randn( aDisChnMGty[1].dG ), [GMIN,GMAX] ) for i in 1:REPFACTOR*NPOP ];

# metagenotypes and genotypes from file
# aIsingMGty, aIsingGty, Npop = mEvoFunc.read_aIsingSigTransGty(isingDTMCprm,"test_#1")

# +
# environmental and evotypic types
aIO = [ [ rand((LOWFLOW,HIGHFLOW),SYSTEMSIZE), rand((LOWFLOW,HIGHFLOW),SYSTEMSIZE) ] for i in 1:NIO ]
for io in aIO io[2] = io[2] ./ sum(io[2]) end

disChnEnv = tCompEnv(aIO, SELSTRENGTH)
disChnEty = tPntMutEty(REPFACTOR,PNTMUTFACTOR,[ aDisChnMGty[1].dG ],NMUTMAX);

for io in disChnEnv.IOidl
    println( replace(e -> e == 1.0 ? 1 : 0, io[1])," → ", replace( x -> x >= 1/SYSTEMSIZE ? 1 : 0, io[2]) )
end

# +
# environmental and evotypic types: (DESIGNED/NONRANDOM) BOOLEAN AND-GATE IO
aIO = [
    [ [HIGHFLOW, HIGHFLOW, HIGHFLOW, HIGHFLOW, HIGHFLOW, HIGHFLOW], [HIGHFLOW, HIGHFLOW, LOWFLOW, LOWFLOW, LOWFLOW, LOWFLOW] ],
    [ [HIGHFLOW, HIGHFLOW, HIGHFLOW, HIGHFLOW, LOWFLOW, LOWFLOW], [LOWFLOW, LOWFLOW, LOWFLOW, LOWFLOW, HIGHFLOW, HIGHFLOW] ],
    [ [LOWFLOW, LOWFLOW, HIGHFLOW, HIGHFLOW, HIGHFLOW, HIGHFLOW], [LOWFLOW, LOWFLOW, LOWFLOW, LOWFLOW, HIGHFLOW, HIGHFLOW] ],
    [ [LOWFLOW, LOWFLOW, HIGHFLOW, HIGHFLOW, LOWFLOW, LOWFLOW], [LOWFLOW, LOWFLOW, LOWFLOW, LOWFLOW, HIGHFLOW, HIGHFLOW] ]
]
for io in aIO io[2] = io[2] ./ sum(io[2]) end

disChnEnv = tCompEnv(aIO, SELSTRENGTH)
disChnEty = tPntMutEty(REPFACTOR,PNTMUTFACTOR,[ aDisChnMGty[1].dG ],NMUTMAX);

for io in disChnEnv.IOidl
    println( replace(e -> e == 1.0 ? 1 : 0, io[1])," → ", replace( x -> x >= 1/SYSTEMSIZE ? 1 : 0, io[2]) )
end
# -

@time disChnPop = mEvoFunc.initLivingPop( NPOP,disChnEty,disChnEnv,aDisChnMGty,aDisChnGty );

aDisChnData = tEvoData[];

# +
# adding input--outputs to computational environment
newIO = [ rand((LOWFLOW,HIGHFLOW),SYSTEMSIZE), rand((LOWFLOW,HIGHFLOW),SYSTEMSIZE) ]
newIO[2] = newIO[2] ./ sum(newIO[2])
push!(disChnEnv.IOidl, newIO);

# deleting input--outputs to computational environment
# deleteat!(disChnEnv.IOidl,<+iio+>)

# reinitializing population with the new environment
disChnPop = mEvoFunc.initLivingPop( NPOP,disChnEty,disChnEnv,aDisChnMGty,aDisChnGty );

for io in disChnEnv.IOidl
    println( replace(e -> e == 1.0 ? 1 : 0, io[1])," → ", replace( x -> x >= 1/SYSTEMSIZE ? 1 : 0, io[2]) )
end
# -

# new selection strength
disChnEnv = tCompEnv(aIO, SELSTRENGTH);
disChnPop = mEvoFunc.initLivingPop( NPOP,disChnEty,disChnEnv,aDisChnMGty,aDisChnGty );

for i in 1:NEPOCHS
	push!( aDisChnData, tEvoData(NGEN, LOGFASMPL) )
    mEvoFunc.evolutionGKP!(disChnPop,aDisChnData[end],elite=false)
end

subplots(3,1,figsize=(7.5,7.5))
subplot(311); title("Average Fitness"); plot( vcat([dataBtc.aveFitness for dataBtc in aDisChnData]...) );
subplot(312); title("Growth Factor"); plot( vcat([dataBtc.growthFactor for dataBtc in aDisChnData]...) );
subplot(313); title("Mutation Factor"); plot( vcat([dataBtc.mutationFactor for dataBtc in aDisChnData]...) );
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
    subplot(421); title("G values histogram");  hist(vcat( [aDisChnData[iBatch].aLivingPop[iPop].aGty[i].G for i in 1:aDisChnData[iBatch].aLivingPop[iPop].pN[2]]... ), collect(aDisChnData[iBatch].aLivingPop[iPop].aGty[1].gbounds[1]:.1:aDisChnData[iBatch].aLivingPop[iPop].aGty[1].gbounds[2]), density=true,rwidth=.9);
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
subplot(421); title("G values histogram");  hist(vcat( [disChnPop.aGty[i].G for i in 1:disChnPop.pN[2]]... ), collect(disChnPop.aGty[1].gbounds[1]:.1:disChnPop.aGty[1].gbounds[2]), density=true,rwidth=.9);
subplot(422); title("G correlations"); imshow(Gstat.cor,cmap="viridis",vmin=-1,vmax=1); colorbar();
subplot(412); title("G eigenvalue distribution"); hist(evGCor,30,density=true,log=false,rwidth=.9); plot(λval,MPval)

subplot(425); title("R values histogram"); hist( reshape(aResp,:,1), [10^i for i in -10:.5:0], log=true, density=true, rwidth=.9); xscale("log");
subplot(426); title("R correlations"); imshow(Rstat.cor,cmap="viridis",vmin=-1,vmax=1); colorbar();
xticks(disChnPop.aMetaGty[1].L-1:disChnPop.aMetaGty[1].L:disChnPop.aMetaGty[1].L2,disChnPop.aMetaGty[1].L:disChnPop.aMetaGty[1].L:disChnPop.aMetaGty[1].L2);
yticks(disChnPop.aMetaGty[1].L-1:disChnPop.aMetaGty[1].L:disChnPop.aMetaGty[1].L2,disChnPop.aMetaGty[1].L:disChnPop.aMetaGty[1].L:disChnPop.aMetaGty[1].L2);
subplot(414); title("R eigenvalue distribution"); hist(evRCor,30,density=true,log=false,rwidth=.9);

suptitle("Generation $gen", fontsize=16)
tight_layout(); subplots_adjust(top=0.9333);
# -

# ## Response Analysis

# +
Nbatch, iPop, Nio = length(aDisChnData), 1, length(disChnEnv.IOidl)

fmax = maximum([ gty.pF[1] for gty in disChnPop.aGty[1:disChnPop.pN[2]] ])
iGty = findall( f -> f == fmax, [ gty.pF[1] for gty in disChnPop.aGty[1:disChnPop.pN[2]] ])[1]

fmin = minimum([ gty.pF[1] for gty in disChnPop.aGty[1:disChnPop.pN[2]] ])
iGtyLow = findall( f -> f == fmin, [ gty.pF[1] for gty in disChnPop.aGty[1:disChnPop.pN[2]] ])[1]

# within batch
# fittestFluxPattern = tFluxPattern(disChnEnv,aDisChnData[Nbatch].aLivingPop[iPop].aGty[iGty])
# fittestCrntPtrn = tCrntPattern(disChnEnv,aDisChnData[Nbatch].aLivingPop[iPop].aGty[iGty])

# @time mEvoFunc.getFluxStat!(disChnEnv,aDisChnData[Nbatch].aLivingPop[iPop].aGty[1],fittestFluxPattern);
# @time mEvoFunc.getCrntPtrn!(aDisChnData[Nbatch].aLivingPop[iPop].aGty[1],fittestFluxPattern,fittestCrntPtrn);

# inputNodes = findall(x -> x == 1.0, disChnEnv.IOidl[io][1])
# outputNodes = aDisChnMGty[1].L2mL .+ findall(x -> x >= 1/SYSTEMSIZE, disChnEnv.IOidl[io][2]);

# last evolution step
fittestFluxPattern = tFluxPattern(disChnEnv,disChnPop.aGty[iGty])
fittestCrntPtrn = tCrntPattern(disChnEnv,disChnPop.aGty[iGty])

@time mEvoFunc.getFluxStat!(disChnEnv,disChnPop.aGty[iGty],fittestFluxPattern);
@time mEvoFunc.getCrntPtrn!(disChnPop.aGty[iGty],fittestFluxPattern,fittestCrntPtrn);

aInputNodes = [ findall(x -> x == 1.0, disChnEnv.IOidl[i][1]) for i in eachindex(disChnEnv.IOidl) ]
aOutputNodes = [ aDisChnMGty[1].L2mL .+ findall(x -> x >= 1/SYSTEMSIZE, disChnEnv.IOidl[i][2]) for i in eachindex(disChnEnv.IOidl) ];;
# -

for io in 1:Nio
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

