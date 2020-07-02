# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.1
#   kernelspec:
#     display_name: Julia (4 threads) 1.4.0
#     language: julia
#     name: julia-(4-threads)-1.4
# ---

# +
using Revise, BenchmarkTools, PyPlot, Distances, Printf
using Statistics, LinearAlgebra
using mEvoTypes
import mEvoFunc, mUtils, mGraphs

Base.show(io::IO, f::Float64) = @printf(io, "%.2f", f)

# +
# evolution and population constants
const NGEN, NPOP, NPOPNICHE = Int32(3*10^3), Int32(10^4), Int32(5)
const REPFACTOR, MUTFACTOR = 0.0, 0.1

# genotypic variables
GRIDSIZE = 5
DIMGSPACE = GRIDSIZE^2
BLOCKSIZE = GRIDSIZE÷3

# environmental constants
const REPSTRENGTH, SELSTRENGTH, MINREPCOEF = 0.0, 0.8, 1.0
const HF, MF, LF = 10.0, 1.0, 0.1;
# -

# ### Potential Shaping

f = [ LF for i in 1:DIMGSPACE ]
for i in eachindex(f)
    # top left block
    if (i-1)÷GRIDSIZE < BLOCKSIZE && (i-1)%GRIDSIZE < BLOCKSIZE
        f[i] = HF
    end
    
    # bottom right line
    if ( (i-1)÷GRIDSIZE == GRIDSIZE - BLOCKSIZE && (i-1)%GRIDSIZE > GRIDSIZE - BLOCKSIZE ) ||
        ( (i-1)÷GRIDSIZE >= GRIDSIZE - BLOCKSIZE && (i-1)%GRIDSIZE == GRIDSIZE - BLOCKSIZE )
        f[i] = HF
    end
    
    # bottom left point
    if (i-1)÷GRIDSIZE == BLOCKSIZE - 1 && (i-1)%GRIDSIZE == GRIDSIZE - BLOCKSIZE
        f[i] = HF
    end
    
    # bottom left line
    if GRIDSIZE - BLOCKSIZE <= (i-1)÷GRIDSIZE <= GRIDSIZE - BLOCKSIZE + 1 && (i-1)%GRIDSIZE == BLOCKSIZE - 1
        f[i] = HF
    end
    
    # top right block
    # if (i-1)÷GRIDSIZE >= GRIDSIZE - BLOCKSIZE && (i-1)%GRIDSIZE < BLOCKSIZE
    #     f[i] = MF
    # end
end

f = [ LF for i in 1:DIMGSPACE ]
for i in eachindex(f)
    # top left block
    if (i-1)÷GRIDSIZE % 2 == 1 && (i-1)%GRIDSIZE % 2 == 1
        f[i] = HF
    end
end

# +
grid = mGraphs.EdgeWeightedSquareLattice(GRIDSIZE,[ MUTFACTOR for i in 1:2GRIDSIZE*(GRIDSIZE-1) ]);

env = TabularEnvironment([f],[REPSTRENGTH, SELSTRENGTH])
ety = TabularEvotype([REPFACTOR],MINREPCOEF,grid);

# +
aS = [ exp(f[g]*env.aSelCoef[2]) for g in 1:DIMGSPACE ]
apS = aS ./ sum(aS)

aP = ( ( Array(mGraphs.transitionMatrix(ety.G)) .+ Diagonal([ 1.0 + MINREPCOEF for i in 1:DIMGSPACE ]) ) * aS ) .* aS
apP = aP ./ sum(aP)

gf = mGraphs.VertexWeightedSquareLattice( GRIDSIZE, f )
fMat = mGraphs.matrixForm(gf);

gpS = mGraphs.VertexWeightedSquareLattice( GRIDSIZE, apS )
pSMat = mGraphs.matrixForm(gpS);

gpP = mGraphs.VertexWeightedSquareLattice( GRIDSIZE, apP )
pPMat = mGraphs.matrixForm(gpP);
# -

aS[7], aS[1] 

subplots(1,3,figsize=(15,7))
subplot(131); title("Bare Fitness Function"); imshow(fMat,cmap="viridis");
subplot(132); title("Selectivity"); imshow(log.(pSMat),cmap="cividis",vmax=0);
subplot(133); title("Evolutionary Potential"); imshow(log.(pPMat),cmap="cividis",vmax=0);

subplots(1,3,figsize=(15,7))
subplot(131); title("Bare Fitness Function"); imshow(fMat,cmap="viridis");
subplot(132); title("Selectivity"); imshow(pSMat,cmap="cividis",vmax=1,vmin=0);
subplot(133); title("Evolutionary Potential"); imshow(pPMat,cmap="cividis",vmax=1,vmin=0);

# ### Population Dynamics

# genotypic types
aGty = [ TabularGenotype( rand( 1:GRIDSIZE^2, 1 ) ) for i in 1:NPOP ];

# +
pop = mEvoFunc.init_Population( NPOP,ety,env,aGty );

evoData = EvoData(NGEN);

mEvoFunc.gmsNicOneED!(pop,evoData,elite=false);
# -

subplots(3,1,figsize=(7.5,7.5))
subplot(311); title("Average Performance: Loss");
    plot( evoData.avePerformance );
subplot(312); title("Growth Factor");
    plot( evoData.growthFactor );
subplot(313); title("Mutation Factor");
    plot( evoData.mutationFactor );
tight_layout();

p = mGraphs.VertexWeightedSquareLattice( GRIDSIZE, mEvoFunc.popStatistics(pop) )
pMat = mGraphs.matrixForm(p);

subplots(1,3,figsize=(15,7))
subplot(131); title("Population Distribution"); imshow(log.(pMat),cmap="cividis",vmax=0); colorbar(shrink=0.5)
subplot(132); title("Selectivity"); imshow(log.(pSMat),cmap="cividis",vmax=0); colorbar(shrink=0.5)
subplot(133); title("Evolutionary Potential"); imshow(log.(pPMat),cmap="cividis",vmax=0); colorbar(shrink=0.5)

subplots(1,3,figsize=(15,7))
subplot(131); title("Population Distribution"); imshow(pMat,cmap="cividis",vmax=.1,vmin=0);
subplot(132); title("Selectivity"); imshow(pSMat,cmap="cividis",vmax=.1,vmin=0);
subplot(133); title("Evolutionary Potential"); imshow(pPMat,cmap="cividis",vmax=.1,vmin=0);

# ### Phase Diagram

# +
# evolution and population constants
const NGEN, NPOP, NPOPNICHE = Int32(3*10^3), Int32(10^4), Int32(5)
const REPFACTOR, MUTFACTOR = 0.0, 0.25

# genotypic variables
GRIDSIZE = 12
DIMGSPACE = GRIDSIZE^2
BLOCKSIZE = GRIDSIZE÷3

# environmental constants
const REPSTRENGTH, SELSTRENGTH, MINREPCOEF = 0.0, 1.0, 4.0
const HF, MF, LF = 10.0, 1.0, 0.1;

f = [ LF for i in 1:DIMGSPACE ]
for i in eachindex(f)
    # top left block
    if (i-1)÷GRIDSIZE < BLOCKSIZE && (i-1)%GRIDSIZE < BLOCKSIZE
        f[i] = HF
    end
    
    # bottom right line
    if ( (i-1)÷GRIDSIZE == GRIDSIZE - BLOCKSIZE && (i-1)%GRIDSIZE > GRIDSIZE - BLOCKSIZE ) ||
        ( (i-1)÷GRIDSIZE >= GRIDSIZE - BLOCKSIZE && (i-1)%GRIDSIZE == GRIDSIZE - BLOCKSIZE )
        f[i] = HF
    end
    
    # bottom left point
    if (i-1)÷GRIDSIZE == BLOCKSIZE - 2 && (i-1)%GRIDSIZE == GRIDSIZE - BLOCKSIZE + 1
        f[i] = HF
    end
    
    # top right point
#     if GRIDSIZE - BLOCKSIZE <= (i-1)÷GRIDSIZE <= GRIDSIZE - BLOCKSIZE + 1 && (i-1)%GRIDSIZE == BLOCKSIZE - 1
#         f[i] = HF
#     end
    
    # top right block
    if (i-1)÷GRIDSIZE >= GRIDSIZE - BLOCKSIZE && (i-1)%GRIDSIZE < BLOCKSIZE
        f[i] = MF
    end
end
# -

mutFactorVals = 0.24 * [ 10^-i for i in 0:1/2:1];
selStrengthVals = [ 10^-i for i in 0:1/2:1];


# +
aEvoData = EvoData[]
aaGty = Vector{Vector{TabularGenotype}}(undef,0)
aPop = Population[]

for β in selStrengthVals, μ in mutFactorVals
    grid = mGraphs.EdgeWeightedSquareLattice(GRIDSIZE,[ μ for i in 1:2GRIDSIZE*(GRIDSIZE-1) ])
    env = TabularEnvironment([f],[REPSTRENGTH, β])
    ety = TabularEvotype([REPFACTOR],MINREPCOEF,grid)

    push!(aEvoData,EvoData(NGEN))
    push!( aaGty, [ TabularGenotype( rand( 1:DIMGSPACE, 1 ) ) for i in 1:NPOP ] );
    push!( aPop, mEvoFunc.init_Population( NPOP,ety,env,aaGty[end] ) );
    
    mEvoFunc.gmsNicOneED!(aPop[end],aEvoData[end],elite=false);
end
# -

[ 10^(i-1) for i in 3/4:1/4:1 ]

# +
aPop[3].aGty .= aPop[6].aGty;

for β in [ 10^(i-1) for i in 3/4:1/4:1 ]
    env = TabularEnvironment([f],[REPSTRENGTH, β])
    mEvoFunc.gmsNicOneED!(aPop[3],aEvoData[3],elite=false);
end

# +
iPop = 8;

subplots(3,1,figsize=(7.5,7.5))
subplot(311); title("Average Performance: Loss");
    plot( aEvoData[iPop].avePerformance );
subplot(312); title("Growth Factor");
    plot( aEvoData[iPop].growthFactor );
subplot(313); title("Mutation Factor");
    plot( aEvoData[iPop].mutationFactor );
tight_layout();

# +
for pop in aPop
    P = mEvoFunc.popStatistics(pop)
    p = mGraphs.VertexWeightedSquareLattice( GRIDSIZE, P )
    pMat = mGraphs.matrixForm(p);
    
    aS = [ exp(f[g]*pop.env.aSelCoef[2]) for g in 1:DIMGSPACE ]
    apS = aS ./ sum(aS)

    aP = ( ( Array(mGraphs.transitionMatrix(pop.ety.G)) .+ Diagonal([ 1.0 + MINREPCOEF for i in 1:DIMGSPACE ]) ) * aS
        ) .* aS
    apP = aP ./ sum(aP)

    gpS = mGraphs.VertexWeightedSquareLattice( GRIDSIZE, apS )
    pSMat = mGraphs.matrixForm(gpS);

    gpP = mGraphs.VertexWeightedSquareLattice( GRIDSIZE, apP )
    pPMat = mGraphs.matrixForm(gpP);
    
    subplots(1,3,figsize=(15,7))
    subplot(131); title("Genotype PMF: β = $(pop.env.aSelCoef[2]), μ = $(pop.ety.G.We[1])");
        imshow(log.(pMat),cmap="cividis",vmax=0);
    subplot(132); title("Selectivity | KL = $(mUtils.KLdivergence(P,apS))");
        imshow(log.(pSMat),cmap="cividis",vmax=0);
    subplot(133); title("Evolutionary Potential | KL = $(mUtils.KLdivergence(P,apP))");
        imshow(log.(pPMat),cmap="cividis",vmax=0);
    
#         imshow(pMat,cmap="cividis",vmax=.15,vmin=0); # colorbar()
#     subplot(132); title("Selectivity"); imshow(pSMat,cmap="cividis",vmax=.15,vmin=0);
#     subplot(133); title("Evolutionary Potential"); imshow(pPMat,cmap="cividis",vmax=.15,vmin=0);
end
# -

