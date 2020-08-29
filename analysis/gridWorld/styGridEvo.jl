# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.5.2
#   kernelspec:
#     display_name: Julia (4 threads) 1.5.0
#     language: julia
#     name: julia-(4-threads)-1.5
# ---

# +
using Revise, BenchmarkTools, PyPlot, Printf, ProgressMeter
import EvolutionaryDynamics, mUtils, mPlot, mGraphs

# Base.show(io::IO, f::Float64) = @printf(io, "%.2f", f)

# +
folderName = "/Users/riccardorao/projects/fisheria/analysis/gridWorld/scratch/sty/"

# jobID = "styGridEvo" # null ID for checking parameters
jobID = "styGridEvo-7785591"  # phase diagram results
# jobID = "styGridEvo-7785776"  # 

include(folderName * jobID * "-parameters.jl")

# aTraj = mUtils.readJLD(folderName * jobID * "_aTraj.jld", "aTraj");

# +
fMat = mGraphs.matrixForm( mGraphs.VertexWeightedSquareLattice( GRIDSIZE, fTbl ) )

# imshow(fMat,cmap="cividis"); # origin="lower"
imshow(fMat,cmap="inferno"); # origin="lower"

xticks( prepend!(collect(3:4:GRIDSIZE-1),0), prepend!(collect(4:4:GRIDSIZE),1) )
yticks( prepend!(collect(3:4:GRIDSIZE-1),0), prepend!(collect(4:4:GRIDSIZE),1) )
colorbar()

# savefig("styLnF.pdf", bbox_inches="tight")
# -

# ## Data Generation/Collection

# +
aTraj = Array{EvolutionaryDynamics.TrajectoryData}(undef, length(aSelStrength), length(aMutFactor))

@time for (i,β) in enumerate(aSelStrength), (j,μ) in enumerate(aMutFactor)
    grid = mGraphs.EdgeWeightedSquareLattice(GRIDSIZE, [ μ for i in 1:2GRIDSIZE*(GRIDSIZE-1) ]);
    ety = EvolutionaryDynamics.TabularEvotype(
#         EvolutionaryDynamics.NeutralReplication(REPCOEF),
        EvolutionaryDynamics.WithoutReplication(),
        EvolutionaryDynamics.StandardMutation(),
        EvolutionaryDynamics.FitnessSelection(),
#         EvolutionaryDynamics.ElitismSelection(),
        grid
    )
    
#     traj = @distributed (+) for i in 1:NSAMPLES
    traj = EvolutionaryDynamics.generateTabularSystemsTrajectories(
            ety, fTbl, β, NPOP, NGENRELAX, NGENSAMPLE
        )
#     end
    aTraj[i,j] = traj
end

# +
using JLD, HDF5, DelimitedFiles

aTraj = mUtils.readJLD(folderName * jobID * "_aTraj.jld", "aTraj");
# -

# ## Analysis

# #### Checks about relaxation

# +
i = 9

subplots(2,1,figsize=(7.5,7.5))
subplot(211); title("Average Fitness");
    plot( aTraj[i].avePerformance );
subplot(212); title("Mutation Factor");
    plot( aTraj[i].mutationFactor );
tight_layout();
# -

# #### Genotype Distributions

# +
Base.show(io::IO, f::Float64) = @printf(io, "%.2f", f)

figNum = figure(figsize=(12,12))

for (i,β) in enumerate(aSelStrength), (j,μ) in enumerate(aMutFactor)

    grid = mGraphs.EdgeWeightedSquareLattice(GRIDSIZE, [ μ for i in 1:2GRIDSIZE*(GRIDSIZE-1) ]);
    Pgty = [ sum( aTraj[i,j].jointProb[g,:] ) for g in 1:size(aTraj[i,j].jointProb)[1] ] / aTraj[i,j].NgenSample

    p = mGraphs.VertexWeightedSquareLattice( GRIDSIZE, Pgty )
    pMat = mGraphs.matrixForm(p);

    ax = figNum.add_subplot(3,3,j + 3*(3-i))
    

    eMin = -10
    for (i,e) in enumerate(pMat)
        if e == 0
           pMat[i] = 10.0^eMin
        end
    end
    
    imshow(log10.(pMat),cmap="inferno",vmax=0, vmin=eMin);
#     title("β = $(β), μ = $(μ)");

    xticks([])
    yticks([])

    if i == 1 xticks( prepend!(collect(3:4:GRIDSIZE-1),0), prepend!(collect(4:4:GRIDSIZE),1) ) end
    if j == 1 yticks( prepend!(collect(3:4:GRIDSIZE-1),0), prepend!(collect(4:4:GRIDSIZE),1) ) end
    if i == 1 && j == 3 colorbar(fraction=0.045) end
end

savefig("pmfPhaseDiagram.pdf", bbox_inches="tight")
# +
figNum = figure(figsize=(12,12))

i, j = 2, 2
β, μ = aSelStrength[i], aMutFactor[j]
    
grid = mGraphs.EdgeWeightedSquareLattice(GRIDSIZE, [ μ for i in 1:2GRIDSIZE*(GRIDSIZE-1) ]);
Pgty = [ sum( aTraj[i,j].jointProb[g,:] ) for g in 1:size(aTraj[i,j].jointProb)[1] ] /
    aTraj[i,j].NgenSample

p = mGraphs.VertexWeightedSquareLattice( GRIDSIZE, Pgty )
pMat = mGraphs.matrixForm(p);

ax = figNum.add_subplot(3,3,j + 3*(3-i))

eMin = -18
for (i,e) in enumerate(pMat)
    if e == 0
       pMat[i] = 10.0^eMin
    end
end

imshow(log.(pMat),cmap="inferno",vmax=-1, vmin=eMin);
# title("β = $(β), μ = $(μ)");

xticks( prepend!(collect(3:4:GRIDSIZE-1),0), prepend!(collect(4:4:GRIDSIZE),1) )
yticks( prepend!(collect(3:4:GRIDSIZE-1),0), prepend!(collect(4:4:GRIDSIZE),1) )
colorbar(fraction=0.045)

savefig("robustness.pdf", bbox_inches="tight")
# -

# #### Genotype Distributions --- Theory

# +
using LinearAlgebra

# array of fitness values
fnaf(β) = [ exp( β * fTbl[g] ) for g in 1:DIMGSPACE ]

# construction of the sensitivity as function of the mutation probability
fnGrid(μ) = mGraphs.EdgeWeightedSquareLattice(GRIDSIZE, [ μ for i in 1:2GRIDSIZE*(GRIDSIZE-1) ])

mMut(μ) = Array(mGraphs.transitionMatrix(fnGrid(μ))) .+ Diagonal( ones(Float64, DIMGSPACE) )

mϕ(β,μ) = fnaf(β)' * mMut(μ)
Z(β,μ) = mϕ(β,μ) * fnaf(β)

prob(β,μ) = fnaf(β)' .* mϕ(β,μ) / Z(β,μ)

# +
using  Distances

Base.show(io::IO, f::Float64) = @printf(io, "%.2f", f)

figNum = figure(figsize=(12,12))

for (i,β) in enumerate(aSelStrength), (j,μ) in enumerate(aMutFactor)

    grid = mGraphs.EdgeWeightedSquareLattice(GRIDSIZE, [ μ for i in 1:2GRIDSIZE*(GRIDSIZE-1) ]);
    PgtyPop = [ sum( aTraj[i,j].jointProb[g,:] ) for g in 1:size(aTraj[i,j].jointProb)[1] ] / aTraj[i,j].NgenSample
    Pgty = prob(β,μ)

    p = mGraphs.VertexWeightedSquareLattice( GRIDSIZE, Pgty )
    pMat = mGraphs.matrixForm(p);

    ax = figNum.add_subplot(3,3,j + 3*(3-i))

    eMin = -18
#     for (i,e) in enumerate(pMat)
#         if e == 0
#            pMat[i] = 10.0^eMin
#         end
#     end

    imshow(log.(pMat),cmap="inferno");
#     title("β = $(β), μ = $(μ)");
    title("DKL = $(kl_divergence(PgtyPop, Pgty))");
    colorbar(fraction=0.045)

    xticks([])
    yticks([])

    if i == 1 xticks( prepend!(collect(3:4:GRIDSIZE-1),0), prepend!(collect(4:4:GRIDSIZE),1) ) end
    if j == 1 yticks( prepend!(collect(3:4:GRIDSIZE-1),0), prepend!(collect(4:4:GRIDSIZE),1) ) end
end

savefig("pmfPhaseDiagramTheo.pdf", bbox_inches="tight")
# -
# ## Sensitivity Analysis

# +
using LinearAlgebra

iβ = 2
β = aSelStrength[iβ]

# array of fitness values
af = [ exp( β * fTbl[g] ) for g in 1:DIMGSPACE ]
aθ = [ r == HF ? 1.0 : 0.0 for r in fTbl ]

mSensitivity = [ (af[go] - af[gt])/af[go] for gt in 1:DIMGSPACE, go in 1:DIMGSPACE ]

# construction of the sensitivity as function of the mutation probability
fnGrid(μ) = mGraphs.EdgeWeightedSquareLattice(GRIDSIZE, [ μ for i in 1:2GRIDSIZE*(GRIDSIZE-1) ])

mMut(μ) = Array(mGraphs.transitionMatrix(fnGrid(μ))) .+ Diagonal( ones(Float64, DIMGSPACE) )

mϕ(μ) = af' * mMut(μ)
Z(μ) = mϕ(μ) * af

prob(μ) = af' .* mϕ(μ) / Z(μ)
probNeutral(μ) = af' .* mϕ(μ) .* aθ' / sum(af' .* mϕ(μ) .* aθ')

Pref = af' / sum(af)
PrefNeutral = af' .* aθ' / sum(af' .* aθ')

Esensitivity(μ) = sum( mSensitivity .* mMut(μ) .* prob(μ) ) 
EsensitivityNeutral(μ) = sum( mSensitivity .* mMut(μ) .* probNeutral(μ) ) 

EsensitivityRef(μ) = sum( mSensitivity .* mMut(μ) .* Pref ) 
EsensitivityRefNeutral(μ) = sum( mSensitivity .* mMut(μ) .* PrefNeutral )
;

# +
expectedSensitivityPop = Vector{Float64}(undef, length(aMutFactor))
expectedSensitivityPopNeutral = Vector{Float64}(undef, length(aMutFactor))

for (iμ,μ) in enumerate(aMutFactor)
    Pgty = [ sum( aTraj[iβ,iμ].jointProb[g,:] ) for g in 1:size(aTraj[iβ,iμ].jointProb)[1] ] / aTraj[iβ,iμ].NgenSample
    PgtyNeutral = [ sum( aTraj[iβ,iμ].jointProb[g,:] ) * aθ[g] for g in 1:size(aTraj[iβ,iμ].jointProb)[1] ] /
    aTraj[iβ,iμ].NgenSample
    expectedSensitivityPop[iμ] = sum( mSensitivity .* mMut(μ) .* Pgty' )
    expectedSensitivityPopNeutral[iμ] = sum( mSensitivity .* mMut(μ) .* PgtyNeutral' )
end

# +
eμIterator = -3.0:0.05:-1.0

aμ = [ 2.5 * 10.0^eμ for eμ in eμIterator ]
as = [ Esensitivity(μ) for μ in aμ ] ;
asNeutral = [ EsensitivityNeutral(μ) for μ in aμ ] ;
asRef = [ EsensitivityRef(μ) for μ in aμ ] ;
asRefNeutral = [ EsensitivityRefNeutral(μ) for μ in aμ ] ;

# +
subplots(1,2, figsize=(10.0, 5.0))
ylims=[-0.01, 0.25]

subplot(121); title("Neutral Network");
    plot( aμ, asNeutral, color="tab:blue", label="genotypic dynamics" )
    plot( aμ, asRefNeutral, color="tab:blue", "--", label="reference dynamics" )
    scatter( aMutFactor, expectedSensitivityPopNeutral, color="tab:blue", marker="+", label="population dynamics" )

    xscale("log")
    # yscale("log")

    ylim(ylims)

    xlabel("mutation probability, μ")
    ylabel("expected sensitivity, s")

    legend()

subplot(122); title("Whole Genotypic Space");
    plot( aμ, as, label="genotypic dynamics" )
    plot( aμ, asRef, color="tab:blue", "--", label="reference dynamics" )
    scatter( aMutFactor, expectedSensitivityPop, color="tab:blue", marker="+", label="population dynamics" )

    xscale("log")
    # yscale("log")

    ylim(ylims)

    xlabel("mutation probability, μ")
#     ylabel("expected sensitivity, s")

    legend(loc="upper left")

tight_layout();

savefig("sensitivity.pdf", bbox_inches="tight")
# -
function squareGridPlot(L; l=1, fsize=7, s0=30, colormap="cividis")
    
    V = [ [ (i-1)÷L + 1, (i-1)%L + 1 ] for i in 1:L^2 ]

    fig = figure(figsize=(fsize,fsize))
    for (i,v) in enumerate(V)
        scatter([v[1]], [v[2]], s=s0, cmap=colormap, c=10, alpha=0.5, marker="s")
    end

    xticks( prepend!(collect(l:l:L),1) )
    yticks( prepend!(collect(l:l:L),1), append!(collect(L:-l:l),1) )
end

# +
squareGridPlot(GRIDSIZE, l=4, fsize=4) 

savefig("genotypicGrid.pdf", bbox_inches="tight")
;
