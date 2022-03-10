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

folderName = "/Users/riccardorao/projects/fisheria/analysis/gridWorld/scratch/sty/";

# +
# null ID for checking parameters
# jobID = "styGridEvo"

# phase diagram results
jobID = "styGridEvo-7785591"
include(folderName * jobID * "-parameters.jl")

# +
# sensitivity analysis results---several population sizes
# aJobID = [ "styGridEvo-7786606", "styGridEvo-7786979", "styGridEvo-7786873", "styGridEvo-7785776" ]

# sensitivity analysis results---N=100 population spanning more mutation rate values
aJobID = [ "styGridEvo-7787646", "styGridEvo-7785776" ] # "styGridEvo-7787631", 

include(folderName * aJobID[end] * "-parameters.jl")

# +
aaMutFactor = Vector{Vector{Float64}}(undef, length(aJobID))
aNpop = Vector{Int32}(undef, length(aJobID))

for (iID, jobID) in enumerate(aJobID)
    include(folderName * jobID * "-parameters.jl")
    aaMutFactor[iID] = aMutFactor
end
;

# +
fMat = mGraphs.matrixForm( mGraphs.VertexWeightedSquareLattice( GRIDSIZE, fTbl ) )

# imshow(fMat,cmap="cividis"); # origin="lower"
imshow(fMat, cmap="copper"); # origin="lower"

xticks( prepend!(collect(3:4:GRIDSIZE-1),0), prepend!(collect(4:4:GRIDSIZE),1) )
yticks( prepend!(collect(3:4:GRIDSIZE-1),0), prepend!(collect(4:4:GRIDSIZE),1) )
colorbar()

# savefig("styLnF.pdf", bbox_inches="tight")
# -

# # Data Generation/Collection

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
# phase diagram results
using JLD, HDF5, DelimitedFiles

aTraj = mUtils.readJLD(folderName * jobID * "_aTraj.jld", "aTraj");

# +
# sensitivity analysis results
using JLD, HDF5, DelimitedFiles

aTraj = Vector{Array{EvolutionaryDynamics.TabularSystems.TrajectoryData,2}}(undef, length(aJobID))
for (i, jID) in enumerate(aJobID)
    aTraj[i] = mUtils.readJLD(folderName * jID * "_aTraj.jld", "aTraj")
end
# -

# # Analysis
#
# ### Checks about relaxation

# +
iβ = 2
iμ = 2

subplots(2,1,figsize=(7.5,7.5))
subplot(211); title("Average Fitness");
    plot( aTraj[iβ,iμ].avePerformance );
    title("β = $(aSelStrength[iβ]), μ = $(aMutFactor[iμ])");
subplot(212); title("Mutation Factor");
    plot( aTraj[iβ,iμ].mutationFactor );
tight_layout();
# -

# ### Genotype Distributions

# +
Base.show(io::IO, f::Float64) = @printf(io, "%.4f", f)

figNum = figure(figsize=(6.4,6.4))

for (i,β) in enumerate(aSelStrength), (j,μ) in enumerate(aMutFactor[1:3])

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
    
    imshow(log10.(pMat), cmap="cividis", vmax=-1, vmin=eMin);
#     title("β = $(β), μ = $(μ)");

    xticks([])
    yticks([])

    if i == 1 xticks( prepend!(collect(3:4:GRIDSIZE-1),0), prepend!(collect(4:4:GRIDSIZE),1) ) end
    if j == 1 yticks( prepend!(collect(3:4:GRIDSIZE-1),0), prepend!(collect(4:4:GRIDSIZE),1) ) end
    if i == 1 && j == 3 colorbar(fraction=0.045) end
end

savefig("styPmfPhaseDiagramNumRaw.pdf", bbox_inches="tight")
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

# ### Theoretical Genotype Distributions

# +
using LinearAlgebra

iβ = 2
β = aSelStrength[iβ]

# array of fitness values
fnaf(β) = [ exp( β * fTbl[g] ) for g in 1:DIMGSPACE ]

# construction of the sensitivity as function of the mutation probability
fnGrid(μ) = mGraphs.EdgeWeightedSquareLattice(GRIDSIZE, [ μ for i in 1:2GRIDSIZE*(GRIDSIZE-1) ])

mMut(μ) = Array(mGraphs.transitionMatrix(fnGrid(μ))) .+ Diagonal( ones(Float64, DIMGSPACE) )
mAdj = Array(mGraphs.adjacencyMatrix(fnGrid(0.1)))


# probability vectors of n-th order genotypic dynamics
mϕ(β, μ, n=1) = sum( ( fnaf(β) .* mMut(μ) )^n, dims=1 )
Z(β, μ, n=1) = mϕ(β, μ, n) * fnaf(β)

prob(β, μ, n=1) = fnaf(β)' .* mϕ(β, μ, n) ./ Z(β, μ, n)

function probInf(β, μ)
    r = real( eigvecs( fnaf(β) .* mMut(μ) )[:,end] )
    map!( x -> x >= 0 ? x : 0.0, r, r )
    return r' ./ sum(r)
end

# +
using  Distances

Base.show(io::IO, f::Float64) = @printf(io, "%.2f", f)

nExpansion = 1

figNum = figure(figsize=(10.7, 10.7))

for (i,β) in enumerate(aSelStrength), (j,μ) in enumerate(aMutFactor)

    grid = mGraphs.EdgeWeightedSquareLattice(GRIDSIZE, [ μ for i in 1:2GRIDSIZE*(GRIDSIZE-1) ]);
    PgtyPop = [ sum( aTraj[i,j].jointProb[g,:] ) for g in 1:size(aTraj[i,j].jointProb)[1] ] / aTraj[i,j].NgenSample
    Pgty = prob(β, μ, nExpansion)
#     Pgty = probInf(β, μ)

    p = mGraphs.VertexWeightedSquareLattice( GRIDSIZE, Pgty )
    pMat = mGraphs.matrixForm(p);

    ax = figNum.add_subplot(3,3,j + 3*(3-i))

    eMin = -18
#     for (i,e) in enumerate(pMat)
#         if e == 0
#            pMat[i] = 10.0^eMin
#         end
#     end

    imshow(log10.(pMat), cmap="cividis");
#     title("β = $(β), μ = $(μ)");
    title("DKL = $(kl_divergence(PgtyPop, Pgty))");
    colorbar(fraction=0.045)

    xticks([])
    yticks([])

    if i == 1 xticks( prepend!(collect(3:4:GRIDSIZE-1),0), prepend!(collect(4:4:GRIDSIZE),1) ) end
    if j == 1 yticks( prepend!(collect(3:4:GRIDSIZE-1),0), prepend!(collect(4:4:GRIDSIZE),1) ) end
end

savefig("styPmfPhaseDiagramTheo$(nExpansion)Raw.pdf", bbox_inches="tight")
# -
# # Sensitivity Analysis

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
mAdj = Array(mGraphs.adjacencyMatrix(fnGrid(0.1)))


# probability vectors of n-th order genotypic dynamics
mϕ(μ, n=1) = sum( ( af .* mMut(μ) )^n, dims=1 )
Z(μ, n=1) = mϕ(μ, n) * af

prob(μ, n=1) = af' .* mϕ(μ, n) ./ Z(μ, n)
probNeutral(μ, n=1) = af' .* mϕ(μ, n) .* aθ' / sum(af' .* mϕ(μ, n) .* aθ')

function probInf(μ)
    r = real( eigvecs( af .* mMut(μ) )[:,end] )
    return r' ./ sum(r)
end

function probInfNeutral(μ)
    r = real( eigvecs( af' .* mMut(μ) )[:,end] )
    return ( r .* aθ )' ./ sum(r .* aθ)
end

# expected sensitivity evaluations
# Esensitivity(μ, n=1) = sum( mSensitivity .* mMut(μ) .* prob(μ, n) )
# EsensitivityNeutral(μ, n=1) = sum( mSensitivity .* mMut(μ) .* probNeutral(μ, n) )

# EsensitivityInf(μ) = sum( mSensitivity .* mMut(μ) .* probInf(μ) )
# EsensitivityInfNeutral(μ) = sum( mSensitivity .* mMut(μ) .* probInfNeutral(μ) )

# expected sensitivity evaluations
Esensitivity(μ, n=1) = sum( mSensitivity .* mMut(μ) .* prob(μ, n) ) / μ
EsensitivityNeutral(μ, n=1) = sum( mSensitivity .* mMut(μ) .* probNeutral(μ, n) ) / μ

EsensitivityInf(μ) = sum( mSensitivity .* mMut(μ) .* probInf(μ) ) / μ
EsensitivityInfNeutral(μ) = sum( mSensitivity .* mMut(μ) .* probInfNeutral(μ) ) / μ

# average fitness values
aveF(μ, n=1) = prob(μ, n) ⋅ af
aveFInf(μ) = probInf(μ) ⋅ af
;

# +
aExpectedSensitivityPop = Vector{Vector{Float64}}(undef, length(aJobID))
aExpectedSensitivityPopNeutral = Vector{Vector{Float64}}(undef, length(aJobID))

aVarianceSensitivityPop = Vector{Vector{Float64}}(undef, length(aJobID))
aVarianceSensitivityPopNeutral = Vector{Vector{Float64}}(undef, length(aJobID))

for iID in eachindex(aJobID)
    aExpectedSensitivityPop[iID] = Vector{Float64}(undef, length(aaMutFactor[iID]))
    aExpectedSensitivityPopNeutral[iID] = Vector{Float64}(undef, length(aaMutFactor[iID]))
    
    aVarianceSensitivityPop[iID] = Vector{Float64}(undef, length(aaMutFactor[iID]))
    aVarianceSensitivityPopNeutral[iID] = Vector{Float64}(undef, length(aaMutFactor[iID]))
    
    for (iμ ,μ) in enumerate(aaMutFactor[iID])
        Pgty = [ sum( aTraj[iID][iβ,iμ].jointProb[g,:] ) for
            g in 1:size(aTraj[iID][iβ,iμ].jointProb)[1] ] / aTraj[iID][iβ,iμ].NgenSample
        PgtyNeutral = [ sum( aTraj[iID][iβ,iμ].jointProb[g,:] ) * aθ[g] for
            g in 1:size(aTraj[iID][iβ,iμ].jointProb)[1] ] / aTraj[iID][iβ,iμ].NgenSample
        
#         aExpectedSensitivityPop[iID][iμ] = sum( mSensitivity .* mMut(μ) .* Pgty' )
#         aExpectedSensitivityPopNeutral[iID][iμ] = sum( mSensitivity .* mMut(μ) .* PgtyNeutral' )
        
        aExpectedSensitivityPop[iID][iμ] = sum( ( mSensitivity .* mMut(μ) ) * Pgty ) / μ
        aExpectedSensitivityPopNeutral[iID][iμ] = sum( ( mSensitivity .* mMut(μ) ) * PgtyNeutral ) / μ
        
        aVarianceSensitivityPop[iID][iμ] = ( sum( mSensitivity .* mMut(μ), dims=1 ).^2 * Pgty / μ^2 )[1] -
            aExpectedSensitivityPop[iID][iμ]^2
        aVarianceSensitivityPopNeutral[iID][iμ] = ( sum( mSensitivity .* mMut(μ), dims=1 ).^2 * PgtyNeutral /
            μ^2 )[1] - aExpectedSensitivityPopNeutral[iID][iμ]^2
    end
end

# +
# sensitivity vectors of theoretical models

an = [0, 1, 100]
eμIterator = -5.0:0.05:-1.0
aμ = [ 2.5 * 10.0^eμ for eμ in eμIterator ]

aSens = Vector{Vector{Float64}}(undef, length(an))
aSensNeutral = Vector{Vector{Float64}}(undef, length(an))
aAveF = Vector{Vector{Float64}}(undef, length(an))

for (i, n) in enumerate(an)
    aSens[i] = [ Esensitivity(μ, n) for μ in aμ ]
    aSensNeutral[i] = [ EsensitivityNeutral(μ, n) for μ in aμ ]
    aAveF[i] = [ aveF(μ, n) for μ in aμ ]
end

aSensInf = [ EsensitivityInf(μ) for μ in aμ ]
aSensInfNeutral = [ EsensitivityInfNeutral(μ) for μ in aμ ]
aAveFInf = [ aveFInf(μ) for μ in aμ ]
;

# +
subplots(1,2, figsize=(10.0, 5.0))
# ylims = [-0.01, 0.25]

subplot(121); title("Neutral Network");
    for (i, n) in enumerate(an)
        plot( aμ, aSensNeutral[i], label="$(n)-th order genotypic dynamics" )
    end 
    plot( aμ, aSensInfNeutral, color="tab:red", "--", label="∞-th order genotypic dynamics" )
    for (iID, expectedSensitivityPopNeutral) in enumerate(aExpectedSensitivityPopNeutral)
        errorbar( aaMutFactor[iID], expectedSensitivityPopNeutral, sqrt.(aVarianceSensitivityPopNeutral[iID]),
            color="tab:green", fmt=".", label="population dynamics" )
    end

    xscale("log")
#     yscale("log")
#     ylim(ylims)

    xlabel("mutation probability, μ")
    ylabel("expected sensitivity, s")

#     legend()

subplot(122); title("Whole Genotypic Space");
    for (i, n) in enumerate(an)
        plot( aμ, aSens[i], label="$(n)-th order genotypic dynamics" )
    end 
    plot( aμ, aSensInf, color="tab:red", "--", label="∞-th order genotypic dynamics" )
    for (iID, expectedSensitivityPop) in enumerate(aExpectedSensitivityPop)
        errorbar( aaMutFactor[iID], expectedSensitivityPop, sqrt.(aVarianceSensitivityPop[iID]),
        color="tab:green", fmt=".", label="population dynamics" )
    end

    xscale("log")
#     yscale("log")
#     ylim(ylims)

    xlabel("mutation probability, μ")
#     ylabel("expected sensitivity, s")

#     legend(loc="upper left")

tight_layout();

savefig("sensitivity.pdf", bbox_inches="tight")
# +
figure()
ylims = [0.2, 2.0]

plot( aμ, aSensNeutral[1], label="homogeneous neutral space" )
plot( aμ, aSensNeutral[2], label="single genotype dynamics" )
plot( aμ, aSensNeutral[3], label="$(an[3]) genotypes dynamics", "--" )
plot( aμ, aSensInfNeutral, label="∞ genotypes dynamics", "--", )
scatter( aaMutFactor[1], aExpectedSensitivityPopNeutral[1], color="tab:green", marker="*" )
scatter( aaMutFactor[2], label="population dynamics, N = 100", aExpectedSensitivityPopNeutral[2],
    color="tab:green", marker="*" )

xscale("log")
ylim(ylims)

xlabel("mutation probability, μ")
ylabel("rescaled average sensitivity, ⟨σ⟩/μ")

legend()

savefig("sensitivity.pdf", bbox_inches="tight")

# +
title("Whole Genotypic Space");
for (i, n) in enumerate(an)
    plot( aμ, aAveF[i], label="$(n)-th order genotypic dynamics" )
end 
plot( aμ, aAveFInf, color="tab:red", "--", label="∞-th order genotypic dynamics" )

xscale("log")

xlabel("mutation probability, μ")
ylabel("average Fitness, ⟨f⟩")

tight_layout();

# savefig("sensitivity.pdf", bbox_inches="tight")
# -
# ### Strength of Selection dependence

# +
using LinearAlgebra

# array of fitness values
fnf(β) = [ exp( β * fTbl[g] ) for g in 1:DIMGSPACE ]
aθ = [ r == HF ? 1.0 : 0.0 for r in fTbl ]

mSensitivityFn(β) = [ (fnf(β)[go] - fnf(β)[gt])/fnf(β)[go] for gt in 1:DIMGSPACE, go in 1:DIMGSPACE ]
mLogSensitivityFn(β) = [ log(fnf(β)[go] / fnf(β)[gt]) for gt in 1:DIMGSPACE, go in 1:DIMGSPACE ]

# construction of the sensitivity as function of the mutation probability
fnGrid(μ) = mGraphs.EdgeWeightedSquareLattice(GRIDSIZE, [ μ for i in 1:2GRIDSIZE*(GRIDSIZE-1) ])

mMut(μ) = Array(mGraphs.transitionMatrix(fnGrid(μ))) .+ Diagonal( ones(Float64, DIMGSPACE) )
mAdj = Array(mGraphs.adjacencyMatrix(fnGrid(0.1))) ;

# probability vectors of n-th order genotypic dynamics
mϕ(β, μ, n=1) = sum( ( fnf(β) .* mMut(μ) )^n, dims=1 )
Z(β, μ, n=1) = mϕ(β, μ, n) * fnf(β)

prob(β, μ, n=1) = fnf(β)' .* mϕ(β, μ, n) ./ Z(β, μ, n)
probNeutral(β, μ, n=1) = fnf(β)' .* mϕ(β, μ, n) .* aθ' / sum(fnf(β)' .* mϕ(β, μ, n) .* aθ')

function probInf(β, μ)
    r = real( eigvecs( fnf(β) .* mMut(μ) )[:,end] )
    return r' ./ sum(r)
end

function probInfNeutral(β, μ)
    r = real( eigvecs( fnf(β)' .* mMut(μ) )[:,end] )
    return ( r .* aθ )' ./ sum(r .* aθ)
end

# expected sensitivity evaluations
Esensitivity(β, μ, n=1) = sum( mSensitivityFn(β) .* mMut(μ) .* prob(β, μ, n) ) / μ
EsensitivityLog(β, μ, n=1) = sum( mLogSensitivityFn(β) .* mMut(μ) .* prob(β, μ, n) ) / μ
EsensitivityNeutral(β, μ, n=1) = sum( mSensitivityFn(β) .* mMut(μ) .* probNeutral(β, μ, n) ) / μ

EsensitivityInf(β, μ) = sum( mSensitivityFn(β) .* mMut(μ) .* probInf(β, μ) ) / μ
EsensitivityInfLog(β, μ) = sum( mLogSensitivityFn(β) .* mMut(μ) .* probInf(β, μ) ) / μ
EsensitivityInfNeutral(β, μ) = sum( mSensitivityFn(β) .* mMut(μ) .* probInfNeutral(β, μ) ) / μ


Vsensitivity(β, μ, Eσ, n=1) = ( sum( mSensitivityFn(β) .* mMut(μ), dims=1 ).^2 * prob(β, μ, n)' / μ^2 )[1] - Eσ^2 
VsensitivityLog(β, μ, Eσ, n=1) = ( sum( mLogSensitivityFn(β) .* mMut(μ), dims=1 ).^2 * prob(β, μ, n)' / μ^2 )[1] - Eσ^2
VsensitivityNeutral(β, μ, Eσ, n=1) = ( sum( mSensitivityFn(β) .* mMut(μ), dims=1 ).^2 * probNeutral(β, μ, n)'
    / μ^2 )[1] - Eσ^2

# +
# sensitivity vectors of theoretical models

aβ = [0.1, 1.0]
n = 1

eμIterator = -5.0:0.05:-1.0
aμ = [ 2.5 * 10.0^eμ for eμ in eμIterator ]

aSens = Vector{Vector{Float64}}(undef, length(aβ))
aSensLog = Vector{Vector{Float64}}(undef, length(aβ))
aSensNeutral = Vector{Vector{Float64}}(undef, length(aβ))

aVarSens = Vector{Vector{Float64}}(undef, length(aβ))
aVarSensLog = Vector{Vector{Float64}}(undef, length(aβ))
aVarSensNeutral = Vector{Vector{Float64}}(undef, length(aβ))

for (i, β) in enumerate(aβ)
    aSens[i] = [ Esensitivity(β, μ, n) for μ in aμ ]
    aSensLog[i] = [ EsensitivityLog(β, μ, n) for μ in aμ ]
    aSensNeutral[i] = [ EsensitivityNeutral(β, μ, n) for μ in aμ ]
    
    aVarSens[i] = [ Vsensitivity(β, μ, aSens[i][iμ], n) for (iμ, μ) in enumerate(aμ) ]
    aVarSensLog[i] = [ VsensitivityLog(β, μ, aSensLog[i][iμ], n) for (iμ, μ) in enumerate(aμ) ]
    aVarSensNeutral[i] = [ VsensitivityNeutral(β, μ, aSensNeutral[i][iμ], n) for (iμ, μ) in enumerate(aμ) ]
end ;

# +
subplots(2, 3, figsize=(12.0, 10.0))
# ylims = [-0.01, 0.25]

subplot(231); title("Neutral Network");
    for (i, β) in enumerate(aβ)
        plot( aμ, aSensNeutral[i] )
        fill_between( aμ, aSensNeutral[i] .- sqrt.(aVarSensNeutral[i]),
            aSensNeutral[i] .+ sqrt.(aVarSensNeutral[i]), alpha=0.2 )
    end 
    xscale("log")
    xlabel("mutation probability, μ")
    ylabel("expected sensitivity, σ")

subplot(232); title("Log Sensitivity");
    for (i, β) in enumerate(aβ)
        plot( aμ, aSensLog[i] )
        fill_between( aμ, aSensLog[i] .- sqrt.(aVarSensLog[i]),
            aSensLog[i] .+ sqrt.(aVarSensLog[i]), alpha=0.2 )
    end 
    xscale("log")
    xlabel("mutation probability, μ")

subplot(233); title("Whole Genotypic Space");
    for (i, β) in enumerate(aβ)
        plot( aμ, aSens[i] )
        fill_between( aμ, aSens[i] .- sqrt.(aVarSens[i]),
            aSens[i] .+ sqrt.(aVarSens[i]), alpha=0.2 )
    end 
    xscale("log")
    xlabel("mutation probability, μ")

subplot(234); title("Neutral Network");
    for (i, β) in enumerate(aβ)
        plot( aμ, aVarSensNeutral[i] )
    end 
    xscale("log")
    xlabel("mutation probability, μ")
    ylabel("sensitivity variance, var{σ}")

subplot(235); title("Log Sensitivity");
    for (i, β) in enumerate(aβ)
        plot( aμ, aVarSensLog[i] )
    end 
    xscale("log")
    xlabel("mutation probability, μ")

subplot(236); title("Whole Genotypic Space");
    for (i, β) in enumerate(aβ)
        plot( aμ, aVarSens[i] )
    end 
    xscale("log")
    xlabel("mutation probability, μ")

tight_layout();

# savefig("sensitivity.pdf", bbox_inches="tight")
# +
subplots(2, 3, figsize=(12.0, 10.0))
# ylims = [-0.01, 0.25]

subplot(231); title("Neutral Network");
    for (i, β) in enumerate(aβ)
        plot( aμ, aSensNeutral[i] )
        fill_between( aμ, aSensNeutral[i] .- sqrt.(aVarSensNeutral[i]),
            aSensNeutral[i] .+ sqrt.(aVarSensNeutral[i]), alpha=0.2 )
    end 
    xscale("log")
    xlabel("mutation probability, μ")
    ylabel("expected sensitivity, σ")

subplot(232); title("Log Sensitivity");
    for (i, β) in enumerate(aβ)
        plot( aμ, aSensLog[i] )
        fill_between( aμ, aSensLog[i] .- sqrt.(aVarSensLog[i]),
            aSensLog[i] .+ sqrt.(aVarSensLog[i]), alpha=0.2 )
    end 
    xscale("log")
    xlabel("mutation probability, μ")

subplot(233); title("Whole Genotypic Space");
    for (i, β) in enumerate(aβ)
        plot( aμ, aSens[i] )
        fill_between( aμ, aSens[i] .- sqrt.(aVarSens[i]),
            aSens[i] .+ sqrt.(aVarSens[i]), alpha=0.2 )
    end 
    xscale("log")
    xlabel("mutation probability, μ")

subplot(234); title("Neutral Network");
    for (i, β) in enumerate(aβ)
        plot( aμ, aVarSensNeutral[i] )
    end 
    xscale("log")
    xlabel("mutation probability, μ")
    ylabel("sensitivity variance, var{σ}")

subplot(235); title("Log Sensitivity");
    for (i, β) in enumerate(aβ)
        plot( aμ, aVarSensLog[i] )
    end 
    xscale("log")
    xlabel("mutation probability, μ")

subplot(236); title("Whole Genotypic Space");
    for (i, β) in enumerate(aβ)
        plot( aμ, aVarSens[i] )
    end 
    xscale("log")
    xlabel("mutation probability, μ")

tight_layout();

# savefig("sensitivity.pdf", bbox_inches="tight")
# -
# ### Fluctuations--Sensitivity Theorem

# +
using LinearAlgebra

iβ = 2
β = aSelStrength[iβ]

# array of fitness values
af = [ exp( β * fTbl[g] ) for g in 1:DIMGSPACE ]
aθ = [ r == HF ? 1.0 : 0.0 for r in fTbl ]

mSensitivity = [ (af[go] - af[gt])/af[go] for gt in 1:DIMGSPACE, go in 1:DIMGSPACE ]
mLogSensitivity = [ log(af[go] / af[gt]) for gt in 1:DIMGSPACE, go in 1:DIMGSPACE ]

# construction of the sensitivity as function of the mutation probability
fnGrid(μ) = mGraphs.EdgeWeightedSquareLattice(GRIDSIZE, [ μ for i in 1:2GRIDSIZE*(GRIDSIZE-1) ])

mMut(μ) = Array(mGraphs.transitionMatrix(fnGrid(μ))) .+ Diagonal( ones(Float64, DIMGSPACE) )
mAdj = Array(mGraphs.adjacencyMatrix(fnGrid(0.1))) ;

# +
aExpectedSensitivityPop        = Vector{Vector{Float64}}(undef, length(aJobID))
aExpectedSensitivityPopNeutral = Vector{Vector{Float64}}(undef, length(aJobID))

aVarianceSensitivityPop        = Vector{Vector{Float64}}(undef, length(aJobID))
aVarianceSensitivityPopNeutral = Vector{Vector{Float64}}(undef, length(aJobID))

aChangeSensitivityPop        = Vector{Vector{Float64}}(undef, length(aJobID))
aChangeSensitivityPopNeutral = Vector{Vector{Float64}}(undef, length(aJobID))

aDivergenceSensitivityPop        = Vector{Vector{Float64}}(undef, length(aJobID))
aDivergenceSensitivityPopNeutral = Vector{Vector{Float64}}(undef, length(aJobID))

for iID in eachindex(aJobID)
    aExpectedSensitivityPop[iID]        = Vector{Float64}(undef, length(aaMutFactor[iID]))
    aExpectedSensitivityPopNeutral[iID] = Vector{Float64}(undef, length(aaMutFactor[iID]))
    
    aVarianceSensitivityPop[iID]        = Vector{Float64}(undef, length(aaMutFactor[iID]))
    aVarianceSensitivityPopNeutral[iID] = Vector{Float64}(undef, length(aaMutFactor[iID]))
    
    aChangeSensitivityPop[iID]        = Vector{Float64}(undef, length(aaMutFactor[iID])-1)
    aChangeSensitivityPopNeutral[iID] = Vector{Float64}(undef, length(aaMutFactor[iID])-1)
    
    aDivergenceSensitivityPop[iID]        = Vector{Float64}(undef, length(aaMutFactor[iID])-1)
    aDivergenceSensitivityPopNeutral[iID] = Vector{Float64}(undef, length(aaMutFactor[iID])-1)
    
    for (iμ ,μ) in enumerate(aaMutFactor[iID])
        Pgty = [ sum( aTraj[iID][iβ,iμ].jointProb[g,:] ) for
            g in 1:size(aTraj[iID][iβ,iμ].jointProb)[1] ] / aTraj[iID][iβ,iμ].NgenSample
        PgtyNeutral = [ sum( aTraj[iID][iβ,iμ].jointProb[g,:] ) * aθ[g] for
            g in 1:size(aTraj[iID][iβ,iμ].jointProb)[1] ]
        PgtyNeutral /= sum(PgtyNeutral)
        
#         aExpectedSensitivityPop[iID][iμ] = sum( mSensitivity .* mMut(μ) .* Pgty' )
#         aExpectedSensitivityPopNeutral[iID][iμ] = sum( mSensitivity .* mMut(μ) .* PgtyNeutral' )
        
        aExpectedSensitivityPop[iID][iμ]        = sum( ( mLogSensitivity .* mMut(μ) ) * Pgty ) / μ
        aExpectedSensitivityPopNeutral[iID][iμ] = sum( ( mSensitivity .* mMut(μ) ) * PgtyNeutral ) / μ
        
        aVarianceSensitivityPop[iID][iμ]        = ( sum( mLogSensitivity .* mMut(μ), dims=1 ).^2 * Pgty /
            μ^2 )[1] - aExpectedSensitivityPop[iID][iμ]^2
        aVarianceSensitivityPopNeutral[iID][iμ] = ( sum( mSensitivity .* mMut(μ), dims=1 ).^2 * PgtyNeutral /
            μ^2 )[1] - aExpectedSensitivityPopNeutral[iID][iμ]^2
    end
        
    for (iμ, μ) in enumerate(aaMutFactor[iID][1:end-1])
        aChangeSensitivityPop[iID][iμ] =
            ( aExpectedSensitivityPop[iID][iμ+1] - aExpectedSensitivityPop[iID][iμ] ) /
            ( - aaMutFactor[iID][iμ+1] + μ )
        aChangeSensitivityPopNeutral[iID][iμ] =
            ( aExpectedSensitivityPopNeutral[iID][iμ+1] - aExpectedSensitivityPopNeutral[iID][iμ] ) /
            ( - aaMutFactor[iID][iμ+1] + μ )
        
        aDivergenceSensitivityPop[iID][iμ] =
            abs( aChangeSensitivityPop[iID][iμ] / aVarianceSensitivityPop[iID][iμ] )
        aDivergenceSensitivityPopNeutral[iID][iμ] =
            abs( aChangeSensitivityPopNeutral[iID][iμ] / aVarianceSensitivityPopNeutral[iID][iμ] )
    end
end

# +
subplots(3, 2, figsize=(10.0, 10.0))

subplot(321); title("Neutral Space: Expected Sensitivity");
    for (iID, expectedSensitivityPopNeutral) in enumerate(aExpectedSensitivityPopNeutral)
        errorbar( aaMutFactor[iID], expectedSensitivityPopNeutral, aVarianceSensitivityPopNeutral[iID]/2,
            color="tab:blue", fmt="." )
    end
    axhline(y=0, linestyle="--")
    xscale("log")
    ylabel("expected sensitivity, ⟨σ⟩")

subplot(322); title("Whole Genotypic Space: Expected Sensitivity");
    for (iID, expectedSensitivityPop) in enumerate(aExpectedSensitivityPop)
        errorbar( aaMutFactor[iID], expectedSensitivityPop, aVarianceSensitivityPop[iID]/2,
            color="tab:blue", fmt="." )
    end
    axhline(y=0, linestyle="--")
    xscale("log")

subplot(323); title("Neutral Space: Sensitivity Variance");
    for (iID, c) in enumerate(aVarianceSensitivityPopNeutral)
        errorbar( aaMutFactor[iID], c, color="tab:blue", fmt="*" )
    end
    xscale("log")
    ylabel("sensitivity variance, ⟨σ²⟩ - ⟨σ⟩²")

subplot(324); title("Whole Genotypic Space: Sensitivity Variance");
    for (iID, c) in enumerate(aVarianceSensitivityPop)
        errorbar( aaMutFactor[iID], c, color="tab:blue", fmt="*" )
    end
    xscale("log")

ylims = [0.01, 10^4]
subplot(325); title("Neutral Network: Deviations from DB");
    for (iID, c) in enumerate(aDivergenceSensitivityPopNeutral)
#     for (iID, c) in enumerate(aChangeSensitivityPopNeutral)
        errorbar( aaMutFactor[iID][1:end-1], c, color="tab:blue", fmt="*", label="population dynamics" )
    end
    axhline(y=1, linestyle="--")
    xscale("log")
    xlabel("mutation probability, μ")
    yscale("log")
    ylabel("log absolute relative deviation")
    ylim(ylims)

subplot(326); title("Whole Genotypic Space: Deviations from DB");
    for (iID, c) in enumerate(aDivergenceSensitivityPop)
#     for (iID, c) in enumerate(aChangeSensitivityPop)
        errorbar( aaMutFactor[iID][1:end-1], c, color="tab:blue", fmt="*", label="population dynamics" )
    end
    axhline(y=1, linestyle="--")
    xscale("log")
    xlabel("mutation probability, μ")
    yscale("log")
    ylim(ylims)

tight_layout();

savefig("sensitivity.pdf", bbox_inches="tight")
# -
# # Information Theoretic Analysis

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
mAdj = Array(mGraphs.adjacencyMatrix(fnGrid(0.1)))


# probability vectors of n-th order genotypic dynamics
mϕ(μ, n=1) = sum( ( af .* mMut(μ) )^n, dims=1 )
Z(μ, n=1) = mϕ(μ, n) * af

prob(μ, n=1) = af' .* mϕ(μ, n) ./ Z(μ, n)
probNeutral(μ, n=1) = af' .* mϕ(μ, n) .* aθ' / sum(af' .* mϕ(μ, n) .* aθ')

function probInf(μ)
    r = real( eigvecs( af .* mMut(μ) )[:,end] )
    return abs.( r' ./ sum(r) )
end

function probInfNeutral(μ)
    r = real( eigvecs( af' .* mMut(μ) )[:,end] )
    return abs.( ( r .* aθ )' ./ sum(r .* aθ) )
end

# expected sensitivity evaluations
entropy(μ, n=1) = - sum( log.(prob(μ, n)) .* prob(μ, n) )
entropyNeutral(μ, n=1) = - sum( mUtils.ReLog.(probNeutral(μ, n)) .* probNeutral(μ, n) )

entropyInf(μ) = - sum( log.(probInf(μ)) .* probInf(μ) )
entropyInfNeutral(μ) = - sum( mUtils.ReLog.(probInfNeutral(μ)) .* probInfNeutral(μ) )

# +
# sensitivity vectors of theoretical models

an = [0, 1, 100]
eμIterator = -5.0:0.05:-1.0
aμ = [ 2.5 * 10.0^eμ for eμ in eμIterator ]

aEntropy = Vector{Vector{Float64}}(undef, length(an))
aEntropyNeutral = Vector{Vector{Float64}}(undef, length(an))

for (i, n) in enumerate(an)
    aEntropy[i] = [ entropy(μ, n) for μ in aμ ]
    aEntropyNeutral[i] = [ entropyNeutral(μ, n) for μ in aμ ]
end

aEntropyInf = [ entropyInf(μ) for μ in aμ ]
aEntropyInfNeutral = [ entropyInfNeutral(μ) for μ in aμ ]
;

# +
subplots(1,2, figsize=(10.0, 5.0))
# ylims = [-0.01, 0.25]

subplot(121); title("Neutral Network");
    for (i, n) in enumerate(an)
        plot( aμ, aEntropyNeutral[i], label="$(n)-th order genotypic dynamics" )
    end 
#     plot( aμ, aEntropyInfNeutral, color="tab:red", "--", label="∞-th order genotypic dynamics" )
#     for (iID, expectedEntropyitivityPopNeutral) in enumerate(aExpectedEntropyitivityPopNeutral)
#         scatter( aaMutFactor[iID], expectedEntropyitivityPopNeutral, color="tab:green", marker="*",
#             label="population dynamics" )
#     end

    xscale("log")
#     yscale("log")
#     ylim(ylims)

    xlabel("mutation probability, μ")
    ylabel("Shannon Entropy, H")

#     legend()

subplot(122); title("Whole Genotypic Space");
    for (i, n) in enumerate(an)
        plot( aμ, aEntropy[i], label="$(n)-th order genotypic dynamics" )
    end 
    plot( aμ, aEntropyInf, color="tab:red", "--", label="∞-th order genotypic dynamics" )
#     for (iID, expectedEntropyitivityPop) in enumerate(aExpectedEntropyitivityPop)
#         scatter( aaMutFactor[iID], expectedEntropyitivityPop, color="tab:green", marker="*",
#             label="population dynamics" )
#     end

    xscale("log")
#     yscale("log")
#     ylim(ylims)

    xlabel("mutation probability, μ")
#     ylabel("expected sensitivity, s")

#     legend(loc="upper left")

tight_layout();

savefig("sensitivity.pdf", bbox_inches="tight")
# -
# # Miscellanea

# +
# fit2clr = Dict( 1.0 => "tab:blue", 5.0 => "tab:red" )
fit2clr = Dict( 1.0 => "tab:grey", 5.0 => "tab:orange" )
fit2shp = Dict( 1.0 => "s", 5.0 => "s" )
fit2siz = Dict( 1.0 => 100, 5.0 => 200 )

function squareGridPlot(L; l=1, fsize=7, s0=210, colormap="cividis")
    
    V = [ [ (i-1)÷L + 1, (i-1)%L + 1 ] for i in 1:L^2 ]

    fig = figure(figsize=(fsize,fsize))
    for (i,v) in enumerate(V)
        scatter([v[1]], [v[2]], s=fit2siz[fTbl[i]], cmap=colormap,
            c=fit2clr[fTbl[i]], alpha=1.0, marker=fit2shp[fTbl[i]])
    end
    
    ylim([12.6,0.4])
    xticks( prepend!(collect(l:l:L),1) )
    yticks( prepend!(collect(l:l:L),1) ) #, append!(collect(L:-l:l),1) )
#     colorbar(fraction=0.02)
end

squareGridPlot(GRIDSIZE, l=4, fsize=4);
savefig("/Users/riccardorao/delocalized/florentia/Finchley/typGrid.svg", bbox_inches="tight")
# savefig("typGridTalk.pdf", bbox_inches="tight")
# -


