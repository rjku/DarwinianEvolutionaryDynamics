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

folderName = "/Users/riccardorao/projects/fisheria/analysis/gridWorld/scratch/rbs/";
# folderName = "/Users/riccardorao/projects/fisheria/analysis/gridWorld/scratch/rbs60/";

# +
# null ID for checking parameters
# aJobID = ["rbsGridEvo"]

# job ID: N = 100
aJobID = ["rbsGridEvo-7796724"] #, "rbsGridEvo-7796725", "rbsGridEvo-7796726", "rbsGridEvo-7796727",
#     "rbsGridEvo-7796729", "rbsGridEvo-7796730"]
# job ID: N = 60
# aJobID = ["rbsGridEvo-7796736", "rbsGridEvo-7796741", "rbsGridEvo-7796742", "rbsGridEvo-7796743",
#     "rbsGridEvo-7796744"]
;

# +
aaMutFactor = Vector{Vector{Float64}}(undef, length(aJobID))

for (iID, jobID) in enumerate(aJobID)
    include(folderName * jobID * "-parameters.jl")
    aaMutFactor[iID] = aMutFactor
end
;
# -

# # Data Generation/Collection

# +
aTraj = [Array{EvolutionaryDynamics.PopulationTrajectoryData}(undef, length(aSelStrength), length(aMutFactor))]
aaMutFactor = [aMutFactor];

@time for (i,β) in enumerate(aSelStrength), (j,μ) in enumerate(aMutFactor)
    grid = mGraphs.EdgeWeightedSquareLattice(GRIDSIZE, [ μ for i in 1:2GRIDSIZE*(GRIDSIZE-1) ]);
    ety = EvolutionaryDynamics.TabularEvotype(
        EvolutionaryDynamics.WithoutReplication(),
        EvolutionaryDynamics.StandardMutation(),
        EvolutionaryDynamics.FitnessSelection(),
        grid
    )

#     traj = @distributed (+) for i in 1:NSAMPLES
    traj = EvolutionaryDynamics.generateTabularSystemsPopulationTrajectories(
            ety=ety, fitnessTbl=fTbl, selCoef=β, Npop=NPOP, nGenRelax=NGENRELAX, nSamples=1 #NSAMPLESPERTRAJ
        )
#     end
    aTraj[1][i,j] = traj
end

using JLD, HDF5, DelimitedFiles

jldopen(folderName * aJobID[1] * "_aTraj.jld", "w") do file
    addrequire(file, EvolutionaryDynamics)
    write(file, "aTraj", aTraj[1])
end

# +
# sensitivity analysis results
using JLD, HDF5, DelimitedFiles

# PopulationTrajectoryData Array: jobID ⊗ (β, μ) parameter values
aTraj = Vector{ Array{EvolutionaryDynamics.TabularSystems.PopulationTrajectoryData, 2} }(undef, length(aJobID) + 1)
for (i, jID) in enumerate(aJobID)
    aTraj[i] = mUtils.readJLD(folderName * jID * "_aTraj.jld", "aTraj")
end
# -

# merging of trajectory data into the last | not to run more than once
aTraj[end] = deepcopy(aTraj[1])
for i in 2:length(aJobID)
    aTraj[end] .+= aTraj[i]
end
push!(aJobID, "merged" * string(NPOP))
push!(aaMutFactor, aaMutFactor[1]);

# +
# tests
# for i in eachindex(aTraj[1])
#     if (aTraj[1][i].mutationFactor != aTraj[2][i].mutationFactor)
#         println(i)
#     end
# end
x = rand(1:length(aJobID)-1)
y = rand( filter( e -> e != x, 1:length(aJobID)-1) )

println( x, ", ", y, " -> ", aTraj[x][1].aPopCmp == aTraj[y][1].aPopCmp)
# -

# #### Multiple Population Sizes Merging

aMergedJobID = [aJobID[end]]
aMergedTraj = [deepcopy(aTraj[end])]
aNpop = [NPOP]

push!( aMergedJobID, aJobID[end] )
push!( aMergedTraj, deepcopy(aTraj[end]) )
push!( aNpop, NPOP)

# # Analysis
#
# ### Checks about relaxation

# +
Base.show(io::IO, f::Float64) = @printf(io, "%.5f", f)

iID = 1
iβ = 1
iμ = 1

subplots(2,1,figsize=(7.5,7.5))
subplot(211); title("Average Fitness");
    plot( aTraj[iID][iβ,iμ].avePerformance );
    title("Average Scaled Growth Coeff. | β = $(aSelStrength[iβ]), μ = $(aMutFactor[iμ])");
subplot(212); title("Mutation Factor");
    plot( aTraj[iID][iβ,iμ].mutationFactor );
tight_layout();
# -

# ### Sensitivity

# +
using LinearAlgebra, Statistics

iβ = 1
β = aSelStrength[iβ]

# vector of fitness values and matrix of fitness changes
vf = [ exp( β * fTbl[g] ) for g in 1:DIMGSPACE ]
mΔf = [ vf[go] - vf[gt] for gt in 1:DIMGSPACE, go in 1:DIMGSPACE ]

# construction of the mutation probability matrix
fnGrid(μ) = mGraphs.EdgeWeightedSquareLattice(GRIDSIZE, [ μ for i in 1:2GRIDSIZE*(GRIDSIZE-1) ])
mMut(μ) = ( Array(mGraphs.transitionMatrix(fnGrid(μ))) .* ( 1.0 .- Diagonal( ones(Float64, DIMGSPACE) ) ) ) / 4μ

# vectors of fitness values after mutations
vfπ(μ) = mMut(μ)' * vf
vΔf(μ) = sum( mΔf .* mMut(μ), dims=1 )[:]

# mσ = [ (vf[go] - vf[gt])/vf[go] for gt in 1:DIMGSPACE, go in 1:DIMGSPACE ]
# mσLog = [ log(vf[go] / vf[gt]) for gt in 1:DIMGSPACE, go in 1:DIMGSPACE ]

# genotype expected loss of growth
# vσ(μ)    = sum( mσ .* mMut(μ), dims=1 )
# vσLog(μ) = sum( mσLog .* mMut(μ), dims=1 )

# +
# sensitivity functions
function σPopVal( vf::AbstractVector, vΔf::AbstractVector, aGenome )
    ΔF, F = 0.0, 0.0

    for genome in aGenome
        ΔF += vΔf[genome]
        F += vf[genome]
    end

    return ΔF / F
end

function σPopStat( vf::AbstractVector, vΔf::AbstractVector,  aaGenome )
    aσ = Vector{Float64}(undef, length(aaGenome))

    for (iSample, aGenome) in enumerate(aaGenome)
        aσ[iSample] = σPopVal( vf, vΔf, aGenome )
    end

    return mean(aσ), var(aσ)
end

# log sensitivity functions
function σLogPopVal( vf::AbstractVector, vfπ::AbstractVector, aGenome )
    Fπ, F = 0.0, 0.0

    for genome in aGenome
        Fπ += vfπ[genome]
        F += vf[genome]
    end

    return log( F / Fπ )
#     return ( F - Fπ )/F
end

function σLogPopStat( vf::AbstractVector, vfπ::AbstractVector,  aaGenome )
    aσ = Vector{Float64}(undef, length(aaGenome))

    for (iSample, aGenome) in enumerate(aaGenome)
        aσ[iSample] = σLogPopVal( vf, vfπ, aGenome )
    end

    return mean(aσ), var(aσ)
end


# +
aExpectedSensitivityPop    = Vector{Vector{Float64}}(undef, length(aJobID))
aExpectedSensitivityPopLog = Vector{Vector{Float64}}(undef, length(aJobID))

aVarianceSensitivityPop    = Vector{Vector{Float64}}(undef, length(aJobID))
aVarianceSensitivityPopLog = Vector{Vector{Float64}}(undef, length(aJobID))

aChangeSensitivityPop    = Vector{Vector{Float64}}(undef, length(aJobID))
aChangeSensitivityPopLog = Vector{Vector{Float64}}(undef, length(aJobID))

aDivergenceSensitivityPop    = Vector{Vector{Float64}}(undef, length(aJobID))
aDivergenceSensitivityPopLog = Vector{Vector{Float64}}(undef, length(aJobID))

for iID in eachindex(aJobID)
    aExpectedSensitivityPop[iID]    = Vector{Float64}(undef, length(aaMutFactor[iID]))
    aExpectedSensitivityPopLog[iID] = Vector{Float64}(undef, length(aaMutFactor[iID]))

    aVarianceSensitivityPop[iID]    = Vector{Float64}(undef, length(aaMutFactor[iID]))
    aVarianceSensitivityPopLog[iID] = Vector{Float64}(undef, length(aaMutFactor[iID]))

    aChangeSensitivityPop[iID]    = Vector{Float64}(undef, length(aaMutFactor[iID])-1)
    aChangeSensitivityPopLog[iID] = Vector{Float64}(undef, length(aaMutFactor[iID])-1)

    aDivergenceSensitivityPop[iID]    = Vector{Float64}(undef, length(aaMutFactor[iID])-1)
    aDivergenceSensitivityPopLog[iID] = Vector{Float64}(undef, length(aaMutFactor[iID])-1)

    for (iμ, μ) in enumerate(aaMutFactor[iID])
        aExpectedSensitivityPop[iID][iμ], aVarianceSensitivityPop[iID][iμ] =
            σPopStat( vf, vΔf(μ), aTraj[iID][iβ,iμ].aPopCmp )
        aExpectedSensitivityPopLog[iID][iμ], aVarianceSensitivityPopLog[iID][iμ] =
            σLogPopStat( vf, vfπ(μ), aTraj[iID][iβ,iμ].aPopCmp )
    end

    for (iμ, μ) in enumerate(aaMutFactor[iID][1:end-1])
        aChangeSensitivityPop[iID][iμ] =
            ( aExpectedSensitivityPop[iID][iμ+1] - aExpectedSensitivityPop[iID][iμ] ) /
            ( - aaMutFactor[iID][iμ+1] + μ ) / NPOP
        aChangeSensitivityPopLog[iID][iμ] =
            ( aExpectedSensitivityPopLog[iID][iμ+1] - aExpectedSensitivityPopLog[iID][iμ] ) /
            ( - aaMutFactor[iID][iμ+1] + μ ) / NPOP

        aDivergenceSensitivityPop[iID][iμ] =
            ( aChangeSensitivityPop[iID][iμ] / aVarianceSensitivityPop[iID][iμ] )
        aDivergenceSensitivityPopLog[iID][iμ] =
            ( aChangeSensitivityPopLog[iID][iμ] / aVarianceSensitivityPopLog[iID][iμ] )
    end
end
# -

using Interpolations, Dierckx

# +
rng = 1:length(aaMutFactor[end])

itpAveSnt = Spline1D( aaMutFactor[end][rng], aExpectedSensitivityPop[end][rng];
    w=ones(length(aaMutFactor[end][rng])), k=5, bc="nearest", s=10^(-4) )

itpVarSnt = Spline1D( aaMutFactor[end][rng], aVarianceSensitivityPop[end][rng];
    w=ones(length(aaMutFactor[end][rng])), k=5, bc="nearest", s=10^(-4.6) )

# itpDevSnt = ( - derivative(itpAveSnt, aaMutFactor[end][rng]) / NPOP ) ./ aVarianceSensitivityPop[end][rng]
devSnt = ( - derivative(itpAveSnt, aaMutFactor[end][rng]) / NPOP ) ./ itpVarSnt(aaMutFactor[end][rng])

# plot( aaMutFactor[end], itpAveSnt(aaMutFactor[end]) )
# plot( aaMutFactor[end], itpVarSnt(aaMutFactor[end]) )

subplots(1, 3, figsize=(15.0, 4.0))
subplot(131); plot( aaMutFactor[end][rng], itpAveSnt(aaMutFactor[end][rng]) ); xscale("log")
subplot(132); plot( aaMutFactor[end][rng], itpVarSnt(aaMutFactor[end][rng]) ); xscale("log")
subplot(133); plot( aaMutFactor[end][rng], devSnt[rng] ); axhline(y=1, linestyle="--"); xscale("log"); ylim([-1,10])

# +
subplots(3, 2, figsize=(10.0, 10.0))

subplot(321); title("(a) Average Population Sensitivity");
    for (iID, EΣ) in enumerate(aExpectedSensitivityPop[end:end])
        errorbar( aaMutFactor[end], EΣ, sqrt.(aVarianceSensitivityPop[end]), fmt="o", fillstyle="none" ) #,
#             color="tab:blue" )
    end
    plot( aaMutFactor[end], itpAveSnt(aaMutFactor[end]), color="tab:blue", linestyle="--" )
#     axhline(y=0, linestyle="--")
    xscale("log")
    ylabel("Average Sensitivity, ⟨Ωₙ⟩")

subplot(322); title("(b) Average Population Log Sensitivity");
    for (iID, EΣ) in enumerate(aExpectedSensitivityPopLog[end:end])
        errorbar( aaMutFactor[end], EΣ, sqrt.(aVarianceSensitivityPop[end]), fmt="o", fillstyle="none" )
    end
#     axhline(y=1, linestyle="--")
    ylabel("Average Log Sensitivity, ⟨ΔlnFₙ⟩")
    xscale("log")

subplot(323); title("(c) Variance of Population Sensitivity");
    for (iID, c) in enumerate(aVarianceSensitivityPop[end:end])
        errorbar( aaMutFactor[end], c, fmt="o", fillstyle="none" )
    end
    plot( aaMutFactor[end], itpVarSnt(aaMutFactor[end]), color="tab:blue", linestyle="--" )
    xscale("log")
    ylabel("Variance Sensitivity, Var{Ωₙ}")

subplot(324); title("(d) Variance Population Log Sensitivity");
    for (iID, c) in enumerate(aVarianceSensitivityPopLog[end:end])
        errorbar( aaMutFactor[end], c, fmt="o", fillstyle="none" )
    end
    xscale("log")
    ylabel("Variance Log Sensitivity, Var{ΔlnFₙ}")

# ylims = [0.01, 10^4]
subplot(325); title("(e) Deviations from Fluctuations Relation");
    for (iID, c) in enumerate(aDivergenceSensitivityPop[end:end])
#     for (iID, c) in enumerate(aChangeSensitivityPopNeutral)
        errorbar( aaMutFactor[end][1:end-1], c, fmt="o", fillstyle="none" )
    end
#     plot( aaMutFactor[end], itpDevSnt, color="tab:red" )
    axhline(y=1, linestyle="--")
    xscale("log")
    xlabel("variation probability, μ̄")
    ylabel("deviation, δ")
#     yscale("log")
#     ylim(ylims)

subplot(326); title("(f) Deviations from Fluctuations Relation");
    for (iID, c) in enumerate(aDivergenceSensitivityPopLog[end:end])
    #     for (iID, c) in enumerate(aChangeSensitivityPop)
        errorbar( aaMutFactor[end][1:end-1], c, fmt="o", fillstyle="none" )
    end
    axhline(y=1, linestyle="--")
    xscale("log")
    xlabel("variation probability, μ̄")
#     ylabel("deviation, δ")
#     yscale("log")
#     ylim(ylims)

# suptitle("Population Sensitivity Statistics: N = 60", fontsize=16)
tight_layout();

savefig("snt100.pdf", bbox_inches="tight")
# -
# #### Multiple Population Size Analysis

# +
aExpectedSensitivityPop    = Vector{Vector{Float64}}(undef, length(aMergedJobID))
aExpectedSensitivityPopLog = Vector{Vector{Float64}}(undef, length(aMergedJobID))

aVarianceSensitivityPop    = Vector{Vector{Float64}}(undef, length(aMergedJobID))
aVarianceSensitivityPopLog = Vector{Vector{Float64}}(undef, length(aMergedJobID))

aChangeSensitivityPop    = Vector{Vector{Float64}}(undef, length(aMergedJobID))
aChangeSensitivityPopLog = Vector{Vector{Float64}}(undef, length(aMergedJobID))

aDivergenceSensitivityPop    = Vector{Vector{Float64}}(undef, length(aMergedJobID))
aDivergenceSensitivityPopLog = Vector{Vector{Float64}}(undef, length(aMergedJobID))

for iID in eachindex(aMergedJobID)
    aExpectedSensitivityPop[iID]    = Vector{Float64}(undef, length(aaMutFactor[iID]))
    aExpectedSensitivityPopLog[iID] = Vector{Float64}(undef, length(aaMutFactor[iID]))

    aVarianceSensitivityPop[iID]    = Vector{Float64}(undef, length(aaMutFactor[iID]))
    aVarianceSensitivityPopLog[iID] = Vector{Float64}(undef, length(aaMutFactor[iID]))

    aChangeSensitivityPop[iID]    = Vector{Float64}(undef, length(aaMutFactor[iID])-1)
    aChangeSensitivityPopLog[iID] = Vector{Float64}(undef, length(aaMutFactor[iID])-1)

    aDivergenceSensitivityPop[iID]    = Vector{Float64}(undef, length(aaMutFactor[iID])-1)
    aDivergenceSensitivityPopLog[iID] = Vector{Float64}(undef, length(aaMutFactor[iID])-1)

    for (iμ, μ) in enumerate(aaMutFactor[iID])
        aExpectedSensitivityPop[iID][iμ], aVarianceSensitivityPop[iID][iμ] =
            σPopStat( vf, vΔf(μ), aTraj[iID][iβ,iμ].aPopCmp )
        aExpectedSensitivityPopLog[iID][iμ], aVarianceSensitivityPopLog[iID][iμ] =
            σLogPopStat( vf, vfπ(μ), aTraj[iID][iβ,iμ].aPopCmp )
    end

    for (iμ, μ) in enumerate(aaMutFactor[iID][1:end-1])
        aChangeSensitivityPop[iID][iμ] =
            ( aExpectedSensitivityPop[iID][iμ+1] - aExpectedSensitivityPop[iID][iμ] ) /
            ( - aaMutFactor[iID][iμ+1] + μ ) / aNpop[iID]
        aChangeSensitivityPopLog[iID][iμ] =
            ( aExpectedSensitivityPopLog[iID][iμ+1] - aExpectedSensitivityPopLog[iID][iμ] ) /
            ( - aaMutFactor[iID][iμ+1] + μ ) / aNpop[iID]

        aDivergenceSensitivityPop[iID][iμ] =
            ( aChangeSensitivityPop[iID][iμ] / aVarianceSensitivityPop[iID][iμ] )
        aDivergenceSensitivityPopLog[iID][iμ] =
            ( aChangeSensitivityPopLog[iID][iμ] / aVarianceSensitivityPopLog[iID][iμ] )
    end
end

# +
subplots(3, 2, figsize=(10.0, 10.0))

subplot(321); title("(a) Average Population Sensitivity");
    for (iID, EΣ) in enumerate(aExpectedSensitivityPop)
        errorbar( aaMutFactor[iID], EΣ, sqrt.(aVarianceSensitivityPop[iID]), fmt="o", fillstyle="none" ) #,
#             color="tab:blue" )
    end
#     axhline(y=0, linestyle="--")
    xscale("log")
    ylabel("Average Sensitivity, ⟨Ωₙ⟩")

subplot(322); title("(b) Average Population Log Sensitivity");
    for (iID, EΣ) in enumerate(aExpectedSensitivityPopLog)
        errorbar( aaMutFactor[iID], EΣ, sqrt.(aVarianceSensitivityPop[iID]), fmt="o", fillstyle="none" )
    end
#     axhline(y=1, linestyle="--")
    ylabel("Average Log Sensitivity, ⟨ΔlnFₙ⟩")
    xscale("log")

subplot(323); title("(c) Variance of Population Sensitivity");
    for (iID, c) in enumerate(aVarianceSensitivityPop)
        errorbar( aaMutFactor[iID], c, fmt="o", fillstyle="none" )
    end
    xscale("log")
    ylabel("Variance Sensitivity, Var{Ωₙ}")

subplot(324); title("(d) Variance Population Log Sensitivity");
    for (iID, c) in enumerate(aVarianceSensitivityPopLog)
        errorbar( aaMutFactor[iID], c, fmt="o", fillstyle="none" )
    end
    xscale("log")
    ylabel("Variance Log Sensitivity, Var{ΔlnFₙ}")

# ylims = [0.01, 10^4]
subplot(325); title("(e) Deviations from Fluctuations Relation");
    for (iID, c) in enumerate(aDivergenceSensitivityPop)
#     for (iID, c) in enumerate(aChangeSensitivityPopNeutral)
        errorbar( aaMutFactor[iID][1:end-1], c, fmt="o", fillstyle="none" )
    end
    axhline(y=1, linestyle="--")
    xscale("log")
    xlabel("variation probability, μ̄")
    ylabel("deviation, δ")
#     yscale("log")
#     ylim(ylims)

subplot(326); title("(f) Deviations from Fluctuations Relation");
    for (iID, c) in enumerate(aDivergenceSensitivityPopLog)
    #     for (iID, c) in enumerate(aChangeSensitivityPop)
        errorbar( aaMutFactor[iID][1:end-1], c, fmt="o", fillstyle="none" )
    end
    axhline(y=1, linestyle="--")
    xscale("log")
    xlabel("variation probability, μ̄")
#     ylabel("deviation, δ")
#     yscale("log")
#     ylim(ylims)

# suptitle("Population Sensitivity Statistics: N = 60", fontsize=16)
tight_layout();

savefig("snt100.pdf", bbox_inches="tight")
# -
aExpectedSensitivityPop[4] - aExpectedSensitivityPop[2]

# +
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
Esensitivity(μ, n=1) = sum( mSensitivity .* mMut(μ) .* prob(μ, n) ) / 4μ
EsensitivityNeutral(μ, n=1) = sum( mSensitivity .* mMut(μ) .* probNeutral(μ, n) ) / 4μ

EsensitivityInf(μ) = sum( mSensitivity .* mMut(μ) .* probInf(μ) ) / 4μ
EsensitivityInfNeutral(μ) = sum( mSensitivity .* mMut(μ) .* probInfNeutral(μ) ) / 4μ

# average fitness values
aveF(μ, n=1) = prob(μ, n) ⋅ af
aveFInf(μ) = probInf(μ) ⋅ af
;

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
