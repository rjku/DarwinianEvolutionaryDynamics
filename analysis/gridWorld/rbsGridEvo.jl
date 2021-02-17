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
#     display_name: Julia 1.5.3
#     language: julia
#     name: julia-1.5
# ---

versioninfo()

# +
using Revise
# using ProgressMeter, BenchmarkTools
using JLD, HDF5 #, DelimitedFiles
using PyPlot, Printf
using Dierckx, LsqFit, ForwardDiff
using DataFrames, GLM

import EvolutionaryDynamics, mUtils, mPlot, mGraphs

folderName = "/Users/riccardorao/projects/fisheria/analysis/gridWorld/scratch/rbs/";
# folderName = "/Users/riccardorao/projects/fisheria/analysis/gridWorld/scratch/rbs60/";


# +
# null ID for checking parameters
# aJobID = ["rbsGridEvo"]

# job ID: N = 100
# aJobID = ["rbsGridEvo-7796724"] #, "rbsGridEvo-7796725", "rbsGridEvo-7796726", "rbsGridEvo-7796727",
#     "rbsGridEvo-7796729", "rbsGridEvo-7796730"]
# job ID: N = 70
aJobID = [
    "rbsGridEvo-7796924", "rbsGridEvo-7796926", "rbsGridEvo-7796928", "rbsGridEvo-7796929", "rbsGridEvo-7796931",
    "rbsGridEvo-7796932", "rbsGridEvo-7796933", "rbsGridEvo-7796934", "rbsGridEvo-7796937", "rbsGridEvo-7796938",
    "rbsGridEvo-7797282", "rbsGridEvo-7797283", "rbsGridEvo-7797284", "rbsGridEvo-7797285", "rbsGridEvo-7797286",
    "rbsGridEvo-7797318", "rbsGridEvo-7797319", "rbsGridEvo-7797320", "rbsGridEvo-7797321", "rbsGridEvo-7797322",
    "rbsGridEvo-7920547", "rbsGridEvo-7920548", "rbsGridEvo-7920549", "rbsGridEvo-7920550", "rbsGridEvo-7920552"
]
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

# ## Data Generation

# +
aTraj = [Array{EvolutionaryDynamics.PopulationTrajectoryData}(undef, length(aSelStrength), length(aMutFactor))]
aaMutFactor = [aMutFactor];

@time for (i, β) in enumerate(aSelStrength), (j, μ) in enumerate(aMutFactor)
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

jldopen(folderName * aJobID[1] * "_aTraj.jld", "w") do file
    addrequire(file, EvolutionaryDynamics)
    write(file, "aTraj", aTraj[1])
end
# -

# ## Data Retrival

# +
# sensitivity analysis results

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

println( "? are (x, y) wrongly equal: (", x, ", ", y, ") -> ", aTraj[x][1].aPopCmp == aTraj[y][1].aPopCmp)
# -

?EvolutionaryDynamics.TabularSystems.PopulationTrajectoryData

# +
# using DelimitedFiles

# for k in eachindex(aMutFactor)
#     open("rbsPopTrajData/rbsGridEvo_aPopCmp_" * string(k) * ".dat", "w") do io
#         writedlm(io, [ aTraj[end][k].aPopCmp[j][i]
#                 for i in eachindex(aTraj[end][k].aPopCmp[1]),
#                     j in eachindex(aTraj[end][k].aPopCmp) ]
#         )
#     end;
# end

# readdlm("rbsPopTrajData/rbsGridEvo_aPopCmp_1.dat", '\t', Int, '\n')
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
subplot(211); title("Average Scaled Growth Coeff. | β = $(aSelStrength[iβ]), μ = $(aMutFactor[iμ])");
    plot( aTraj[iID][iβ,iμ].avePerformance );
subplot(212); title("Mutation Factor: ⟨mf⟩ = $(mean(aTraj[iID][iβ,iμ].mutationFactor))");
    plot( aTraj[iID][iβ,iμ].mutationFactor );
#     yscale("log")
tight_layout();
# -

# ## Sensitivity

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

    return log(F/Fπ)
#     return 1.0 - Fπ/F
#     return F/Fπ
end

function σLogPopStat( vf::AbstractVector, vfπ::AbstractVector,  aaGenome )
    aσ = Vector{Float64}(undef, length(aaGenome))

    for (iSample, aGenome) in enumerate(aaGenome)
        aσ[iSample] = σLogPopVal( vf, vfπ, aGenome )
    end

    return mean(aσ), var(aσ)
end


# +
# evaluation of sensitivity and related statistics for all single samples plus the merged one

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

# +
# evaluation of sensitivity and related statistics for all samples combined

aExpectedSensitivityPop    = [Vector{Float64}(undef, length(aaMutFactor[end]))]
aExpectedSensitivityPopLog = [Vector{Float64}(undef, length(aaMutFactor[end]))]

aVarianceSensitivityPop    = [Vector{Float64}(undef, length(aaMutFactor[end]))]
aVarianceSensitivityPopLog = [Vector{Float64}(undef, length(aaMutFactor[end]))]

aChangeSensitivityPop    = [Vector{Float64}(undef, length(aaMutFactor[end])-1)]
aChangeSensitivityPopLog = [Vector{Float64}(undef, length(aaMutFactor[end])-1)]

aDivergenceSensitivityPop    = [Vector{Float64}(undef, length(aaMutFactor[end])-1)]
aDivergenceSensitivityPopLog = [Vector{Float64}(undef, length(aaMutFactor[end])-1)]

for (iμ, μ) in enumerate(aaMutFactor[end])
    aExpectedSensitivityPop[end][iμ], aVarianceSensitivityPop[end][iμ] =
        σPopStat( vf, vΔf(μ), aTraj[end][iβ,iμ].aPopCmp )
    aExpectedSensitivityPopLog[end][iμ], aVarianceSensitivityPopLog[end][iμ] =
        σLogPopStat( vf, vfπ(μ), aTraj[end][iβ,iμ].aPopCmp )
end

for (iμ, μ) in enumerate(aaMutFactor[end][1:end-1])
    aChangeSensitivityPop[end][iμ] =
        ( aExpectedSensitivityPop[end][iμ+1] - aExpectedSensitivityPop[end][iμ] ) /
        ( - aaMutFactor[end][iμ+1] + μ ) / NPOP
    aChangeSensitivityPopLog[end][iμ] =
        ( aExpectedSensitivityPopLog[end][iμ+1] - aExpectedSensitivityPopLog[end][iμ] ) /
        ( - aaMutFactor[end][iμ+1] + μ ) / NPOP

    aDivergenceSensitivityPop[end][iμ] =
        ( aChangeSensitivityPop[end][iμ] / aVarianceSensitivityPop[end][iμ] )
    aDivergenceSensitivityPopLog[end][iμ] =
        ( aChangeSensitivityPopLog[end][iμ] / aVarianceSensitivityPopLog[end][iμ] )
end
# -

# #### Interpolation Approach

# +
iμMin, iμMax, iμΔ = 2, length(aaMutFactor[end]), 1 # <--- 2, 0, 1
rng, rngDev, rngDevMu = iμMin:iμΔ:iμMax, iμMin:iμΔ:iμMax-1, iμMin:iμΔ:iμMax-1

# neutral weights
wVec = ones(Float64,length(aaMutFactor[end][rng])) # <---
# more weight for low μ
# wVec = [ 1.0 - 0.1 * ( x - rng[1] ) / ( rng[end] - rng[1] ) for x in rng ]
# more weight for high μ
# wVec = [ 2.0( x - rng[1] ) / ( rng[end] - rng[1] ) for x in rng ] 
# more weight for middle-valued μ
# wVec = [ exp(- ( x - rng[ trunc(Int32, length(rng)/2) ] )^2 / 50 ) for x in rng ] 

# wVec[7] = 0.0

# wVec = 1 ./ sqrt.(aVarianceSensitivityPop[end][rng])
wVec[10:11] .= 10    # <--- 10

itpAveSnt = Spline1D( aaMutFactor[end][rng], aExpectedSensitivityPop[end][rng];
    w=wVec, k=1, bc="error", s=10^(-4.0) )

# itpAveSnt = Spline1D( aaMutFactor[end][rng], aExpectedSensitivityPop[end][rng];
#     w=wVec, k=1, bc="nearest", s=10^(-3) )

itpVarSnt = Spline1D( aaMutFactor[end][rng], aVarianceSensitivityPop[end][rng];
    w=wVec, k=1, bc="nearest", s=10^(-5.0) )

devSnt = ( - derivative(itpAveSnt, aaMutFactor[end][rng]) / NPOP ) ./ aVarianceSensitivityPop[end][rng]
devSntItp = ( - derivative(itpAveSnt, aaMutFactor[end][rng]) / NPOP ) ./ itpVarSnt(aaMutFactor[end][rng])

# plot( aaMutFactor[end], itpAveSnt(aaMutFactor[end]) )
# plot( aaMutFactor[end], itpVarSnt(aaMutFactor[end]) )

subplots(1, 3, figsize=(15.0, 4.0))
subplot(131);
    errorbar( 4aaMutFactor[end], aExpectedSensitivityPop[end], fmt="o", fillstyle="none" )
    plot( 4aaMutFactor[end][rng], itpAveSnt(aaMutFactor[end][rng]), "-", color="tab:red" );
for x in get_knots( itpAveSnt )
    axvline(x=4x, linestyle="dotted", color="tab:grey")
end
    xscale("log")
subplot(132);
#     plot( 4aaMutFactor[end][rng], itpVarSnt(aaMutFactor[end][rng]) );
    plot( 4aaMutFactor[end][rng], derivative(itpAveSnt, aaMutFactor[end][rng]), "-" );
    xscale("log")
subplot(133);
    errorbar( 4aaMutFactor[end][rngDevMu], aDivergenceSensitivityPop[end][rngDev], fmt="o", fillstyle="none" )
    plot( 4aaMutFactor[end][rng], devSnt, "*", color="tab:red");
#     errorbar( aaMutFactor[end][rngDev], aDivergenceSensitivityPop[end][rngDev], fmt="o", fillstyle="none" )
#     plot( 4aaMutFactor[end][rng], devSntItp, ".");
    axvline(x=1/NPOP, ymax=1.0, linestyle="dotted", color="tab:grey")
    axhline(y=1, linestyle="--", color="tab:grey");
    xscale("log");
ylim([-4,10]);

[ get_knots( itpAveSnt ) get_coeffs( itpAveSnt ) ]
# -

for (x, y) in zip( get_knots( itpAveSnt ), get_coeffs( itpAveSnt ) )
    println("{$(log10(x)), $(10y)},")
end

for (x, y) in zip(log10.(aaMutFactor[end][rng]), aExpectedSensitivityPop[end][rng])
    println("{$x, $(10y)},")
end

# +
iμMin, iμMax, iμΔ = 2, length(aaMutFactor[end]), 1
rng, rngDev, rngDevMu = iμMin:iμΔ:iμMax, iμMin:iμΔ:iμMax-1, iμMin:iμΔ:iμMax-1

# neutral weights
wVec = ones(Float64,length(aaMutFactor[end][rng]))
# more weight for low μ
# wVec = [ 1.0 - 0.1 * ( x - rng[1] ) / ( rng[end] - rng[1] ) for x in rng ]
# more weight for high μ
# wVec = [ 2.0( x - rng[1] ) / ( rng[end] - rng[1] ) for x in rng ] 
# more weight for middle-valued μ
# wVec = [ exp(- ( x - rng[ trunc(Int32, length(rng)/2) ] )^2 / 50 ) for x in rng ] 

# wVec[7] = 0.0
wVec[10:11] .= 10

itpAveSntLog = Spline1D( aaMutFactor[end][rng], aExpectedSensitivityPopLog[end][rng];
    w=wVec, k=1, bc="nearest", s=10^(-3.6) )

itpVarSntLog = Spline1D( aaMutFactor[end][rng], aVarianceSensitivityPopLog[end][rng];
    w=wVec, k=1, bc="nearest", s=10^(-3.0) )

devSntLog = ( - derivative(itpAveSntLog, aaMutFactor[end][rng]) / NPOP ) ./ aVarianceSensitivityPopLog[end][rng]
devSntItpLog = ( - derivative(itpAveSntLog, aaMutFactor[end][rng]) / NPOP ) ./ itpVarSntLog(aaMutFactor[end][rng])

# plot( aaMutFactor[end], itpAveSnt(aaMutFactor[end]) )
# plot( aaMutFactor[end], itpVarSnt(aaMutFactor[end]) )

subplots(1, 3, figsize=(15.0, 4.0))
subplot(131);
    errorbar( 4aaMutFactor[end], aExpectedSensitivityPopLog[end], fmt="o", fillstyle="none" )
    plot( 4aaMutFactor[end][rng], itpAveSntLog(aaMutFactor[end][rng]), "--" );
    xscale("log")
subplot(132); plot( 4aaMutFactor[end][rng], itpVarSntLog(aaMutFactor[end][rng]) ); xscale("log")
subplot(133);
    errorbar( 4aaMutFactor[end][rngDevMu], aDivergenceSensitivityPopLog[end][rngDev], fmt="o", fillstyle="none" )
    plot( 4aaMutFactor[end][rng], devSntLog, "*", color="tab:red");
#     errorbar( aaMutFactor[end][rngDev], aDivergenceSensitivityPop[end][rngDev], fmt="o", fillstyle="none" )
#     plot( aaMutFactor[end][rng], devSntItpLog, ".");
    axhline(y=1, linestyle="--", color="tab:grey");
    xscale("log");
# ylim([-4,10]);

# +
dataDiv = DataFrame(X=log.(aaMutFactor[end][rngDev]), Y=aDivergenceSensitivityPop[end][rngDev])
lmDiv = lm(@formula(Y ~ X), dataDiv)

dataDivItp = DataFrame(X=log.(aaMutFactor[end][rng]), Y=devSntItp)
lmDivItp = lm(@formula(Y ~ X), dataDivItp)

r2(lmDiv), r2(lmDivItp)
# -

# #### Sigmoidal Fit Approach

# +
iμMin, iμMax, iμΔ = 1, length(aaMutFactor[end]), 1 # length(aaMutFactor[end])
rng, rngDev, rngDevMu = iμMin:iμΔ:iμMax, iμMin:iμΔ:iμMax-1, iμMin:iμΔ:iμMax-1

μVals = aaMutFactor[end][rng]
ΩVals = aExpectedSensitivityPop[end][rng]
varΩVals = aVarianceSensitivityPop[end][rng]

wVec = (varΩVals).^(-1/2)
wVec[9:12] .= 100.0
# wVec = [ 1.0 - 0.1 * ( x - rng[1] ) / ( rng[end] - rng[1] ) for x in rng ]
# wVec = ones(Float64,length(aaMutFactor[end][rng]))
# wVec = zeros(Float64,length(aaMutFactor[end][rng]))

@. polyRatio3Model(μ, p) = ( p[1] + p[2] / μ + p[3] / μ^2 ) / ( p[4] + p[5] / μ + p[6] / μ^2 )
@. polyRatioModel(μ, p) = ( p[1] + p[2] / μ ) / ( p[3] + p[4] / μ )

@. tanhModel(μ, p) = p[1] - p[2] * tanh( p[3] + p[4] * μ )
@. arctanModel(μ, p) = p[1] - p[2] * atan( p[3] + p[4] * μ )
@. gdModel(μ, p) = p[1] - p[2] * atan( tanh( p[3] + p[4] * μ ) )
    algebraicSigmoid(x) = x / sqrt( 1 + x^2 )
@. algebraicModel(μ, p) = p[1] - p[2] * algebraicSigmoid( p[3] * μ )

prmPolyRatio3 = [10.0, 0.3/10.0^2, 0.1, 10.0, 10.0^-2, 0.1]
prmPolyRatio  = [10.0, 0.3/10.0^2, 10.0, 10.0^2]

prmTanh   = [ 0.15, 0.30, 0.01, 0.1 ]
prmArctan = [ 0.15, 0.30, 0.01, 0.1 ]

# fitAveSnt = LsqFit.curve_fit(polyRatioModel, μVals, ΩVals, wVec, prmPolyRatio)
# fitΩ(μ) = polyRatioModel(μ, coef(fitAveSnt))

# fitAveSnt = LsqFit.curve_fit(tanhModel, μVals, ΩVals, wVec, prmTanh)
# fitΩ(μ) = tanhModel(μ, coef(fitAveSnt))

# fitAveSnt = LsqFit.curve_fit(arctanModel, μVals, ΩVals, wVec, prmTanh)
# fitΩ(μ) = arctanModel(μ, coef(fitAveSnt))

fitAveSnt = LsqFit.curve_fit(gdModel, μVals, ΩVals, wVec, prmTanh)
fitΩ(μ) = gdModel(μ, coef(fitAveSnt))

# fitAveSnt = LsqFit.curve_fit(algebraicModel, μVals, ΩVals, wVec, prmTanh)
# fitΩ(μ) = algebraicModel(μ, coef(fitAveSnt))


prmArctanVar = [ 0.04, 0.04, 0.01, 0.1 ]

# fitVarSnt = LsqFit.curve_fit(sigmoidModel2, μVals, varΩVals, p2)
# fitVarΩ(μ) = sigmoidModel2(μ, coef(fitVarSnt))

# fitVarSnt = LsqFit.curve_fit(arctanModel, μVals, varΩVals, prmArctanVar)
# fitVarΩ(μ) = arctanModel(μ, coef(fitVarSnt))


devSnt = ( - [ ForwardDiff.derivative(fitΩ, μ) for μ in μVals ] / NPOP ) ./ varΩVals
# devSntFit = ( - [ ForwardDiff.derivative(fitΩ, μ) for μ in μVals ] / NPOP ) ./ fitVarΩ.(μVals)

subplots(1, 3, figsize=(15.0, 4.0))
subplot(131);
    errorbar( 4aaMutFactor[end], aExpectedSensitivityPop[end], fmt="o", fillstyle="none" )
    plot( 4μVals, fitΩ.(μVals), "--");
    xscale("log")
subplot(132);
    errorbar( aaMutFactor[end], aVarianceSensitivityPop[end], fmt="o", fillstyle="none" );
#     plot( μVals, fitVarΩ.(μVals), "--" );
    xscale("log")
subplot(133); axhline(y=1, linestyle="--");
    errorbar( 4aaMutFactor[end][rngDevMu], aDivergenceSensitivityPop[end][rngDev], fmt="o", fillstyle="none" )
#     plot( μVals, devSntFit, ".");
    plot( 4μVals, devSnt, "+");
    axhline(y=1, linestyle="--", color="tab:grey");
    axvline(x=1/NPOP, ymax=1.0, linestyle="dotted", color="tab:grey")
    xscale("log");
    ylim([-2,8]);

coef(fitAveSnt) #, coef(fitVarSnt)
# -

# #### Plots

# +
rngPlt = rngDev
rngPltFit = rngDev
rngPltDev = rngDev
rngPltItp = 1:iμMax-iμMin

# subplots(3, 2, figsize=(10.0, 10.0))
subplots(2, 2, figsize=(10.0, 7.0))


subplot(221); title("(a) Average Population Sensitivity");
for (iID, EΣ) in enumerate(aExpectedSensitivityPop[end:end])
    errorbar( 4aaMutFactor[end][rngPlt], EΣ[rngPlt], sqrt.(aVarianceSensitivityPop[end][rngPltDev]),
        fmt="o", fillstyle="none" )
end
# plot( 4aaMutFactor[end][rngPltFit], fitΩ.(aaMutFactor[end][rngPltFit]), color="tab:red", linestyle="--" )
plot( 4aaMutFactor[end][rngPltFit], itpAveSnt(aaMutFactor[end][rngPltFit]), color="tab:red", linestyle="--" )
xscale("log")
xlabel("variation probability, μ")
ylabel("Average Sensitivity, ⟨Ωₙ⟩")

subplot(224); title("(d) Average Population Log Sensitivity");
for (iID, EΣ) in enumerate(aExpectedSensitivityPopLog[end:end])
    errorbar( 4aaMutFactor[end][rngPlt], EΣ[rngPlt], sqrt.(aVarianceSensitivityPopLog[end][rngPlt]),
        fmt="o", fillstyle="none" )
end
xlabel("variation probability, μ")
ylabel("Average Log Sensitivity, ⟨ΔlnFₙ⟩")
xscale("log")


subplot(222); title("(b) Variance of Population Sensitivity");
for (iID, c) in enumerate(aVarianceSensitivityPop[end:end])
    errorbar( 4aaMutFactor[end][rngPlt], c[rngPlt], fmt="o", fillstyle="none" )
end
# plot( aaMutFactor[end][rngPltDev], fitVarΩ.(aaMutFactor[end][rngPltDev]), color="tab:red", linestyle="--" )
xscale("log")
yscale("log")
xlabel("variation probability, μ")
ylabel("Variance Sensitivity, Var{Ωₙ}")
# 
# subplot(324); title("(d) Variance Population Log Sensitivity");
# for (iID, c) in enumerate(aVarianceSensitivityPopLog[end:end])
#     errorbar( aaMutFactor[end][rng], c[rng], fmt="o", fillstyle="none" )
# end
# xscale("log")
# yscale("log")
# ylabel("Variance Log Sensitivity, Var{ΔlnFₙ}")


subplot(223); title("(c) Deviations from Fluctuations Relation");
for (iID, c) in enumerate(aDivergenceSensitivityPop[end:end])
    errorbar( 4aaMutFactor[end][rngPltDev], c[rngPltDev], fmt="o", fillstyle="none" )
end
# # plot( aaMutFactor[end][rngDev], predict(lmDiv), color="tab:blue", linestyle="-." )
plot( 4aaMutFactor[end][rngPltDev], devSnt[rngPltItp], color="tab:red", "*" )
# plot( aaMutFactor[end][rng], predict(lmDivItp), color="tab:red", linestyle="-." )
axhline(y=1, linestyle="--", color="tab:grey")
axvline(x=1/NPOP, ymax=1.0, linestyle="dotted", color="tab:grey")
# text(1/(NPOP+15), 7.7, "1/N")
# # axvline(x=10/NPOP, linestyle="--", color="tab:grey")
xscale("log")
xlabel("variation probability, μ")
ylabel("deviation, δ")
ylim([-1,6]);

# subplot(222); title("(f) Deviations from Fluctuations Relation");
# for (iID, c) in enumerate(aDivergenceSensitivityPopLog[end:end])
#     errorbar( 4aaMutFactor[end][1:end-1], c, fmt="o", fillstyle="none" )
# end
# axhline(y=1, linestyle="--", color="tab:grey")
# xscale("log")
# xlabel("variation probability, μ̄")
# ylabel("deviation, δ")

# suptitle("Population Sensitivity Statistics: N = 60", fontsize=16)
tight_layout();

# savefig("snt" * string(NPOP) * ".pdf", bbox_inches="tight")
# +
fit2clr = Dict( 1.0 => "tab:blue", 5.0 => "tab:red" )
# fit2clr = Dict( 1.0 => "tab:grey", 5.0 => "tab:orange" )
fit2shp = Dict( 1.0 => "s", 5.0 => "s" )
fit2siz = Dict( 1.0 => 130, 5.0 => 230 )

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
# savefig("typGrid.pdf", bbox_inches="tight")
# savefig("typGridTalk.pdf", bbox_inches="tight")
# -

using PyPlot

fig, ax = subplots(1,1)
ax.plot(collect(1:10))
axin1 = ax.inset_axes([0.8, 0.1, 0.15, 0.15])

# +
rngPlt = rngDev
rngPltFit = rngDev
rngPltDev = rngDev
rngPltItp = 1:iμMax-iμMin

# subplots(3, 2, figsize=(10.0, 10.0))
fig, ax = subplots(1, 1, figsize=(5.0, 3.5))

# subplot(111); title("Average Population Sensitivity");
for (iID, EΣ) in enumerate(aExpectedSensitivityPop[end:end])
    ax.errorbar( 4aaMutFactor[end][rngPlt], EΣ[rngPlt], sqrt.(aVarianceSensitivityPop[end][rngPltDev]),
        fmt="o", fillstyle="none" )
end
# plot( 4aaMutFactor[end][rngPltFit], fitΩ.(aaMutFactor[end][rngPltFit]), color="tab:red", linestyle="--" )
# plot( 4aaMutFactor[end][rngPltFit], itpAveSnt(aaMutFactor[end][rngPltFit]), color="tab:red", linestyle="--" )
xscale("log")
xlabel("variation probability, μ")
ylabel("Average Sensitivity, ⟨Ωₙ⟩")

axin = ax.inset_axes([0.62, 0.62, 0.35, 0.33])
# for (iID, c) in enumerate(aDivergenceSensitivityPop[end:end])
#     errorbar( 4aaMutFactor[end][rngPltDev], c[rngPltDev], fmt="o", fillstyle="none" )
# end
# # plot( aaMutFactor[end][rngDev], predict(lmDiv), color="tab:blue", linestyle="-." )
axin.plot( 4aaMutFactor[end][rngPltDev], devSnt[rngPltItp], color="tab:red", "*" )
# # plot( aaMutFactor[end][rng], predict(lmDivItp), color="tab:red", linestyle="-." )
axin.axhline(y=1, linestyle="--", color="tab:grey")
# axvline(x=1/NPOP, ymax=1.0, linestyle="dotted", color="tab:grey")
# text(1/(NPOP-8), 0.2, "← 1/N")
# # axvline(x=10/NPOP, linestyle="--", color="tab:grey")
axin.semilogx()
axin.set_xlabel("μ", fontsize=8, labelpad=-3)
axin.set_ylabel("deviation, δ", fontsize=8, labelpad=0)
axin.set_ylim([0,4.5]);
axin.tick_params(axis="x", labelsize=8)
axin.tick_params(axis="y", labelsize=8)

tight_layout();
savefig("snt" * string(NPOP) * "mt.pdf", bbox_inches="tight")
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
# -

# ### Diversity

# +
# sensitivity functions
function dPopVal( Gsize::Integer, aGenome )
    
    θpop = zeros(Int8, Gsize)

    for genome in aGenome
        θpop[genome] = Int8(1)
    end

    return sum(θpop)
end

function dPopStat( Gsize::Integer,  aaGenome )
    ad = Vector{Int32}(undef, length(aaGenome))

    for (iSample, aGenome) in enumerate(aaGenome)
        ad[iSample] = convert(Int32, dPopVal(Gsize, aGenome))
    end

    return median(ad), var(ad)
end


# +
aExpectedDiversityPop    = Vector{Vector{Float64}}(undef, length(aJobID))
aVarianceDiversityPop    = Vector{Vector{Float64}}(undef, length(aJobID))

for iID in eachindex(aJobID)
    aExpectedDiversityPop[iID]    = Vector{Float64}(undef, length(aaMutFactor[iID]))
    aVarianceDiversityPop[iID]    = Vector{Float64}(undef, length(aaMutFactor[iID]))

    for (iμ, μ) in enumerate(aaMutFactor[iID])
        aExpectedDiversityPop[iID][iμ], aVarianceDiversityPop[iID][iμ] =
            dPopStat( DIMGSPACE, aTraj[iID][iβ,iμ].aPopCmp )
    end
end

# +
# title("Average Population Diversity");
for (iID, dty) in enumerate(aExpectedDiversityPop[end:end])
    errorbar( 4aaMutFactor[end][rngPlt], dty[rngPlt], #sqrt.(aVarianceDiversityPop[end][rngPlt]),
        fmt="o", fillstyle="none" )
end
axvline(x=1/NPOP, linestyle="dotted", color="tab:grey")
xscale("log")
xlabel("variation probability, μ")
ylabel("median of no. of distinct types, median{d}")
ylim([0,14])

savefig("dty" * string(NPOP) * ".pdf", bbox_inches="tight")
;
