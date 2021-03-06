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
#     display_name: Julia 1.6.1
#     language: julia
#     name: julia-1.6
# ---

versioninfo()

# +
using Revise
using JLD, HDF5
using PyPlot, Printf
using Statistics
# using ProgressMeter, BenchmarkTools
# using Dierckx, LsqFit, ForwardDiff
# using DataFrames, GLM

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

@time for (i, ??) in enumerate(aSelStrength), (j, ??) in enumerate(aMutFactor)
    grid = mGraphs.EdgeWeightedSquareLattice(GRIDSIZE, [ ?? for i in 1:2GRIDSIZE*(GRIDSIZE-1) ]);
    ety = EvolutionaryDynamics.TabularEvotype(
        EvolutionaryDynamics.WithoutReplication(),
        EvolutionaryDynamics.StandardMutation(),
        EvolutionaryDynamics.FitnessSelection(),
        grid
    )

#     traj = @distributed (+) for i in 1:NSAMPLES
    traj = EvolutionaryDynamics.generateTabularSystemsPopulationTrajectories(
            ety=ety, fitnessTbl=fTbl, selCoef=??, Npop=NPOP, nGenRelax=NGENRELAX, nSamples=1 #NSAMPLESPERTRAJ
        )
#     end
    aTraj[1][i,j] = traj
end

# jldopen(folderName * aJobID[1] * "_aTraj.jld", "w") do file
#     addrequire(file, EvolutionaryDynamics)
#     write(file, "aTraj", aTraj[1])
# end
# -

# ## Data Retrival

# +
# sensitivity analysis results

# PopulationTrajectoryData Array: jobID ??? (??, ??) parameter values
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

# +
# using DelimitedFiles

# writing delimitedFiles
# for k in eachindex(aMutFactor)
#     open("rbsPopTrajData/rbsGridEvo_aPopCmp_" * string(k) * ".dat", "w") do io
#         writedlm(io, [ aTraj[end][k].aPopCmp[j][i]
#                 for i in eachindex(aTraj[end][k].aPopCmp[1]),
#                     j in eachindex(aTraj[end][k].aPopCmp) ]
#         )
#     end;
# end

# reading delimitedFiles
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
i?? = 1
i?? = 1

subplots(2,1,figsize=(7.5,7.5))
subplot(211); title("Average Scaled Growth Coeff. | ?? = $(aSelStrength[i??]), ?? = $(aMutFactor[i??])");
    plot( aTraj[iID][i??,i??].avePerformance );
subplot(212); title("Mutation Factor: ???mf??? = $(mean(aTraj[iID][i??,i??].mutationFactor))");
    plot( aTraj[iID][i??,i??].mutationFactor );
#     yscale("log")
tight_layout();
# -

# ## Sensitivity

# +
using LinearAlgebra

i?? = 1
?? = aSelStrength[i??]

# vector of fitness values and matrix of fitness changes
vf = [ exp( ?? * fTbl[g] ) for g in 1:DIMGSPACE ]
m??f = [ vf[go] - vf[gt] for gt in 1:DIMGSPACE, go in 1:DIMGSPACE ]

# construction of the mutation probability matrix
fnGrid(??) = mGraphs.EdgeWeightedSquareLattice(GRIDSIZE, [ ?? for i in 1:2GRIDSIZE*(GRIDSIZE-1) ])
mMut(??) = ( Array(mGraphs.transitionMatrix(fnGrid(??))) .* ( 1.0 .- Diagonal( ones(Float64, DIMGSPACE) ) ) ) / 4??

# vectors of fitness values after mutations
vf??(??) = mMut(??)' * vf
v??f(??) = sum( m??f .* mMut(??), dims=1 )[:]

# m?? = [ (vf[go] - vf[gt])/vf[go] for gt in 1:DIMGSPACE, go in 1:DIMGSPACE ]
# m??Log = [ log(vf[go] / vf[gt]) for gt in 1:DIMGSPACE, go in 1:DIMGSPACE ]

# genotype expected loss of growth
# v??(??)    = sum( m?? .* mMut(??), dims=1 )
# v??Log(??) = sum( m??Log .* mMut(??), dims=1 )

# +
# sensitivity functions
function ??PopVal( vf::AbstractVector, v??f::AbstractVector, aGenome )
    ??F, F = 0.0, 0.0

    for genome in aGenome
        ??F += v??f[genome]
        F += vf[genome]
    end

    return ??F / F
end

function ??PopStat( vf::AbstractVector, v??f::AbstractVector,  aaGenome )
    a?? = Vector{Float64}(undef, length(aaGenome))

    for (iSample, aGenome) in enumerate(aaGenome)
        a??[iSample] = ??PopVal( vf, v??f, aGenome )
    end

    return mean(a??), var(a??)
end

# log sensitivity functions
function ??LogPopVal( vf::AbstractVector, vf??::AbstractVector, aGenome )
    F??, F = 0.0, 0.0

    for genome in aGenome
        F?? += vf??[genome]
        F += vf[genome]
    end

    return log(F/F??)
#     return 1.0 - F??/F
#     return F/F??
end

function ??LogPopStat( vf::AbstractVector, vf??::AbstractVector,  aaGenome )
    a?? = Vector{Float64}(undef, length(aaGenome))

    for (iSample, aGenome) in enumerate(aaGenome)
        a??[iSample] = ??LogPopVal( vf, vf??, aGenome )
    end

    return mean(a??), var(a??)
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

    for (i??, ??) in enumerate(aaMutFactor[iID])
        aExpectedSensitivityPop[iID][i??], aVarianceSensitivityPop[iID][i??] =
            ??PopStat( vf, v??f(??), aTraj[iID][i??,i??].aPopCmp )
        aExpectedSensitivityPopLog[iID][i??], aVarianceSensitivityPopLog[iID][i??] =
            ??LogPopStat( vf, vf??(??), aTraj[iID][i??,i??].aPopCmp )
    end

    for (i??, ??) in enumerate(aaMutFactor[iID][1:end-1])
        aChangeSensitivityPop[iID][i??] =
            ( aExpectedSensitivityPop[iID][i??+1] - aExpectedSensitivityPop[iID][i??] ) /
            ( - aaMutFactor[iID][i??+1] + ?? ) / NPOP
        aChangeSensitivityPopLog[iID][i??] =
            ( aExpectedSensitivityPopLog[iID][i??+1] - aExpectedSensitivityPopLog[iID][i??] ) /
            ( - aaMutFactor[iID][i??+1] + ?? ) / NPOP

        aDivergenceSensitivityPop[iID][i??] =
            ( aChangeSensitivityPop[iID][i??] / aVarianceSensitivityPop[iID][i??] )
        aDivergenceSensitivityPopLog[iID][i??] =
            ( aChangeSensitivityPopLog[iID][i??] / aVarianceSensitivityPopLog[iID][i??] )
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

for (i??, ??) in enumerate(aaMutFactor[end])
    aExpectedSensitivityPop[end][i??], aVarianceSensitivityPop[end][i??] =
        ??PopStat( vf, v??f(??), aTraj[end][i??,i??].aPopCmp )
    aExpectedSensitivityPopLog[end][i??], aVarianceSensitivityPopLog[end][i??] =
        ??LogPopStat( vf, vf??(??), aTraj[end][i??,i??].aPopCmp )
end

for (i??, ??) in enumerate(aaMutFactor[end][1:end-1])
    aChangeSensitivityPop[end][i??] =
        ( aExpectedSensitivityPop[end][i??+1] - aExpectedSensitivityPop[end][i??] ) /
        ( - aaMutFactor[end][i??+1] + ?? ) / NPOP
    aChangeSensitivityPopLog[end][i??] =
        ( aExpectedSensitivityPopLog[end][i??+1] - aExpectedSensitivityPopLog[end][i??] ) /
        ( - aaMutFactor[end][i??+1] + ?? ) / NPOP

    aDivergenceSensitivityPop[end][i??] =
        ( aChangeSensitivityPop[end][i??] / aVarianceSensitivityPop[end][i??] )
    aDivergenceSensitivityPopLog[end][i??] =
        ( aChangeSensitivityPopLog[end][i??] / aVarianceSensitivityPopLog[end][i??] )
end
# -

# #### Interpolation Approach

using Dierckx

# +
i??Min, i??Max, i???? = 2, length(aaMutFactor[end]), 1 # <--- 2, 0, 1
rng, rngDev, rngDevMu = i??Min:i????:i??Max, i??Min:i????:i??Max-1, i??Min:i????:i??Max-1

# neutral weights
wVec = ones(Float64,length(aaMutFactor[end][rng])) # <---
# more weight for low ??
# wVec = [ 1.0 - 0.1 * ( x - rng[1] ) / ( rng[end] - rng[1] ) for x in rng ]
# more weight for high ??
# wVec = [ 2.0( x - rng[1] ) / ( rng[end] - rng[1] ) for x in rng ] 
# more weight for middle-valued ??
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
i??Min, i??Max, i???? = 2, length(aaMutFactor[end]), 1
rng, rngDev, rngDevMu = i??Min:i????:i??Max, i??Min:i????:i??Max-1, i??Min:i????:i??Max-1

# neutral weights
wVec = ones(Float64,length(aaMutFactor[end][rng]))
# more weight for low ??
# wVec = [ 1.0 - 0.1 * ( x - rng[1] ) / ( rng[end] - rng[1] ) for x in rng ]
# more weight for high ??
# wVec = [ 2.0( x - rng[1] ) / ( rng[end] - rng[1] ) for x in rng ] 
# more weight for middle-valued ??
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
i??Min, i??Max, i???? = 1, length(aaMutFactor[end]), 1 # length(aaMutFactor[end])
rng, rngDev, rngDevMu = i??Min:i????:i??Max, i??Min:i????:i??Max-1, i??Min:i????:i??Max-1

??Vals = aaMutFactor[end][rng]
??Vals = aExpectedSensitivityPop[end][rng]
var??Vals = aVarianceSensitivityPop[end][rng]

wVec = (var??Vals).^(-1/2)
wVec[9:12] .= 100.0
# wVec = [ 1.0 - 0.1 * ( x - rng[1] ) / ( rng[end] - rng[1] ) for x in rng ]
# wVec = ones(Float64,length(aaMutFactor[end][rng]))
# wVec = zeros(Float64,length(aaMutFactor[end][rng]))

@. polyRatio3Model(??, p) = ( p[1] + p[2] / ?? + p[3] / ??^2 ) / ( p[4] + p[5] / ?? + p[6] / ??^2 )
@. polyRatioModel(??, p) = ( p[1] + p[2] / ?? ) / ( p[3] + p[4] / ?? )

@. tanhModel(??, p) = p[1] - p[2] * tanh( p[3] + p[4] * ?? )
@. arctanModel(??, p) = p[1] - p[2] * atan( p[3] + p[4] * ?? )
@. gdModel(??, p) = p[1] - p[2] * atan( tanh( p[3] + p[4] * ?? ) )
    algebraicSigmoid(x) = x / sqrt( 1 + x^2 )
@. algebraicModel(??, p) = p[1] - p[2] * algebraicSigmoid( p[3] * ?? )

prmPolyRatio3 = [10.0, 0.3/10.0^2, 0.1, 10.0, 10.0^-2, 0.1]
prmPolyRatio  = [10.0, 0.3/10.0^2, 10.0, 10.0^2]

prmTanh   = [ 0.15, 0.30, 0.01, 0.1 ]
prmArctan = [ 0.15, 0.30, 0.01, 0.1 ]

# fitAveSnt = LsqFit.curve_fit(polyRatioModel, ??Vals, ??Vals, wVec, prmPolyRatio)
# fit??(??) = polyRatioModel(??, coef(fitAveSnt))

# fitAveSnt = LsqFit.curve_fit(tanhModel, ??Vals, ??Vals, wVec, prmTanh)
# fit??(??) = tanhModel(??, coef(fitAveSnt))

# fitAveSnt = LsqFit.curve_fit(arctanModel, ??Vals, ??Vals, wVec, prmTanh)
# fit??(??) = arctanModel(??, coef(fitAveSnt))

fitAveSnt = LsqFit.curve_fit(gdModel, ??Vals, ??Vals, wVec, prmTanh)
fit??(??) = gdModel(??, coef(fitAveSnt))

# fitAveSnt = LsqFit.curve_fit(algebraicModel, ??Vals, ??Vals, wVec, prmTanh)
# fit??(??) = algebraicModel(??, coef(fitAveSnt))


prmArctanVar = [ 0.04, 0.04, 0.01, 0.1 ]

# fitVarSnt = LsqFit.curve_fit(sigmoidModel2, ??Vals, var??Vals, p2)
# fitVar??(??) = sigmoidModel2(??, coef(fitVarSnt))

# fitVarSnt = LsqFit.curve_fit(arctanModel, ??Vals, var??Vals, prmArctanVar)
# fitVar??(??) = arctanModel(??, coef(fitVarSnt))


devSnt = ( - [ ForwardDiff.derivative(fit??, ??) for ?? in ??Vals ] / NPOP ) ./ var??Vals
# devSntFit = ( - [ ForwardDiff.derivative(fit??, ??) for ?? in ??Vals ] / NPOP ) ./ fitVar??.(??Vals)

subplots(1, 3, figsize=(15.0, 4.0))
subplot(131);
    errorbar( 4aaMutFactor[end], aExpectedSensitivityPop[end], fmt="o", fillstyle="none" )
    plot( 4??Vals, fit??.(??Vals), "--");
    xscale("log")
subplot(132);
    errorbar( aaMutFactor[end], aVarianceSensitivityPop[end], fmt="o", fillstyle="none" );
#     plot( ??Vals, fitVar??.(??Vals), "--" );
    xscale("log")
subplot(133); axhline(y=1, linestyle="--");
    errorbar( 4aaMutFactor[end][rngDevMu], aDivergenceSensitivityPop[end][rngDev], fmt="o", fillstyle="none" )
#     plot( ??Vals, devSntFit, ".");
    plot( 4??Vals, devSnt, "+");
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
rngPltItp = 1:i??Max-i??Min

# subplots(3, 2, figsize=(10.0, 10.0))
subplots(2, 2, figsize=(10.0, 7.0))


subplot(221); title("(a) Average Population Sensitivity");
for (iID, E??) in enumerate(aExpectedSensitivityPop[end:end])
    errorbar( 4aaMutFactor[end][rngPlt], E??[rngPlt], sqrt.(aVarianceSensitivityPop[end][rngPltDev]),
        fmt="o", fillstyle="none" )
end
# plot( 4aaMutFactor[end][rngPltFit], fit??.(aaMutFactor[end][rngPltFit]), color="tab:red", linestyle="--" )
plot( 4aaMutFactor[end][rngPltFit], itpAveSnt(aaMutFactor[end][rngPltFit]), color="tab:red", linestyle="--" )
xscale("log")
xlabel("variation probability, ??")
ylabel("Average Sensitivity, ???????????")

subplot(224); title("(d) Average Population Log Sensitivity");
for (iID, E??) in enumerate(aExpectedSensitivityPopLog[end:end])
    errorbar( 4aaMutFactor[end][rngPlt], E??[rngPlt], sqrt.(aVarianceSensitivityPopLog[end][rngPlt]),
        fmt="o", fillstyle="none" )
end
xlabel("variation probability, ??")
ylabel("Average Log Sensitivity, ?????lnF??????")
xscale("log")


subplot(222); title("(b) Variance of Population Sensitivity");
for (iID, c) in enumerate(aVarianceSensitivityPop[end:end])
    errorbar( 4aaMutFactor[end][rngPlt], c[rngPlt], fmt="o", fillstyle="none" )
end
# plot( aaMutFactor[end][rngPltDev], fitVar??.(aaMutFactor[end][rngPltDev]), color="tab:red", linestyle="--" )
xscale("log")
yscale("log")
xlabel("variation probability, ??")
ylabel("Variance Sensitivity, Var{?????}")
# 
# subplot(324); title("(d) Variance Population Log Sensitivity");
# for (iID, c) in enumerate(aVarianceSensitivityPopLog[end:end])
#     errorbar( aaMutFactor[end][rng], c[rng], fmt="o", fillstyle="none" )
# end
# xscale("log")
# yscale("log")
# ylabel("Variance Log Sensitivity, Var{??lnF???}")


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
xlabel("variation probability, ??")
ylabel("deviation, ??")
ylim([-1,6]);

# subplot(222); title("(f) Deviations from Fluctuations Relation");
# for (iID, c) in enumerate(aDivergenceSensitivityPopLog[end:end])
#     errorbar( 4aaMutFactor[end][1:end-1], c, fmt="o", fillstyle="none" )
# end
# axhline(y=1, linestyle="--", color="tab:grey")
# xscale("log")
# xlabel("variation probability, ????")
# ylabel("deviation, ??")

# suptitle("Population Sensitivity Statistics: N = 60", fontsize=16)
tight_layout();

# savefig("snt" * string(NPOP) * ".pdf", bbox_inches="tight")
# +
fit2clr = Dict( 1.0 => "tab:blue", 5.0 => "tab:red" )
# fit2clr = Dict( 1.0 => "tab:grey", 5.0 => "tab:orange" )
fit2shp = Dict( 1.0 => "s", 5.0 => "s" )
fit2siz = Dict( 1.0 => 130, 5.0 => 230 )

function squareGridPlot(L; l=1, fsize=7, s0=210, colormap="cividis")
    V = [ [ (i-1)??L + 1, (i-1)%L + 1 ] for i in 1:L^2 ]

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
rngPltItp = 1:i??Max-i??Min

# subplots(3, 2, figsize=(10.0, 10.0))
fig, ax = subplots(1, 1, figsize=(5.0, 3.5))

# subplot(111); title("Average Population Sensitivity");
for (iID, E??) in enumerate(aExpectedSensitivityPop[end:end])
    ax.errorbar( 4aaMutFactor[end][rngPlt], E??[rngPlt], sqrt.(aVarianceSensitivityPop[end][rngPltDev]),
        fmt="o", fillstyle="none" )
end
# plot( 4aaMutFactor[end][rngPltFit], fit??.(aaMutFactor[end][rngPltFit]), color="tab:red", linestyle="--" )
# plot( 4aaMutFactor[end][rngPltFit], itpAveSnt(aaMutFactor[end][rngPltFit]), color="tab:red", linestyle="--" )
xscale("log")
xlabel("variation probability, ??")
ylabel("Average Sensitivity, ???????????")

axin = ax.inset_axes([0.62, 0.62, 0.35, 0.33])
# for (iID, c) in enumerate(aDivergenceSensitivityPop[end:end])
#     errorbar( 4aaMutFactor[end][rngPltDev], c[rngPltDev], fmt="o", fillstyle="none" )
# end
# # plot( aaMutFactor[end][rngDev], predict(lmDiv), color="tab:blue", linestyle="-." )
axin.plot( 4aaMutFactor[end][rngPltDev], devSnt[rngPltItp], color="tab:red", "*" )
# # plot( aaMutFactor[end][rng], predict(lmDivItp), color="tab:red", linestyle="-." )
axin.axhline(y=1, linestyle="--", color="tab:grey")
# axvline(x=1/NPOP, ymax=1.0, linestyle="dotted", color="tab:grey")
# text(1/(NPOP-8), 0.2, "??? 1/N")
# # axvline(x=10/NPOP, linestyle="--", color="tab:grey")
axin.semilogx()
axin.set_xlabel("??", fontsize=8, labelpad=-3)
axin.set_ylabel("deviation, ??", fontsize=8, labelpad=0)
axin.set_ylim([0,4.5]);
axin.tick_params(axis="x", labelsize=8)
axin.tick_params(axis="y", labelsize=8)

tight_layout();
# savefig("snt" * string(NPOP) * "mt.pdf", bbox_inches="tight")
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

    for (i??, ??) in enumerate(aaMutFactor[iID])
        aExpectedSensitivityPop[iID][i??], aVarianceSensitivityPop[iID][i??] =
            ??PopStat( vf, v??f(??), aTraj[iID][i??,i??].aPopCmp )
        aExpectedSensitivityPopLog[iID][i??], aVarianceSensitivityPopLog[iID][i??] =
            ??LogPopStat( vf, vf??(??), aTraj[iID][i??,i??].aPopCmp )
    end

    for (i??, ??) in enumerate(aaMutFactor[iID][1:end-1])
        aChangeSensitivityPop[iID][i??] =
            ( aExpectedSensitivityPop[iID][i??+1] - aExpectedSensitivityPop[iID][i??] ) /
            ( - aaMutFactor[iID][i??+1] + ?? ) / aNpop[iID]
        aChangeSensitivityPopLog[iID][i??] =
            ( aExpectedSensitivityPopLog[iID][i??+1] - aExpectedSensitivityPopLog[iID][i??] ) /
            ( - aaMutFactor[iID][i??+1] + ?? ) / aNpop[iID]

        aDivergenceSensitivityPop[iID][i??] =
            ( aChangeSensitivityPop[iID][i??] / aVarianceSensitivityPop[iID][i??] )
        aDivergenceSensitivityPopLog[iID][i??] =
            ( aChangeSensitivityPopLog[iID][i??] / aVarianceSensitivityPopLog[iID][i??] )
    end
end

# +
subplots(3, 2, figsize=(10.0, 10.0))

subplot(321); title("(a) Average Population Sensitivity");
    for (iID, E??) in enumerate(aExpectedSensitivityPop)
        errorbar( aaMutFactor[iID], E??, sqrt.(aVarianceSensitivityPop[iID]), fmt="o", fillstyle="none" ) #,
#             color="tab:blue" )
    end
#     axhline(y=0, linestyle="--")
    xscale("log")
    ylabel("Average Sensitivity, ???????????")

subplot(322); title("(b) Average Population Log Sensitivity");
    for (iID, E??) in enumerate(aExpectedSensitivityPopLog)
        errorbar( aaMutFactor[iID], E??, sqrt.(aVarianceSensitivityPop[iID]), fmt="o", fillstyle="none" )
    end
#     axhline(y=1, linestyle="--")
    ylabel("Average Log Sensitivity, ?????lnF??????")
    xscale("log")

subplot(323); title("(c) Variance of Population Sensitivity");
    for (iID, c) in enumerate(aVarianceSensitivityPop)
        errorbar( aaMutFactor[iID], c, fmt="o", fillstyle="none" )
    end
    xscale("log")
    ylabel("Variance Sensitivity, Var{?????}")

subplot(324); title("(d) Variance Population Log Sensitivity");
    for (iID, c) in enumerate(aVarianceSensitivityPopLog)
        errorbar( aaMutFactor[iID], c, fmt="o", fillstyle="none" )
    end
    xscale("log")
    ylabel("Variance Log Sensitivity, Var{??lnF???}")

# ylims = [0.01, 10^4]
subplot(325); title("(e) Deviations from Fluctuations Relation");
    for (iID, c) in enumerate(aDivergenceSensitivityPop)
#     for (iID, c) in enumerate(aChangeSensitivityPopNeutral)
        errorbar( aaMutFactor[iID][1:end-1], c, fmt="o", fillstyle="none" )
    end
    axhline(y=1, linestyle="--")
    xscale("log")
    xlabel("variation probability, ????")
    ylabel("deviation, ??")
#     yscale("log")
#     ylim(ylims)

subplot(326); title("(f) Deviations from Fluctuations Relation");
    for (iID, c) in enumerate(aDivergenceSensitivityPopLog)
    #     for (iID, c) in enumerate(aChangeSensitivityPop)
        errorbar( aaMutFactor[iID][1:end-1], c, fmt="o", fillstyle="none" )
    end
    axhline(y=1, linestyle="--")
    xscale("log")
    xlabel("variation probability, ????")
#     ylabel("deviation, ??")
#     yscale("log")
#     ylim(ylims)

# suptitle("Population Sensitivity Statistics: N = 60", fontsize=16)
tight_layout();

savefig("snt100.pdf", bbox_inches="tight")
# -
aExpectedSensitivityPop[4] - aExpectedSensitivityPop[2]

# +
# probability vectors of n-th order genotypic dynamics
m??(??, n=1) = sum( ( af .* mMut(??) )^n, dims=1 )
Z(??, n=1) = m??(??, n) * af

prob(??, n=1) = af' .* m??(??, n) ./ Z(??, n)
probNeutral(??, n=1) = af' .* m??(??, n) .* a??' / sum(af' .* m??(??, n) .* a??')

function probInf(??)
    r = real( eigvecs( af .* mMut(??) )[:,end] )
    return r' ./ sum(r)
end

function probInfNeutral(??)
    r = real( eigvecs( af' .* mMut(??) )[:,end] )
    return ( r .* a?? )' ./ sum(r .* a??)
end

# expected sensitivity evaluations
Esensitivity(??, n=1) = sum( mSensitivity .* mMut(??) .* prob(??, n) ) / 4??
EsensitivityNeutral(??, n=1) = sum( mSensitivity .* mMut(??) .* probNeutral(??, n) ) / 4??

EsensitivityInf(??) = sum( mSensitivity .* mMut(??) .* probInf(??) ) / 4??
EsensitivityInfNeutral(??) = sum( mSensitivity .* mMut(??) .* probInfNeutral(??) ) / 4??

# average fitness values
aveF(??, n=1) = prob(??, n) ??? af
aveFInf(??) = probInf(??) ??? af
;

# +
# sensitivity vectors of theoretical models

an = [0, 1, 100]
e??Iterator = -5.0:0.05:-1.0
a?? = [ 2.5 * 10.0^e?? for e?? in e??Iterator ]

aSens = Vector{Vector{Float64}}(undef, length(an))
aSensNeutral = Vector{Vector{Float64}}(undef, length(an))
aAveF = Vector{Vector{Float64}}(undef, length(an))

for (i, n) in enumerate(an)
    aSens[i] = [ Esensitivity(??, n) for ?? in a?? ]
    aSensNeutral[i] = [ EsensitivityNeutral(??, n) for ?? in a?? ]
    aAveF[i] = [ aveF(??, n) for ?? in a?? ]
end

aSensInf = [ EsensitivityInf(??) for ?? in a?? ]
aSensInfNeutral = [ EsensitivityInfNeutral(??) for ?? in a?? ]
aAveFInf = [ aveFInf(??) for ?? in a?? ]
;
# -

# ### Expected Growth Rate and Diversity

# +
function fPopVal( vf::AbstractVector, aGenome )
    F = 0.0

    for genome in aGenome
        F += vf[genome]
    end

    return F
end

function fPopStat( vf::AbstractVector,  aaGenome )
    af = Vector{Float64}(undef, length(aaGenome))

    for (iSample, aGenome) in enumerate(aaGenome)
        af[iSample] = log( fPopVal( vf, aGenome ) ) / log( maximum(vf) * length(aGenome) )
    end

    return mean(af), var(af)
end

# +
function dPopVal( Gsize::Integer, aGenome )
    ??pop = zeros(Int8, Gsize)

    for genome in aGenome
        ??pop[genome] = Int8(1)
    end

    return sum(??pop)
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

aExpectedReproductionPop    = Vector{Vector{Float64}}(undef, length(aJobID))
aVarianceReproductionPop    = Vector{Vector{Float64}}(undef, length(aJobID))

for iID in eachindex(aJobID)
    aExpectedDiversityPop[iID]    = Vector{Float64}(undef, length(aaMutFactor[iID]))
    aVarianceDiversityPop[iID]    = Vector{Float64}(undef, length(aaMutFactor[iID]))

    aExpectedReproductionPop[iID]    = Vector{Float64}(undef, length(aaMutFactor[iID]))
    aVarianceReproductionPop[iID]    = Vector{Float64}(undef, length(aaMutFactor[iID]))

    for (i??, ??) in enumerate(aaMutFactor[iID])
        aExpectedDiversityPop[iID][i??], aVarianceDiversityPop[iID][i??] =
            dPopStat( DIMGSPACE, aTraj[iID][i??,i??].aPopCmp )
        
        aExpectedReproductionPop[iID][i??], aVarianceReproductionPop[iID][i??] =
            fPopStat( vf, aTraj[iID][i??,i??].aPopCmp )
    end
end
# -

aExpectedReproductionPop[end:end][1];

# +
using LinearAlgebra

i?? = 1
?? = aSelStrength[i??]

# vector of fitness values and matrix of fitness changes
vf = [ exp( ?? * fTbl[g] ) for g in 1:DIMGSPACE ]
mSnt = [ ( vf[go] - vf[gt] ) / vf[go] for gt in 1:DIMGSPACE, go in 1:DIMGSPACE ]

# construction of the mutation probability matrix
fnGrid(??) = mGraphs.EdgeWeightedSquareLattice(GRIDSIZE, [ ?? for i in 1:2GRIDSIZE*(GRIDSIZE-1) ])
????(??) = ( Array(mGraphs.transitionMatrix(fnGrid(??))) .* ( 1.0 .- Diagonal( ones(Float64, DIMGSPACE) ) ) )

vSnt(??) = sum( mSnt .* ????(??), dims=1 )

# vectors of fitness values after mutations
vf????(??) = ????(??)' * vf

subplots(1, 2, figsize=(12.0, 4.5))

subplot(122); title("(b) Population Reproduction Rate");
for (iID, dty) in enumerate(aExpectedReproductionPop[end:end])
    errorbar( 4aaMutFactor[end][rngPlt], (dty[rngPlt]), #sqrt.(aVarianceDiversityPop[end][rngPlt]),
        fmt="o", fillstyle="none" )
end
axvline(x=1/NPOP, linestyle="dotted", color="tab:grey")
xscale("log")
xlabel("variation probability, ??")
# yscale("log")
ylim([0.6,1.1])
ylabel("rescaled log cumulative growth, ???lnF??????/ln(Nf*)")

subplot(121); title("(a) Population Diversity");
for (iID, dty) in enumerate(aExpectedDiversityPop[end:end])
    errorbar( 4aaMutFactor[end][rngPlt], dty[rngPlt], #sqrt.(aVarianceDiversityPop[end][rngPlt]),
        fmt="o", fillstyle="none" )
end
axvline(x=1/NPOP, linestyle="dotted", color="tab:grey")
xscale("log")
xlabel("variation probability, ??")
ylabel("median of no. of distinct types, median{d}")
ylim([0,14])

savefig("dty" * string(NPOP) * ".pdf", bbox_inches="tight")
# savefig("further" * string(NPOP) * ".pdf", bbox_inches="tight")
;
# -

