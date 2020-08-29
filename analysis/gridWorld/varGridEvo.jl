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
using Revise, BenchmarkTools, PyPlot, Printf
import EvolutionaryDynamics, mUtils, mPlot, mGraphs

# Base.show(io::IO, f::Float64) = @printf(io, "%.2f", f)

# +
folderName = "/Users/riccardorao/projects/fisheria/analysis/gridWorld/scratch/var/"

jobID = "varGridEvo-7785617"

include(folderName * jobID * "-parameters.jl")

# +
fMat = mGraphs.matrixForm( mGraphs.VertexWeightedSquareLattice( GRIDSIZE, aFitnessTbl[1] ) )

# imshow(fMat, cmap="cividis"); # origin="lower"
imshow(fMat, cmap="inferno"); # origin="lower"

xticks( prepend!(collect(3:4:GRIDSIZE-1),0), prepend!(collect(4:4:GRIDSIZE),1) )
yticks( prepend!(collect(3:4:GRIDSIZE-1),0), prepend!(collect(4:4:GRIDSIZE),1) )
# colorbar()

savefig("varLnF1.pdf", bbox_inches="tight")
;
# -

# ### Data Generation/Retrival

# +
aTraj = Array{EvolutionaryDynamics.VaryingEnvironmentTrajectoryData}(undef,
    length(aSelStrength), length(aMutFactor), length(aEnvTransRate))

@time for (i,β) in enumerate(aSelStrength), (j,μ) in enumerate(aMutFactor), (k,w) in enumerate(aEnvTransRate)
    
    grid = mGraphs.EdgeWeightedSquareLattice(GRIDSIZE, [ μ for i in 1:2GRIDSIZE*(GRIDSIZE-1) ]);
    ety = EvolutionaryDynamics.TabularEvotype(
#         EvolutionaryDynamics.NeutralReplication(REPCOEF),
        EvolutionaryDynamics.WithoutReplication(),
        EvolutionaryDynamics.StandardMutation(),
        EvolutionaryDynamics.FitnessSelection(),
        grid
    )
    transMtx = [ [ 1.0 - w, w ] [ w, 1.0 - w ] ]
    
#     traj = @distributed (+) for i in 1:NSAMPLES
    traj = EvolutionaryDynamics.generateTabularSystemsTrajectories(
            ety=ety, aFitnessTbl=aFitnessTbl, selCoef=β, transMtx=transMtx, Npop=NPOP,
            NgenRelax=convert(Int32,RELNGENRELAX/w), NgenSample=convert(Int64,RELNGENSAMPLE/w)
        )
#     end
    aTraj[i,j,k] = traj
end
# -

using JLD, HDF5, DelimitedFiles
aTraj = mUtils.readJLD(folderName * jobID * "_aTraj.jld", "aTraj");

# ### Analysis

# +
i = 1

subplots(3,1,figsize=(7.5,7.5))
subplot(311); title("Average Fitness");
    plot( aTraj[i].avePerformance );
subplot(312); title("Environmental State");
    plot( aTraj[i].envState );
subplot(313); title("Mutation Factor");
    plot( aTraj[i].mutationFactor );
tight_layout();

# +
Base.show(io::IO, f::Float64) = @printf(io, "%.2f", f)

figNum = figure(figsize=(12,12))

iβ = 3
β = aSelStrength[iβ]

# for (iμ,μ) in enumerate(aMutFactor), (iw,w) in enumerate(aEnvTransRate)
# for (iμ,μ) in zip(2:3,aMutFactor[2:3]), (iw,w) in zip(1:2,aEnvTransRate[1:2])
for (iμ,μ) in zip(2:2,aMutFactor[2:2]), (iw,w) in zip(2:2,aEnvTransRate[2:2])

    grid = mGraphs.EdgeWeightedSquareLattice(GRIDSIZE, [ μ for i in 1:2GRIDSIZE*(GRIDSIZE-1) ]);
    Pgty = [ sum( aTraj[iβ,iμ,iw].jointProb[g,:,:,:] ) for g in 1:size(aTraj[iβ,iμ,iw].jointProb)[1] ] /
        aTraj[iβ,iμ,iw].NgenSample

    p = mGraphs.VertexWeightedSquareLattice( GRIDSIZE, Pgty )
    pMat = mGraphs.matrixForm(p);

    ax = figNum.add_subplot(3,3,iμ + 3*(3-iw))

    imshow(log.(pMat), cmap="inferno", vmax=-3, vmin=-12);
#     title("β = $(β), μ = $(μ), τ = $(1.0/w)");

#     eMin = -13
#     for (i,e) in enumerate(pNicMat)
#         if e == 0
#            pNicMat[i] = 10.0^eMin
#         end
#     end

    xticks([])
    yticks([])

    if iw == 1 xticks( prepend!(collect(5:6:GRIDSIZE-1),0), prepend!(collect(6:6:GRIDSIZE),1) ) end
    if iμ == 2 yticks( prepend!(collect(5:6:GRIDSIZE-1),0), prepend!(collect(6:6:GRIDSIZE),1) ) end
    if iw == 1 && iμ == 3 colorbar(fraction=0.045) end
end

# savefig("varPmfPhaseDigramRaw.pdf", bbox_inches="tight")
savefig("evolvability.pdf", bbox_inches="tight")
# -
# ### Relaxation Time-scale Anaysis

# +
using LinearAlgebra

iβ = 3
β = aSelStrength[iβ]

# array of fitness values
af = [ exp( β * aFitnessTbl[1][g] ) for g in 1:DIMGSPACE ]

# construction of the transition matrix for genotype dynamics as function of the mutation probability
gridu(μ) = mGraphs.EdgeWeightedSquareLattice(GRIDSIZE, [ μ for i in 1:2GRIDSIZE*(GRIDSIZE-1) ])
mMut(μ) = Array(mGraphs.transitionMatrix(gridu(μ))) .+ Diagonal([ 1.0 for i in 1:DIMGSPACE ])
mW(μ) = af .* mMut(μ) ./ sum(af .* mMut(μ), dims=1)

# relaxation time scales as function of the mutation probability
τ1(μ) = - 1.0 / log( eigvals( mW(μ) )[end-1] )
τ2(μ) = - 1.0 / log( eigvals( mW(μ) )[end-2] )
τ3(μ) = - 1.0 / log( eigvals( mW(μ) )[end-3] )
τ4(μ) = - 1.0 / log( eigvals( mW(μ) )[end-4] ) ;

# +
eμIterator = -2.0:0.1:-1.0

aμ = [ 2.0 * 10.0^eμ for eμ in eμIterator ]
aτ1 = [ τ1(μ) for μ in aμ ]
aτ2 = [ τ2(μ) for μ in aμ ]
aτ3 = [ τ3(μ) for μ in aμ ]
aτ4 = [ τ4(μ) for μ in aμ ] ;

# +
plot( aμ, aτ1 )
plot( aμ, aτ2 )
plot( aμ, aτ3 )
plot( aμ, aτ4 )

xscale("log")
yscale("log")

xlabel("mutation rate, μ")
ylabel("relaxation time-scale, τ")

for (μ, w) in zip(aMutFactor, aEnvTransRate)
    axvline( x=μ, color="black", linewidth=1, linestyle="--" )
    axhline( y=1.0/w, xmin=aμ[1], color="black", linewidth=1, linestyle="--" )
end
# -

# The conclusion that I would draw is the following.
# The second and faster time scales better capture the relaxation time scale of the population.
# When the environmental change time scale is faster than the population relaxation one (horizontal dashed lines below the orange, green, and red lines) then "evolvability" is possible.

# ### Information Theoretical Analysis

# +
Jγ(P) = sum(P, dims=4)
Jε(P) = sum(P, dims=1)

Pγ(P) = sum(P, dims=[2,3,4])
Pε(P) = sum(P, dims=[1,3,4])
Pγε(P) = sum(P, dims=[3,4])

PγPast(P) = sum(P, dims=[1,2,4])
PεPast(P) = sum(P, dims=[1,2,3])
PγεPast(P) = sum(P, dims=[1,2])

F(β) = [ β*(afTb[ie][(iγ-1)%DIMGSPACE + 1] - afTb[ie][(iγPast-1)%DIMGSPACE + 1]) for
        iγ in 1:5DIMGSPACE, ie in 1:NENV, iγPast in 1:5DIMGSPACE, iePast in 1:1 ]

fitnessGain(P,β) = sum( Jγ(P) .* F(β) )

infoGain(P) = sum( Jγ(P) .* mUtils.ReLog.( ( Pγε(P) .* PγPast(P) ) ./ ( sum(P, dims=[1,4]) .* Pγ(P) ) ) )

uselessInfo(P) = - sum( Jε(P) .* mUtils.ReLog.( ( sum(P, dims=[1,4]) .* PεPast(P) ) ./ ( PγεPast(P) .* Pε(P) ) ) )

# function EP(P)
#     EP = 0.0
#     for iγ in 1:size(P)[1], iε in 1:size(P)[2], iγpast in 1:size(P)[3]
#         EP += Jγ(P)[iγ, iε, iγpast, 1] * mUtils.ReLog( Jγ(P)[iγ, iε, iγpast, 1] / Jγ(P)[iγpast, iε, iγ, 1] )
#     end
#     return EP
# end

# +
Base.show(io::IO, f::Float64) = @printf(io, "%.6f", f)

k = 2
β = aSelStrength[k]

IGmat = Matrix{Float64}(undef, 3, 3)
UImat = Matrix{Float64}(undef, 3, 3)
FGmat = Matrix{Float64}(undef, 3, 3)
# EPmat = Matrix{Float64}(undef, 3, 3)

#     aS = [ exp(f[g]*pop.env.aSelCoef[2]) for g in 1:DIMGSPACE ]
#     apS = aS ./ sum(aS)

#     aP = ( ( Array(mGraphs.transitionMatrix(pop.ety.G)) .+ Diagonal([ 1.0 + MINREPCOEF for i in 1:DIMGSPACE ]) ) * aS
#         ) .* aS
#     apP = aP ./ sum(aP)

#     gpS = mGraphs.VertexWeightedSquareLattice( GRIDSIZE, apS )
#     pSMat = mGraphs.matrixForm(gpS);

#     gpP = mGraphs.VertexWeightedSquareLattice( GRIDSIZE, apP )
#     pPMat = mGraphs.matrixForm(gpP);

for (i,λM) in enumerate(aλM), (j,W) in enumerate(aW)
    P = aTraj[i,j,k].jointProb / aTraj[i,j,k].NgenSample

    IGmat[i,j] = infoGain(P)
    UImat[i,j] = uselessInfo(P)
    FGmat[i,j] = fitnessGain(P,β)
#     EPmat[i,j] = EP(P)
end
# -

# j = 3 [:,j]
display( IGmat )
display( UImat )
display( FGmat )
# EPmat
