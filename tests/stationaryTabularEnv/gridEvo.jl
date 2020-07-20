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
using Revise, BenchmarkTools, PyPlot, Printf, JLD, HDF5, DelimitedFiles, ProgressMeter
import EvolutionaryDynamics, mUtils, mPlot, mGraphs

Base.show(io::IO, f::Float64) = @printf(io, "%.2f", f)

# +
folderName = "/home/rjku/projects/fisheria/tests/stationaryTabularEnv/"

jobID = "gridEvo"

include(folderName * jobID * "-parameters.jl")

# aTraj = mUtils.readJLD(folderName * jobID * "_aTraj.jld", "aTraj");

# +
fMat = mGraphs.matrixForm( mGraphs.VertexWeightedSquareLattice( GRIDSIZE, fTbl ) )

imshow(fMat,cmap="cividis"); # origin="lower"

xticks( prepend!(collect(3:4:GRIDSIZE-1),0), prepend!(collect(4:4:GRIDSIZE),1) )
yticks( prepend!(collect(3:4:GRIDSIZE-1),0), prepend!(collect(4:4:GRIDSIZE),1) )
colorbar()

# savefig("lnf.pdf", bbox_inches="tight")
# -

# ### Data Generation

# +
aTraj = Array{EvolutionaryDynamics.TrajectoryData}(undef, length(aSelStrength), length(aMutFactor))

@time for (i,β) in enumerate(aSelStrength), (j,μ) in enumerate(aMutFactor)
    grid = mGraphs.EdgeWeightedSquareLattice(GRIDSIZE, [ μ for i in 1:2GRIDSIZE*(GRIDSIZE-1) ]);
#     traj = @distributed (+) for i in 1:NSAMPLES
    traj = EvolutionaryDynamics.generateTabularSystemsTrajectories(
            REPFACTOR, MINREPCOEF, μ, grid, fTbl, REPSTRENGTH, β, NPOP, NGENRELAX, NGENSAMPLE
        )
#     end
    aTraj[i,j] = traj
end
# -

# ### Analysis

# +
i = 1

subplots(2,1,figsize=(7.5,7.5))
subplot(211); title("Average Fitness");
    plot( aTraj[i].avePerformance );
subplot(212); title("Mutation Factor");
    plot( aTraj[i].mutationFactor );
tight_layout();

# +
Base.show(io::IO, f::Float64) = @printf(io, "%.2f", f)

figNum = figure(figsize=(12,12))

for (i,μ) in enumerate(aMutFactor), (j,β) in enumerate(aSelStrength)

    grid = mGraphs.EdgeWeightedSquareLattice(GRIDSIZE, [ μ for i in 1:2GRIDSIZE*(GRIDSIZE-1) ]);
    Pgty = [ sum( aTraj[i,j].jointProb[g,:] ) for g in 1:size(aTraj[i,j].jointProb)[1] ] /
        aTraj[i,j].NgenSample

    p = mGraphs.VertexWeightedSquareLattice( GRIDSIZE, Pgty )
    pMat = mGraphs.matrixForm(p);

    ax = figNum.add_subplot(3,3,(4-i) + 3*(j-1))

    imshow(log.(pMat),cmap="cividis",vmax=0, vmin=-10);
    title("β = $(β), μ = $(μ)");

#     eMin = -13
#     for (i,e) in enumerate(pNicMat)
#         if e == 0
#            pNicMat[i] = 10.0^eMin
#         end
#     end

    xticks([])
    yticks([])

    if j == 3 xticks( prepend!(collect(3:4:GRIDSIZE-1),0), prepend!(collect(4:4:GRIDSIZE),1) ) end
    if i == 1 yticks( prepend!(collect(3:4:GRIDSIZE-1),0), prepend!(collect(4:4:GRIDSIZE),1) ) end
    if i == 1 && j == 3 colorbar(fraction=0.045) end
end

# savefig("plots/PMFphaseDiagramTheory.pdf", bbox_inches="tight")
# -
# #### Information Theoretical Analysis

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

# #### Control Dynamics Analysis

# +
Base.show(io::IO, f::Float64) = @printf(io, "%.2f", f)

figNum = figure(figsize=(12,12))

i = 1
λM = aλM[i]

k = 3
β = aSelStrength[k]

for (m,μ) in enumerate(aMutFactor), (j,W) in enumerate(aW)

    Pgty = [ sum( aTrajCtrl[i,j,k,m].jointProb[g,:,:,:] ) for g in 1:size(aTrajCtrl[i,j,k,m].jointProb)[1] ] /
        aTrajCtrl[i,j,k,m].NgenSample
    Pgty1 = [ sum( [ Pgty[ g1 + DIMGSPACE * g2 ] for g2 in 0:4 ] ) for g1 in 1:DIMGSPACE ]
    Bgty1 = [ [ Pgty[ g1 + DIMGSPACE * g2 ] for g2 in 1:4 ] / Pgty1[g1] for g1 in 1:DIMGSPACE ];

    p = mGraphs.VertexWeightedSquareLattice( GRIDSIZE, Pgty1 )
    pMat = mGraphs.matrixForm(p);

    ax = figNum.add_subplot(3,3,m + 3*(j-1))

    title("μ = $(μ), w = $(log10(aEnvTransRate[j])), β = $(β)");

    imshow(log.(pMat), cmap="cividis", vmax=0, vmin=-10, origin="lower");

#     eMin = -13
#     for (i,e) in enumerate(pNicMat)
#         if e == 0
#            pNicMat[i] = 10.0^eMin
#         end
#     end

    xticks([])
    yticks([])

    if j == 3
        xticks( prepend!(collect(3:4:GRIDSIZE-1),0), prepend!(collect(4:4:GRIDSIZE),1) )
    end

    if m == 1
        yticks( prepend!(collect(3:4:GRIDSIZE-1),0), prepend!(collect(4:4:GRIDSIZE),1) )
    end

    if m == 3 && j == 3
       colorbar(fraction=0.045)
    end
end

# savefig("plots/PMFphaseDiagramTheory.pdf", bbox_inches="tight")
# +
Base.show(io::IO, f::Float64) = @printf(io, "%.2f", f)

figNum = figure(figsize=(20,20))

i = 1
λM = aλM[i]

k = 3
β = aSelStrength[k]

for (m,μ) in enumerate(aMutFactor), (j,W) in enumerate(aW)

    Pgty = [ sum( aTrajCtrl[i,j,k,m].jointProb[g,:,:,:] ) for g in 1:size(aTrajCtrl[i,j,k,m].jointProb)[1] ] /
        aTrajCtrl[i,j,k,m].NgenSample
    Pgty1 = [ sum( [ Pgty[ g1 + DIMGSPACE * g2 ] for g2 in 0:4 ] ) for g1 in 1:DIMGSPACE ]
    Bgty1 = [ [ Pgty[ g1 + DIMGSPACE * g2 ] for g2 in 1:4 ] / Pgty1[g1] for g1 in 1:DIMGSPACE ];

    p = mGraphs.VertexWeightedSquareLattice( GRIDSIZE, Pgty1 )
    pMat = mGraphs.matrixForm(p);

    ax = figNum.add_subplot(3,3,m + 3*(j-1))

    title("μ = $(μ), w = $(log10(aEnvTransRate[j])), β = $(log10(β))");

    mPlot.squareGridArrowPlot(GRIDSIZE, Pgty1, Bgty1, fsize=10, s0=100 )

    xticks([])
    yticks([])

    if j == 3
        xticks( prepend!(collect(3:4:GRIDSIZE-1),0), prepend!(collect(4:4:GRIDSIZE),1) )
    end

    if m == 1
        yticks( prepend!(collect(3:4:GRIDSIZE-1),0), prepend!(collect(4:4:GRIDSIZE),1) )
    end
end

# savefig("plots/" *  jobID * "-ArrowPlot.pdf", bbox_inches="tight")
# +
Base.show(io::IO, f::Float64) = @printf(io, "%.6f", f)

IGmatCtrl = Matrix{Float64}(undef, 3, 3)
UImatCtrl = Matrix{Float64}(undef, 3, 3)
FGmatCtrl = Matrix{Float64}(undef, 3, 3)
# EPmat = Matrix{Float64}(undef, 3, 3)

i = 1
λM = aλM[i]

k = 2
β = aSelStrength[k]

for (m,μ) in enumerate(aMutFactor), (j,W) in enumerate(aW)

    P = aTrajCtrl[i,j,k,m].jointProb / aTrajCtrl[i,j,k,m].NgenSample

    IGmatCtrl[m,j] = infoGain(P)
    UImatCtrl[m,j] = uselessInfo(P)
    FGmatCtrl[m,j] = fitnessGain(P,β)
#     EPmat[i,j] = EP(P)
end
# -

# j = 3 [:,j]
display( IGmatCtrl )
display( UImatCtrl )
display( FGmatCtrl )
# EPmat
