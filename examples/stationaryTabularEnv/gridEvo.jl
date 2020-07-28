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
folderName = "/home/rjku/projects/fisheria/examples/stationaryTabularEnv/"

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
    ety = EvolutionaryDynamics.TabularEvotype(
        EvolutionaryDynamics.NeutralReplication(REPCOEF),
        EvolutionaryDynamics.StandardMutation(),
        EvolutionaryDynamics.FitnessSelection(),
        grid
    )
    
#     traj = @distributed (+) for i in 1:NSAMPLES
    traj = EvolutionaryDynamics.generateTabularSystemsTrajectories(
            ety, fTbl, β, NPOP, NGENRELAX, NGENSAMPLE
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

for (i,β) in enumerate(aSelStrength), (j,μ) in enumerate(aMutFactor)

    grid = mGraphs.EdgeWeightedSquareLattice(GRIDSIZE, [ μ for i in 1:2GRIDSIZE*(GRIDSIZE-1) ]);
    Pgty = [ sum( aTraj[i,j].jointProb[g,:] ) for g in 1:size(aTraj[i,j].jointProb)[1] ] /
        aTraj[i,j].NgenSample

    p = mGraphs.VertexWeightedSquareLattice( GRIDSIZE, Pgty )
    pMat = mGraphs.matrixForm(p);

    ax = figNum.add_subplot(3,3,j + 3*(3-i))

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
