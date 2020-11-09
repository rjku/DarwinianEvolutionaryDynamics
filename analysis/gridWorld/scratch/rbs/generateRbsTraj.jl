# -*- coding: utf-8 -*-

using ClusterManagers, Distributed
using BenchmarkTools, JLD, HDF5, DelimitedFiles

jobID = ARGS[3]
jobTag = ARGS[4] * "-" * jobID
println("\njobTag: " * jobTag )

folderName = ARGS[2]
println("folder: " * folderName )

np = parse(Int64,ARGS[1])
addprocs_slurm(np, t="48:00:00")
println("\n")

# spawn fisheria project path push at all workers
push!(LOAD_PATH, "/home/riccardorao/projects/julietta")
push!(LOAD_PATH, "/home/riccardorao/projects/fisheria/src")

for i in workers()
	@spawnat i ( push!(LOAD_PATH, "/home/riccardorao/projects/julietta") )
	@spawnat i ( push!(LOAD_PATH, "/home/riccardorao/projects/fisheria/src") )
end

@everywhere using mGraphs, mUtils
@everywhere import EvolutionaryDynamics

# copy parameter files for the records
parametersFileName = folderName * jobTag * "-parameters.jl"
Base.run(`cp rbsGridEvo-parameters.jl $parametersFileName`)

# print code version for the records
Base.run(`git --git-dir=/home/riccardorao/projects/fisheria/.git log --format=format:%h -1`)

# ############################################################################

include("rbsGridEvo-parameters.jl")

aTraj = Array{EvolutionaryDynamics.PopulationTrajectoryData}(undef, length(aSelStrength), length(aMutFactor))

@time for (i,β) in enumerate(aSelStrength), (j,μ) in enumerate(aMutFactor)
    grid = mGraphs.EdgeWeightedSquareLattice(GRIDSIZE, [ μ for i in 1:2GRIDSIZE*(GRIDSIZE-1) ]);
    ety = EvolutionaryDynamics.TabularEvotype(
        EvolutionaryDynamics.WithoutReplication(),
        EvolutionaryDynamics.StandardMutation(),
        EvolutionaryDynamics.FitnessSelection(),
        grid
    )

	traj = @distributed (+) for i in 1:NSAMPLES
    	EvolutionaryDynamics.generateTabularSystemsPopulationTrajectories(
			ety=ety, fitnessTbl=fTbl, selCoef=β, Npop=NPOP, nGenRelax=trunc(Int32,DIMGSPACE/μ), nSamples=NSAMPLESPERTRAJ
        )
	end
    aTraj[i,j] = traj
end

jldopen(folderName * jobTag * "_aTraj.jld", "w") do file
	addrequire(file, EvolutionaryDynamics)
	write(file, "aTraj", aTraj)
end

# ############################################################################

# possibly copy slurmLog.out somewhere else
# slurmLogFileName = jobTag * "-slurmLog.out"
# Base.run(`cp $slurmLogFileName $folderName`)
