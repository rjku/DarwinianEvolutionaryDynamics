
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
Base.run(`cp varGridEvo-parameters.jl $parametersFileName`)

# print code version for the records
Base.run(`git --git-dir=/home/riccardorao/projects/fisheria/.git log --format=format:%h -1`)

#############################################################################

include("varGridEvo-parameters.jl")

aTraj = Array{EvolutionaryDynamics.VaryingEnvironmentTrajectoryData}(undef,
    length(aSelStrength), length(aMutFactor), length(aEnvTransRate))

@time for (i,β) in enumerate(aSelStrength), (j,μ) in enumerate(aMutFactor), (k,w) in enumerate(aEnvTransRate)

    grid = mGraphs.EdgeWeightedSquareLattice(GRIDSIZE, [ μ for i in 1:2GRIDSIZE*(GRIDSIZE-1) ]);
    ety = EvolutionaryDynamics.TabularEvotype(
        # EvolutionaryDynamics.NeutralReplication(REPCOEF),
        EvolutionaryDynamics.WithoutReplication(),
        EvolutionaryDynamics.StandardMutation(),
        EvolutionaryDynamics.FitnessSelection(),
        grid
    )
    transMtx = [ [ 1.0 - w, w ] [ w, 1.0 - w ] ]

	traj = @distributed (+) for s in 1:NSAMPLES
    	EvolutionaryDynamics.generateTabularSystemsTrajectories(
            ety=ety, aFitnessTbl=aFitnessTbl, selCoef=β, transMtx=transMtx, Npop=NPOP,
            # NgenRelax=convert(Int32,RELNGENRELAX/w), NgenSample=convert(Int64,RELNGENSAMPLE/w)
            NgenRelax=convert(Int32,RELNGENRELAX/w), NgenSample=NGENSAMPLE
        )
	end
    aTraj[i,j,k] = traj
end

jldopen(folderName * jobTag * "_aTraj.jld", "w") do file
	addrequire(file, EvolutionaryDynamics)
	write(file, "aTraj", aTraj)
end

#############################################################################

# possibly copy slurmLog.out somewhere else
# slurmLogFileName = jobTag * "-slurmLog.out"
# Base.run(`cp $slurmLogFileName $folderName`)
