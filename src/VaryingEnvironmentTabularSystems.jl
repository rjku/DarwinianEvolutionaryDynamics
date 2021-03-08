
module VaryingEnvironmentTabularSystems

# broken: dependencies must be fixed

#-------

struct VaryingTabularEvotypeI{T<:AbstractGraph} <: AbstractTabularEvotype
	pRepFactor::Vector{Float64}
	minRepCoef::Float64
	mutCoef::Float64
	λM::Float64						# mutliplication factor for mutation on the grid
	graph::T
end

export VaryingTabularEvotypeI
#-------

struct VaryingTabularEnvironment{T<:Vector{<:Vector{<:Real}}} <: AbstractTabularEnvironment
	fitnessTbl::T
	selCoef::Float64

	transMtx::Matrix{Float64}		# transition matrix
	envState::Vector{Int32}			# index current environment
end
	VaryingTabularEnvironment(fitnessTbl, selCoef, transMtx) = VaryingTabularEnvironment(fitnessTbl, selCoef, transMtx, [Int32(1)])

IsVaryingEnvironment(::Type{<:VaryingTabularEnvironment}) = MarkovianEnvironment()

export VaryingTabularEnvironment
#-------

struct VaryingEnvironmentTrajectoryData <: AbstractEvolutionData
	nGenRelax::Int32
	NgenSample::Int64

	avePerformance::Array{Float64,1}
	growthFactor::Array{Float64,1}
	mutationFactor::Array{Float64,1}
	envState::Array{Int32,1}

	jointProb::Array{Int64,4}
end
	VaryingEnvironmentTrajectoryData(nGenRelax::Integer,NgenSample::Integer,cardG::Integer,cardE::Integer) = VaryingEnvironmentTrajectoryData(
		Int32(nGenRelax), Int64(NgenSample),
		Array{Float64}(undef, NgenSample), Array{Float64}(undef, NgenSample), Array{Float64}(undef, NgenSample),
		Array{Int32}(undef, NgenSample), zeros(Int64, cardG, cardE, cardG, cardE)
	)

import Base: +

+(trj1::VaryingEnvironmentTrajectoryData, trj2::VaryingEnvironmentTrajectoryData) = VaryingEnvironmentTrajectoryData(
	trj1.nGenRelax, trj1.NgenSample + trj2.NgenSample,
	trj1.avePerformance, trj1.growthFactor, trj1.mutationFactor, trj1.envState,
	trj1.jointProb .+ trj2.jointProb
)

export VaryingEnvironmentTrajectoryData
#-------


function _fitness!(::MarkovianEnvironment, gty::AbstractGenotype, env::AbstractTabularEnvironment)
	gty.aFitness .= [ exp(env.selCoef * env.fitnessTbl[envState(env)][gty.genome[1]]), env.fitnessTbl[envState(env)][gty.genome[1]] ]
end

#-------

function evolve!(pop::AbstractPopulation, traj::VaryingEnvironmentTrajectoryData)
	for gen in 1:traj.nGenRelax
		if evolveEnvironment!(pop.env)
			fitness!(pop)
		end
		evoStep!(pop)
	end

	ancestryVec = ancestry(pop.ety, pop.pN[2])

	for gen in 1:traj.NgenSample
		aGtyPast = deepcopy(pop.aGty)
		envStatePast = envState(pop.env)

		if evolveEnvironment!(pop.env)
			fitness!(pop)
		end

		traj.growthFactor[gen], traj.mutationFactor[gen], iGtyPast = evoStep!(pop, ancestryVec)
		traj.avePerformance[gen] = mean([ scaledLogFitness(pop.aGty[i]) for i in 1:pop.pN[2] ])
		traj.envState[gen] = envState(pop.env)

		traj.jointProb[ pop.aGty[1].genome[1], traj.envState[gen], aGtyPast[iGtyPast].genome[1], envStatePast ] += 1
	end
end

function generateTabularSystemsTrajectories(;
		ety::AbstractTabularEvotype, aFitnessTbl::Vector{<:AbstractArray{<:Real}}, selCoef::Real, transMtx::Matrix{<:Real},
		Npop::Integer, nGenRelax::Integer, NgenSample::Integer
	)

	env = VaryingTabularEnvironment(aFitnessTbl, selCoef, transMtx)
	# aGty = [ Genotype([convert(Int32, rand(THREADRNG[threadid()], 1:ety.graph.Nv ))]) for i in 1:Npop ];
	aGty = [ Genotype([convert(Int32, rand(1:ety.graph.Nv ))]) for i in 1:Npop ];
	pop = init_Population(convert(Int32, Npop), ety, env, aGty)

	trajData = VaryingEnvironmentTrajectoryData(nGenRelax, NgenSample, ety.graph.Nv, length(aFitnessTbl))
	evolve!(pop, trajData)

	return trajData
end

function evolveOld!(pop::AbstractPopulation, traj::VaryingEnvironmentTrajectoryData)
	for gen in 1:traj.nGenRelax
		gmsEvoStep!(pop)
	end

	ancestry = [ i for i in 1:length(pop.aGty) ]

	for gen in 1:traj.NgenSample
		aGtyPast = deepcopy(pop.aGty)
		envStatePast = pop.env.envState[1]

		traj.growthFactor[gen], traj.mutationFactor[gen], iGtyPast = gmsEvoStep!(pop, ancestry)
		traj.envState[gen] = pop.env.envState[1]
		traj.avePerformance[gen] = mean( [pop.aGty[i].aF[1] for i in 1:pop.pN[2]] )

		traj.jointProb[ coord(pop.aGty[1], pop.ety.graph.Nv), traj.envState[gen], coord(aGtyPast[iGtyPast], pop.ety.graph.Nv), envStatePast ] += 1
	end
end

function mutation!(gty::AbstractGenotype, ety::VaryingTabularEvotypeI, env::AbstractTabularEnvironment)
	aAdjV, aMutProb = neighbors(ety.G, gty.genome[1])
	sqrtρhat = sqrt( 1.0 + ety.minRepCoef + ety.pRepFactor[1]*gty.aFitness[2] )
	aMutProb .= ( aMutProb .* [ i == gty.genome[2] ? ety.λM : 1.0 for i in eachindex(aMutProb) ] ) / sqrtρhat

	Nmut::Int32 = 0

	# site mutation
	# i = draw(aMutProb, rand(THREADRNG[threadid()]))
	i = draw(aMutProb, rand())
	if i > 0
		gty.genome[1] = aAdjV[i]
		fitness!(gty,env)
		Nmut += 1
	end

	# mutation boost direction
	# if rand(THREADRNG[threadid()]) <= ety.mutCoef / sqrtρhat
	if rand() <= ety.mutCoef / sqrtρhat
		# gty.genome[2] = typeof(gty.genome[2])(rand(THREADRNG[threadid()],0:4))		# <--- you are assuming the square grid!
		gty.genome[2] = typeof(gty.genome[2])(rand(0:4))								# <--- you are assuming the square grid!
		Nmut += 1
	end

	return Nmut
end

coord(gty::AbstractGenotype, cardGty1::Integer) = gty.genome[1] + cardGty1 * gty.genome[2]

function gmsEvoStep!(pop::Population{<:AbstractEvotype,<:VaryingTabularEnvironment,<:Vector{<:AbstractGenotype}}, ancestry)
	hasChanged = evolveEnvironment!(pop.env)
	if hasChanged
		fitness!(pop)
	end
	gf = replication!(pop, ancestry)
	Nmut = mutation!(pop)
	iGtySelected = _selection!(FitnessSelection(), pop)

	return gf, Nmut, ancestry[iGtySelected]
end

function generateVaryingTabularSystemsITrajectories(
		repFactor, minRepCoef, mutCoef, λM, graph::AbstractGraph, aFtb, repStrength, selStrength, transitionMatrix,
		Npop::Integer, nGenRelax::Integer, NgenSample::Integer
	)

	ety = VaryingTabularEvotypeI([Float64(repFactor)], Float64(minRepCoef), Float64(mutCoef), Float64(λM), graph)
	env = VaryingTabularEnvironment(aFtb, [repStrength,selStrength], transitionMatrix)
	# aGty = [ EvolutionaryDynamics.Genotype( Int32[rand(THREADRNG[threadid()], 1:graph.Nv ), 0] ) for i in 1:Npop ];
	aGty = [ EvolutionaryDynamics.Genotype( Int32[rand(1:graph.Nv ), 0] ) for i in 1:Npop ];
	pop = init_Population( Int32(Npop), ety, env, aGty )

	#  assuming that you deal with a square grid
	traj = VaryingEnvironmentTrajectoryData( nGenRelax, NgenSample, 5graph.Nv, length(aFtb) )
	evolveOld!(pop, traj)

	return traj
end

export generateVaryingTabularSystemsITrajectories

end # module VaryingEnvironmentTabularSystems
