
module TabularSystems

# using ..EvolutionaryDynamics
using  ..Types
import ..Types: IsVaryingEnvironment

using  ..Methods: ancestry, replication!, selection!, fitness!, evolveEnvironment!, draw
using  ..Methods: THREADRNG
import ..Methods: mutation!, _mutation!, _fitness!

using mUtils, mGraphs

import Statistics: mean
import Distributions: Categorical

using Base.Threads, Random
import Future

# threads random generators initialization
# const THREADRNG = let m = MersenneTwister(1)
#             [m; accumulate(Future.randjump, fill(big(10)^20, nthreads()-1), init=m)]
#         end;

abstract type AbstractTabularEvotype <: AbstractEvotype end
abstract type AbstractTabularEnvironment <: AbstractEnvironment end
# -------

struct TabularEvotype{Trep<:ReplicationType,Tmut<:MutationType,Tsel<:SelectionType,Tgrh<:AbstractGraph} <: AbstractTabularEvotype
	repType::Trep
	mutType::Tmut
	selType::Tsel
	graph::Tgrh
end

struct VaryingTabularEvotypeI{T<:AbstractGraph} <: AbstractTabularEvotype
	pRepFactor::Vector{Float64}
	minRepCoef::Float64
	mutCoef::Float64
	λM::Float64						# mutliplication factor for mutation on the grid
	graph::T
end

export TabularEvotype, VaryingTabularEvotypeI
# -------

struct TabularEnvironment{T<:Vector{<:Real}} <: AbstractTabularEnvironment
	fitnessTbl::T
	selCoef::Float64
end

struct VaryingTabularEnvironment{T<:Vector{<:Vector{<:Real}}} <: AbstractTabularEnvironment
	fitnessTbl::T
	selCoef::Float64

	transMtx::Matrix{Float64}		# transition matrix
	envState::Vector{Int32}			# index current environment
end
	VaryingTabularEnvironment(fitnessTbl, selCoef, transMtx) = VaryingTabularEnvironment(fitnessTbl, selCoef, transMtx, [Int32(1)])

	IsVaryingEnvironment(::Type{<:TabularEnvironment}) = StationaryEnvironment()
	IsVaryingEnvironment(::Type{<:VaryingTabularEnvironment}) = MarkovianEnvironment()

export TabularEnvironment, VaryingTabularEnvironment
# -------

scaledLogFitness(gty::AbstractGenotype) = gty.aFitness[2]
export scaledLogFitness
# -------

struct TrajectoryData <: AbstractEvolutionData
	nGenRelax::Int32
	NgenSample::Int64

	avePerformance::Array{Float64,1}
	growthFactor::Array{Float64,1}
	mutationFactor::Array{Float64,1}

	jointProb::Array{Int64,2}
end
	TrajectoryData(nGenRelax::Integer, NgenSample::Integer, cardG::Integer) = TrajectoryData( Int32(nGenRelax), Int64(NgenSample),
		Array{Float64}(undef,NgenSample), Array{Float64}(undef,NgenSample), Array{Float64}(undef,NgenSample), zeros(Int64,cardG,cardG)
	)

struct PopulationTrajectoryData{TaG<:Vector{<:AbstractGenotype}} <: AbstractEvolutionData
	nGenRelax::Int32
	nSamples::Int32
	nGenSamples::Int64

	avePerformance::Vector{Float64}
	growthFactor::Vector{Float64}
	mutationFactor::Vector{Float64}

	aPopCmp::Vector{TaG}
end
	function PopulationTrajectoryData(nGenRelax::Integer, nSamples::Integer, gtyType::Type{<:AbstractGenotype})
		nGenSamples = Int64(nGenRelax * nSamples)

		return PopulationTrajectoryData(
			convert(Int32, nGenRelax), convert(Int32, nSamples), nGenSamples,
			Vector{Float64}(undef, nGenSamples), Vector{Float64}(undef, nGenSamples), Vector{Float64}(undef, nGenSamples),
			Vector{Vector{gtyType}}(undef, nSamples)
		)
	end

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

+(trj1::TrajectoryData, trj2::TrajectoryData) = TrajectoryData(
	trj1.nGenRelax, trj1.NgenSample + trj2.NgenSample,
	trj1.avePerformance, trj1.growthFactor, trj1.mutationFactor,
	trj1.jointProb .+ trj2.jointProb
)

+(trj1::PopulationTrajectoryData, trj2::PopulationTrajectoryData) = PopulationTrajectoryData(
	trj1.nGenRelax, trj1.nSamples + trj2.nSamples, trj1.nGenSamples + trj2.nGenSamples,
	trj1.avePerformance, trj1.growthFactor, trj1.mutationFactor,
	vcat( trj1.aPopCmp, trj2.aPopCmp )
)

+(trj1::VaryingEnvironmentTrajectoryData, trj2::VaryingEnvironmentTrajectoryData) = VaryingEnvironmentTrajectoryData(
	trj1.nGenRelax, trj1.NgenSample + trj2.NgenSample,
	trj1.avePerformance, trj1.growthFactor, trj1.mutationFactor, trj1.envState,
	trj1.jointProb .+ trj2.jointProb
)

export TrajectoryData, PopulationTrajectoryData, VaryingEnvironmentTrajectoryData
# =========


function init_Population(N::Int32, ety::Tety, env::Tenv, aGty::Vector{TG}) where {
		Tety<:AbstractEvotype, Tenv<:AbstractTabularEnvironment, TG<:AbstractGenotype
	}
	N <= length(aGty) || throw(DimensionMismatch("Inconsistent Dimensions: N must be greater than or equal to length(aGty)"))

	@threads for i in 1:N fitness!(aGty[i], env) end
	return Population{Tety,Tenv,Vector{TG}}(Int32[N,N], ety, env, aGty)
end

export init_Population

function _fitness!(::StationaryEnvironment, gty::AbstractGenotype, env::AbstractTabularEnvironment)
	gty.aFitness .= [ exp(env.selCoef * env.fitnessTbl[gty.genome[1]]), env.fitnessTbl[gty.genome[1]] ]
end

function _fitness!(::MarkovianEnvironment, gty::AbstractGenotype, env::AbstractTabularEnvironment)
	gty.aFitness .= [ exp(env.selCoef * env.fitnessTbl[envState(env)][gty.genome[1]]), env.fitnessTbl[envState(env)][gty.genome[1]] ]
end

#---------
# MUTATION

function _mutation!(::StandardMutation, gty::AbstractGenotype, ety::TabularEvotype)
	aAdjV, aMutProb = neighbors(ety.graph, gty.genome[1])

	i = draw(aMutProb, rand(THREADRNG[threadid()]))

	if i > 0
		gty.genome[1] = aAdjV[i]
		return Int32(1)
	else
		return Int32(0)
	end
end

_mutation!(::RescaledMutation, gty::AbstractGenotype, ety::TabularEvotype) = _mutation!(ety.repType, gty, ety)

function _mutation!(repType::NeutralReplication, gty::AbstractGenotype, ety::TabularEvotype)
	aAdjV, aMutProb = neighbors(ety.graph, gty.genome[1])
	aMutProb .= aMutProb / ( 1.0 + ety.repType.repCoef )

	i = draw(aMutProb, rand(THREADRNG[threadid()]))
	if i > 0
		gty.genome[1] = aAdjV[i]
		return Int32(1)
	else
		return Int32(0)
	end
end

function _mutation!(repType::FitnessReplication, gty::AbstractGenotype, ety::TabularEvotype)
	aAdjV, aMutProb = neighbors(ety.graph, gty.genome[1])
	aMutProb .= aMutProb / ( 1.0 + ety.repType.repFactor * fitness(gty) )

	i = draw(aMutProb, rand(THREADRNG[threadid()]))
	if i > 0
		gty.genome[1] = aAdjV[i]
		return Int32(1)
	else
		return Int32(0)
	end
end

#------

function evoStep!(pop::AbstractPopulation)
	gf = replication!(pop)
	nMut = mutation!(pop)
	selection!(pop)

	return gf, nMut
end

function evoStep!(pop::AbstractPopulation, ancestryVec::Vector{<:Integer})
	gf = replication!(pop, ancestryVec)
	nMut = mutation!(pop)
	iGtySelected = selection!(pop)

	return gf, nMut, ancestryVec[iGtySelected]
end

function evolve!(pop::AbstractPopulation, traj::TrajectoryData)
	for gen in 1:traj.nGenRelax
		evoStep!(pop)
	end

	ancestryVec = ancestry(pop.ety, pop.pN[2])

	for gen in 1:traj.NgenSample
		aGtyPast = deepcopy(pop.aGty)

		traj.growthFactor[gen], traj.mutationFactor[gen], iGtyPast = evoStep!(pop, ancestryVec)
		traj.avePerformance[gen] = mean([ scaledLogFitness(pop.aGty[i]) for i in 1:pop.pN[2] ])

		traj.jointProb[ pop.aGty[1].genome[1], aGtyPast[iGtyPast].genome[1] ] += 1
	end
end

function evolve!(pop::AbstractPopulation, traj::PopulationTrajectoryData)
	for iSample in 0:traj.nSamples-1
		for gen in 1:traj.nGenRelax
			traj.growthFactor[gen + iSample * traj.nGenRelax], traj.mutationFactor[gen + iSample * traj.nGenRelax] = evoStep!(pop)
			traj.avePerformance[gen + iSample * traj.nGenRelax] = mean([ scaledLogFitness(pop.aGty[i]) for i in 1:pop.pN[2] ])
		end

		traj.aPopCmp[iSample + 1] = deepcopy(pop.aGty)
	end
end

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

function generateTabularSystemsTrajectories(
		ety::AbstractTabularEvotype, fitnessTbl::AbstractArray{<:Real}, selCoef::Real, Npop::Integer, nGenRelax::Integer, NgenSample::Integer
	)

	env = TabularEnvironment(fitnessTbl, selCoef)
	aGty = [ Genotype([convert(Int32, rand(THREADRNG[threadid()], 1:ety.graph.Nv))]) for i in 1:Npop ]
	pop = init_Population(convert(Int32, Npop), ety, env, aGty)

	trajData = TrajectoryData(nGenRelax, NgenSample, ety.graph.Nv)
	evolve!(pop, trajData)

	return trajData
end

function generateTabularSystemsPopulationTrajectories(;
		ety::AbstractTabularEvotype, fitnessTbl::AbstractVector{<:Real}, selCoef::Real, Npop::Integer, nGenRelax::Integer, nSamples::Integer
	)

	env = TabularEnvironment(fitnessTbl, selCoef)

	# genotypes are generated according to fitness
	aSelectionProb = fitnessTbl ./ sum(fitnessTbl)
	aGty = [ Genotype( [ convert( Int32, rand(THREADRNG[threadid()], Categorical(aSelectionProb)) ) ] ) for i in 1:Npop ]

	pop = init_Population(convert(Int32, Npop), ety, env, aGty)

	trajData = PopulationTrajectoryData(nGenRelax, nSamples, Genotype)
	evolve!(pop, trajData)

	return trajData
end

function generateTabularSystemsTrajectories(;
		ety::AbstractTabularEvotype, aFitnessTbl::Vector{<:AbstractArray{<:Real}}, selCoef::Real, transMtx::Matrix{<:Real},
		Npop::Integer, nGenRelax::Integer, NgenSample::Integer
	)

	env = VaryingTabularEnvironment(aFitnessTbl, selCoef, transMtx)
	aGty = [ Genotype([convert(Int32, rand(THREADRNG[threadid()], 1:ety.graph.Nv ))]) for i in 1:Npop ];
	pop = init_Population(convert(Int32, Npop), ety, env, aGty)

	trajData = VaryingEnvironmentTrajectoryData(nGenRelax, NgenSample, ety.graph.Nv, length(aFitnessTbl))
	evolve!(pop, trajData)

	return trajData
end

export generateTabularSystemsTrajectories, generateTabularSystemsPopulationTrajectories

#-----

function mutation!(gty::AbstractGenotype, ety::VaryingTabularEvotypeI, env::AbstractTabularEnvironment)
	aAdjV, aMutProb = neighbors(ety.G, gty.genome[1])
	sqrtρhat = sqrt( 1.0 + ety.minRepCoef + ety.pRepFactor[1]*gty.aFitness[2] )
	aMutProb .= ( aMutProb .* [ i == gty.genome[2] ? ety.λM : 1.0 for i in eachindex(aMutProb) ] ) / sqrtρhat

	Nmut::Int32 = 0

	# site mutation
	r::Float64 = rand(THREADRNG[threadid()])
	i = draw(aMutProb,r)
	if i > 0
		gty.genome[1] = aAdjV[i]
		fitness!(gty,env)
		Nmut += 1
	end

	# mutation boost direction
	r = rand(THREADRNG[threadid()])
	if r <= ety.mutCoef / sqrtρhat
		gty.genome[2] = typeof(gty.genome[2])(rand(THREADRNG[threadid()], 0:4))		# <--- you are assuming the square grid!
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

function generateVaryingTabularSystemsITrajectories(
		repFactor, minRepCoef, mutCoef, λM, graph::AbstractGraph, aFtb, repStrength, selStrength, transitionMatrix,
		Npop::Integer, nGenRelax::Integer, NgenSample::Integer
	)

	ety = VaryingTabularEvotypeI([Float64(repFactor)], Float64(minRepCoef), Float64(mutCoef), Float64(λM), graph)
	env = VaryingTabularEnvironment(aFtb, [repStrength,selStrength], transitionMatrix)
	aGty = [ EvolutionaryDynamics.Genotype( Int32[rand(THREADRNG[threadid()], 1:graph.Nv ), 0] ) for i in 1:Npop ];
	pop = init_Population( Int32(Npop), ety, env, aGty )

	#  assuming that you deal with a square grid
	traj = VaryingEnvironmentTrajectoryData( nGenRelax, NgenSample, 5graph.Nv, length(aFtb) )
	evolveOld!(pop, traj)

	return traj
end

export generateVaryingTabularSystemsITrajectories

end # module TabularSystems
