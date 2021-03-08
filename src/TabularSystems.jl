
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

export TabularEvotype
# -------

struct TabularEnvironment{T<:Vector{<:Real}} <: AbstractTabularEnvironment
	fitnessTbl::T
	selCoef::Float64
end

IsVaryingEnvironment(::Type{<:TabularEnvironment}) = StationaryEnvironment()

export TabularEnvironment
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

"""
	PopulationTrajectoryData <: AbstractEvolutionData

	type encoding population (rather than genotype) trajectory data

	# arguments
	- nGenRelax::Int32						number of generations for the relaxation period
	- nSamples::Int32						number of samples
	- nGenSamples::Int64					number of generations between one sampling and the next

	- avePerformance::Vector{Float64}		average fitness along the trajectory
	- growthFactor::Vector{Float64}			growth coefficients along the trajectory
	- mutationFactor::Vector{Float64}		mutation coefficients along the trajectory

	- aPopCmp::Vector{Vector{Int32}}		sampled population compositions
"""
struct PopulationTrajectoryData <: AbstractEvolutionData
	nGenRelax::Int32
	nSamples::Int32
	nGenSamples::Int64

	avePerformance::Vector{Float64}
	growthFactor::Vector{Float64}
	mutationFactor::Vector{Float64}

	aPopCmp::Vector{Vector{Int32}}
end
	"""
		PopulationTrajectoryData(nGenRelax::Integer, nSamples::Integer)

		initialization of a type PopulationTrajectoryData

		# arguments
		- nGenRelax::Int32						number of generations for the relaxation period
		- nSamples::Int32						number of samples

		# notes
		- nGenSamples::Int64					set to nGenRelax * nSamples ~ 1 sample after any relaxation period
	"""
	function PopulationTrajectoryData(nGenRelax::Integer, nSamples::Integer)
		nGenSamples = Int64(nGenRelax * nSamples)

		return PopulationTrajectoryData(
			convert(Int32, nGenRelax), convert(Int32, nSamples), nGenSamples,
			Vector{Float64}(undef, nGenSamples), Vector{Float64}(undef, nGenSamples), Vector{Float64}(undef, nGenSamples),
			Vector{Vector{Int32}}(undef, nSamples)
		)
	end

import Base: +

+(trj1::TrajectoryData, trj2::TrajectoryData) = TrajectoryData(
	trj1.nGenRelax, trj1.NgenSample + trj2.NgenSample,
	trj1.avePerformance, trj1.growthFactor, trj1.mutationFactor,
	trj1.jointProb .+ trj2.jointProb
)

"""
	+(trj1::PopulationTrajectoryData, trj2::PopulationTrajectoryData)

	merges two PopulationTrajectoryData types together into trj

	- trj.nGenRelax 		= trj1.nGenRelax
	- trj.nSamples 			= trj1.nSamples + trj2.nSamples
	- trj.nGenSamples 		= trj1.nGenSamples + trj2.nGenSamples,
	- trj.avePerformance	= trj1.avePerformance
	- trj.growthFactor		= trj1.growthFactor
	- trj.mutationFactor 	= trj1.mutationFactor
	- trj.aPopCmp			= vcat( trj1.aPopCmp, trj2.aPopCmp )
"""
+(trj1::PopulationTrajectoryData, trj2::PopulationTrajectoryData) = PopulationTrajectoryData(
	trj1.nGenRelax, trj1.nSamples + trj2.nSamples, trj1.nGenSamples + trj2.nGenSamples,
	trj1.avePerformance, trj1.growthFactor, trj1.mutationFactor,
	vcat( trj1.aPopCmp, trj2.aPopCmp )
)

export TrajectoryData, PopulationTrajectoryData
# =========


function init_Population(N::Int32, ety::Tety, env::Tenv, aGty::Vector{TG}) where {
		Tety<:AbstractEvotype, Tenv<:AbstractTabularEnvironment, TG<:AbstractGenotype
	}
	N <= length(aGty) || throw(DimensionMismatch("Inconsistent Dimensions: N must be greater than or equal to length(aGty)"))

	@threads for i in 1:N fitness!(aGty[i], env) end
	return Population{Tety,Tenv,Vector{TG}}(Int32[N,N], ety, env, aGty)
end

function _fitness!(::StationaryEnvironment, gty::AbstractGenotype, env::AbstractTabularEnvironment)
	gty.aFitness .= [ exp(env.selCoef * env.fitnessTbl[gty.genome[1]]), env.fitnessTbl[gty.genome[1]] ]
end

export init_Population
#---------
# MUTATION

function _mutation!(::StandardMutation, gty::AbstractGenotype, ety::TabularEvotype)
	aAdjV, aMutProb = neighbors(ety.graph, gty.genome[1])

	# i = draw(aMutProb, rand(THREADRNG[threadid()]))
	i = draw(aMutProb, rand())

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

	# i = draw(aMutProb, rand(THREADRNG[threadid()]))
	i = draw(aMutProb, rand())
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

	# i = draw(aMutProb, rand(THREADRNG[threadid()]))
	i = draw(aMutProb, rand())
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

		# traj.aPopCmp[iSample + 1] = deepcopy(pop.aGty[1:pop.pN[2]])
		traj.aPopCmp[iSample + 1] = [ pop.aGty[i].genome[1] for i in 1:pop.pN[2] ]
	end
end

function generateTabularSystemsTrajectories(
		ety::AbstractTabularEvotype, fitnessTbl::AbstractArray{<:Real}, selCoef::Real, Npop::Integer, nGenRelax::Integer, NgenSample::Integer
	)

	env = TabularEnvironment(fitnessTbl, selCoef)
	# aGty = [ Genotype([convert(Int32, rand(THREADRNG[threadid()], 1:ety.graph.Nv))]) for i in 1:Npop ]
	aGty = [ Genotype([convert(Int32, rand(1:ety.graph.Nv))]) for i in 1:Npop ]
	pop = init_Population(convert(Int32, Npop), ety, env, aGty)

	trajData = TrajectoryData(nGenRelax, NgenSample, ety.graph.Nv)
	evolve!(pop, trajData)

	return trajData
end

"""
	generateTabularSystemsPopulationTrajectories(; <keywords arguments>)

	generates a PopulationTrajectoryData type

	# arguments
	- ety::AbstractTabularEvotype				evotype
	- fitnessTbl::AbstractVector{<:Real}		fitness value table for each genotype
	- selCoef::Real								selection strength coefficient
	- Npop::Integer								population size
	- nGenRelax::Integer						number of generations for the relaxation period
	- nSamples::Integer							number of samples
"""
function generateTabularSystemsPopulationTrajectories(;
		ety::AbstractTabularEvotype, fitnessTbl::AbstractVector{<:Real}, selCoef::Real, Npop::Integer, nGenRelax::Integer, nSamples::Integer
	)

	env = TabularEnvironment(fitnessTbl, selCoef)

	# initial genotypes are sampled according to fitness
	aSelectionProb = fitnessTbl ./ sum(fitnessTbl)
	# aGty = [ Genotype( [ convert( Int32, rand(THREADRNG[threadid()], Categorical(aSelectionProb)) ) ] ) for i in 1:Npop ]
	aGty = [ Genotype( [ convert( Int32, rand(Categorical(aSelectionProb)) ) ] ) for i in 1:Npop ]

	pop = init_Population(convert(Int32, Npop), ety, env, aGty)

	trajData = PopulationTrajectoryData(nGenRelax, nSamples)
	evolve!(pop, trajData)

	return trajData
end

export generateTabularSystemsTrajectories, generateTabularSystemsPopulationTrajectories

#-----

end # module TabularSystems
