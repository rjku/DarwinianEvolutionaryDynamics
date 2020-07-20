
module TabularSystems

# using ..EvolutionaryDynamics
using ..Types
import ..Methods: replication!, mutation!, _selection!

using mUtils, mGraphs

import Statistics: mean

using Base.Threads, Random
import Future

# threads random generators initialization
const THREADRNG = let m = MersenneTwister(1)
            [m; accumulate(Future.randjump, fill(big(10)^20, nthreads()-1), init=m)]
        end;

abstract type AbstractTabularEvotype <: AbstractEvotype end
abstract type AbstractTabularEnvironment <: AbstractEnvironment end
# -------

struct TabularEvotype{T<:AbstractGraph} <: AbstractTabularEvotype
	pRepFactor::Vector{Float64}
	minRepCoef::Float64
	graph::T
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
	aSelCoef::Vector{Float64}
end

struct VaryingTabularEnvironment{T<:Vector{<:Vector{<:Real}}} <: AbstractTabularEnvironment
	fitnessTbl::T
	aSelCoef::Vector{Float64}

	transMtx::Matrix{Float64}		# transition matrix
	envState::Vector{Int32}			# index current environment
end
	VaryingTabularEnvironment(fitnessTbl, aSelCoef, transMtx) = VaryingTabularEnvironment(fitnessTbl, aSelCoef, transMtx, [Int32(1)])

	IsVaryingEnvironment(::Type{<:TabularEnvironment}) = StationaryEnvironment()
	IsVaryingEnvironment(::Type{<:VaryingTabularEnvironment}) = MarkovianEnvironment()

export TabularEnvironment, VaryingTabularEnvironment
# -------

scaledLogFitness(gty::AbstractGenotype) = gty.aFitness[2]
export scaledLogFitness
# -------

struct TrajectoryData <: AbstractEvolutionData
	NgenRelax::Int32
	NgenSample::Int64

	avePerformance::Array{Float64,1}
	growthFactor::Array{Float64,1}
	mutationFactor::Array{Float64,1}

	jointProb::Array{Int64,2}
end
	TrajectoryData(NgenRelax::Integer, NgenSample::Integer, cardG::Integer) = TrajectoryData( Int32(NgenRelax), Int64(NgenSample),
		Array{Float64}(undef,NgenSample), Array{Float64}(undef,NgenSample), Array{Float64}(undef,NgenSample), zeros(Int64,cardG,cardG)
	)

struct VaryingEnvironmentTrajectoryData <: AbstractEvolutionData
	NgenRelax::Int32
	NgenSample::Int64

	avePerformance::Array{Float64,1}
	growthFactor::Array{Float64,1}
	mutationFactor::Array{Float64,1}
	envState::Array{Int32,1}

	jointProb::Array{Int64,4}
end
	VaryingEnvironmentTrajectoryData(NgenRelax::Integer,NgenSample::Integer,cardG::Integer,cardE::Integer) = VaryingEnvironmentTrajectoryData( Int32(NgenRelax), Int64(NgenSample),
		Array{Float64}(undef,NgenSample), Array{Float64}(undef,NgenSample), Array{Float64}(undef,NgenSample), Array{Int32}(undef,NgenSample), zeros(Int64,cardG,cardE,cardG,cardE)
	)

import Base: +

+(trj1::TrajectoryData,trj2::TrajectoryData) = TrajectoryData(
	trj1.NgenRelax, trj1.NgenSample + trj2.NgenSample, trj1.avePerformance, trj1.growthFactor, trj1.mutationFactor, trj1.jointProb .+ trj2.jointProb
)

+(trj1::VaryingEnvironmentTrajectoryData,trj2::VaryingEnvironmentTrajectoryData) = VaryingEnvironmentTrajectoryData(
	trj1.NgenRelax, trj1.NgenSample + trj2.NgenSample, trj1.avePerformance, trj1.growthFactor, trj1.mutationFactor, trj1.envState, trj1.jointProb .+ trj2.jointProb
)

export TrajectoryData, VaryingEnvironmentTrajectoryData
# =========


function init_Population(N::Int32, ety::Tety, env::Tenv, aGty::Vector{TG}) where { Tety<:AbstractEvotype, Tenv<:AbstractTabularEnvironment, TG<:AbstractGenotype }
	N <= length(aGty) || throw(DimensionMismatch("Inconsistent Dimensions: N must be greater than or equal to length(aGty)"))
	# @threads
	for i in 1:N fitness!(aGty[i], env) end
	return Population{Tety,Tenv,Vector{TG}}(Int32[N,N], ety, env, aGty)
end

export init_Population

fitness!(gty::AbstractGenotype, env::Tenv) where Tenv <: AbstractTabularEnvironment = _fitness!(IsVaryingEnvironment(Tenv), gty, env)

function _fitness!(::StationaryEnvironment, gty, env)
	gty.aFitness .= [ exp(env.fitnessTbl[gty.genome[1]] * env.aSelCoef[1]), env.fitnessTbl[gty.genome[1]] ]
end

function _fitness!(::MarkovianEnvironment, gty, env)
	gty.aFitness .= [ exp(env.fitnessTbl[envState(env)][gty.genome[1]]*env.aSelCoef[1]), env.fitnessTbl[envState(env)][gty.genome[1]] ]
end

# -----------
#  MUTATION

function mutation!(gty::AbstractGenotype, ety::TabularEvotype, env::AbstractTabularEnvironment)
	aAdjV, aMutProb = neighbors(ety.graph, gty.genome[1])
	aMutProb .= aMutProb / ( 1.0 + ety.minRepCoef + ety.pRepFactor[1]*gty.aFitness[2] )

	r::Float64 = rand(THREADRNG[threadid()])
	i = draw(aMutProb,r)

	if i > 0
		gty.genome[1] = aAdjV[i]
		fitness!(gty,env)
		return Int32(1)
	else
		return Int32(0)
	end
end

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
		gty.genome[2] = typeof(gty.genome[2])(rand(0:4))		# <--- you are assuming the square grid!
		Nmut += 1
	end

	return Nmut
end

function updateFitness!(pop::AbstractPopulation)
	for i in 1:pop.pN[2]
		fitness!(pop.aGty[i], pop.env)
	end
end

coord(gty::AbstractGenotype, cardGty1::Integer) = gty.genome[1] + cardGty1 * gty.genome[2]

function gmsEvoStep!(pop::Population{<:AbstractEvotype,<:TabularEnvironment,<:Vector{<:AbstractGenotype}})
	replication!(pop)
	mutation!(pop)
	_selection!(FitnessSelection(), pop)
end

function gmsEvoStep!(pop::Population{<:AbstractEvotype,<:TabularEnvironment,<:Vector{<:AbstractGenotype}}, ancestry)
	gf = replication!(pop)
	Nmut = mutation!(pop)
	iGtySelected = _selection!(FitnessSelection(), pop)

	return gf, Nmut, ancestry[iGtySelected]
end

function evolve!(pop::AbstractPopulation, traj::TrajectoryData)
	for gen in 1:traj.NgenRelax
		gmsEvoStep!(pop)
	end

	ancestry = [ i for i in 1:length(pop.aGty) ]

	for gen in 1:traj.NgenSample
		aGtyPast = deepcopy(pop.aGty)

		traj.growthFactor[gen], traj.mutationFactor[gen], iGtyPast = gmsEvoStep!(pop, ancestry)
		traj.avePerformance[gen] = mean( [ scaledLogFitness(pop.aGty[i]) for i in 1:pop.pN[2]] )

		traj.jointProb[ pop.aGty[1].genome[1], aGtyPast[iGtyPast].genome[1] ] += 1
	end
end

function generateTabularSystemsTrajectories(
		repFactor, minRepCoef, mutCoef, graph::AbstractGraph, fitnessTbl, repStrength, selStrength, Npop::Integer, NgenRelax::Integer, NgenSample::Integer
	)

	ety = TabularEvotype([Float64(repFactor)], Float64(minRepCoef), graph)
	env = TabularEnvironment(fitnessTbl, [repStrength, selStrength])
	aGty = [ Genotype( Int32[rand( 1:graph.Nv )] ) for i in 1:Npop ];
	pop = init_Population( Int32(Npop), ety, env, aGty )

	#  assuming that you deal with a square grid
	traj = TrajectoryData( NgenRelax, NgenSample, graph.Nv )
	evolve!(pop, traj)

	return traj
end

function gmsEvoStep!(pop::Population{<:AbstractEvotype,<:VaryingTabularEnvironment,<:Vector{<:AbstractGenotype}}, ancestry)
	hasChanged = evolveEnvironment!(pop.env)
	if hasChanged
		updateFitness(pop)
	end
	gf = replication!(pop, ancestry)
	Nmut = mutation!(pop)
	iGtySelected = _selection!(FitnessSelection(), pop)

	return gf, Nmut, ancestry[iGtySelected]
end

function evolve!(pop::AbstractPopulation, traj::VaryingEnvironmentTrajectoryData)
	for gen in 1:traj.NgenRelax
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
		repFactor, minRepCoef, mutCoef, λM, graph::AbstractGraph, aFtb, repStrength, selStrength, transitionMatrix, Npop::Integer, NgenRelax::Integer, NgenSample::Integer
	)

	ety = VaryingTabularEvotypeI([Float64(repFactor)], Float64(minRepCoef), Float64(mutCoef), Float64(λM), graph)
	env = VaryingTabularEnvironment(aFtb, [repStrength,selStrength], transitionMatrix)
	aGty = [ EvolutionaryDynamics.Genotype( Int32[rand( 1:graph.Nv ), 0] ) for i in 1:Npop ];
	pop = init_Population( Int32(Npop), ety, env, aGty )

	#  assuming that you deal with a square grid
	traj = VaryingEnvironmentTrajectoryData( NgenRelax, NgenSample, 5graph.Nv, length(aFtb) )
	evolve!(pop, traj)

	return traj
end

export generateTabularSystemsTrajectories, generateVaryingTabularSystemsITrajectories

end # module TabularSystems
