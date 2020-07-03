
# =======
#  TYPES

abstract type AbstractTabularEnvironment <: AbstractEnvironment end

struct TabularEnvironment{T<:Vector{<:Vector{<:Real}}} <: AbstractTabularEnvironment
	fTb::T
	aSelCoef::Vector{Float64}
end

struct VaryingTabularEnvironment{T<:Vector{<:Vector{<:Real}}} <: AbstractTabularEnvironment
	fTb::T
	aSelCoef::Vector{Float64}

	W::Matrix{Float64}			# transition matrix
	e::Vector{Int32}			# index current environment
end
	VaryingTabularEnvironment(fTb,aSelCoef,W) = VaryingTabularEnvironment(fTb,aSelCoef,W,[Int32(1)])

export VaryingTabularEnvironment, TabularEnvironment
# -------

using mGraphs

struct TabularEvotype{T<:AbstractGraph} <: atGtyMutEty
	pRepFactor::Vector{Float64}
	minRepCoef::Float64
	G::T
end

struct VaryingTabularEvotypeI{T<:AbstractGraph} <: atGtyMutEty
	pRepFactor::Vector{Float64}
	minRepCoef::Float64
	mutCoef::Float64
	λM::Float64						# mutliplication factor for mutation on the grid
	G::T
end

export TabularEvotype, VaryingTabularEvotypeI
# -------

struct TabularGenotype{TG<:AbstractArray} <: AbstractGenotype
	G::TG
	aF::Vector{Float64}
end
	TabularGenotype(G) = TabularGenotype(G,[0.0,0.0,0.0])

Base.copy(gty::TabularGenotype) = TabularGenotype( deepcopy(gty.G), deepcopy(gty.aF) )

struct tVarGty{TG<:AbstractArray} <: AbstractGenotype
	G::TG
	aF::Vector{Float64}
end
	tVarGty(G) = tVarGty(G,[0.0,0.0,0.0])

Base.copy(gty::tVarGty) = tVarGty( deepcopy(gty.G), deepcopy(gty.aF) )

export TabularGenotype
# -------

struct EvoData <: AbstractEvolutionData
	Ngen::Int32
	avePerformance::Array{Float64,1}
	growthFactor::Array{Float64,1}
	mutationFactor::Array{Float64,1}
end
	EvoData(Ngen::Integer) = EvoData( Int32(Ngen), Array{Float64}(undef,Ngen), Array{Float64}(undef,Ngen), Array{Float64}(undef,Ngen) )

struct TrajectoryData <: AbstractEvolutionData
	NgenRelax::Int32
	NgenSample::Int64

	avePerformance::Array{Float64,1}
	growthFactor::Array{Float64,1}
	mutationFactor::Array{Float64,1}
	envState::Array{Int32,1}

	jointProb::Array{Int64,4}
end
	TrajectoryData(NgenRelax::Integer,NgenSample::Integer,cardG::Integer,cardE::Integer) = TrajectoryData( Int32(NgenRelax), Int64(NgenSample),
		Array{Float64}(undef,NgenSample), Array{Float64}(undef,NgenSample), Array{Float64}(undef,NgenSample), Array{Int32}(undef,NgenSample),
		zeros(Int64,cardG,cardE,cardG,cardE)
	)

import Base: +

+(trj1::TrajectoryData,trj2::TrajectoryData) = TrajectoryData(
	trj1.NgenRelax, trj1.NgenSample + trj2.NgenSample, trj1.avePerformance, trj1.growthFactor, trj1.mutationFactor, trj1.envState, trj1.jointProb .+ trj2.jointProb
)

export EvoData, TrajectoryData


# =========
# METHODS

using Base.Threads, Random
import Future

# threads random generators initialization
const THREADRNG = let m = MersenneTwister(1)
            [m; accumulate(Future.randjump, fill(big(10)^20, nthreads()-1), init=m)]
        end;

function init_Population( N::Int32,ety::Tety,env::Tenv,aGty::Vector{TG} ) where { Tety<:AbstractEvotype,Tenv<:AbstractTabularEnvironment,TG<:AbstractGenotype }
	N <= length(aGty) || throw(DimensionMismatch("Inconsistent Dimensions: N must be > length(aGty)"))
	# @threads
	for i in 1:N fitness!(aGty[i],env) end
	return Population{Tety,Tenv,Vector{TG}}( Int32[N,N,length(aGty)],ety,env,aGty )
end

function fitness!(gty::AbstractGenotype,env::AbstractTabularEnvironment)
	gty.aF .= [ env.fTb[env.e[1]][gty.G[1]], exp(env.fTb[env.e[1]][gty.G[1]]*env.aSelCoef[1]), exp(env.fTb[env.e[1]][gty.G[1]]*env.aSelCoef[2]) ]
end

# -----------
#  MUTATION

function mutation!(gty::TabularGenotype,ety::TabularEvotype,env::AbstractTabularEnvironment)
	aAdjV, aMutProb = neighbors(ety.G, gty.G[1])
	aMutProb .= aMutProb / ( 1.0 + ety.minRepCoef + ety.pRepFactor[1]*gty.aF[2] )

	r::Float64 = rand(THREADRNG[threadid()])
	i = draw(aMutProb,r)

	if i > 0
		gty.G[1] = aAdjV[i]
		fitness!(gty,env)
		return Int32(1)
	else
		return Int32(0)
	end
end

function mutation!(gty::TabularGenotype,ety::VaryingTabularEvotypeI,env::AbstractTabularEnvironment)
	aAdjV, aMutProb = neighbors(ety.G, gty.G[1])
	sqrtρhat = sqrt( 1.0 + ety.minRepCoef + ety.pRepFactor[1]*gty.aF[2] )
	aMutProb .= ( aMutProb .* [ i == gty.G[2] ? ety.λM : 1.0 for i in eachindex(aMutProb) ] ) / sqrtρhat

	Nmut::Int32 = 0

	# site mutation
	r::Float64 = rand(THREADRNG[threadid()])
	i = draw(aMutProb,r)
	if i > 0
		gty.G[1] = aAdjV[i]
		fitness!(gty,env)
		Nmut += 1
	end

	# mutation boost direction
	r = rand(THREADRNG[threadid()])
	if r <= ety.mutCoef / sqrtρhat
		gty.G[2] = typeof(gty.G[2])(rand(0:4))		# <--- you are assuming the square grid!
		Nmut += 1
	end

	return Nmut
end

function evolveEnvironment!(env::VaryingTabularEnvironment,aGty::Vector{<:AbstractGenotype})
	d = Categorical( env.W[:,env.e[1]] )
	i = rand(THREADRNG[threadid()], d)

	if env.e[1] != i
		env.e[1] = i

		for gty in aGty
			fitness!(gty,env)
		end
	end
end

function gmsEvoStep!(pop::Population{<:AbstractEvotype,<:VaryingTabularEnvironment,<:Vector{<:AbstractGenotype}}, elite::Bool=false)
	evolveEnvironment!(pop.env,pop.aGty)
	gf = replication!(pop)
	Nmut = mutation!(pop)
	selection!(pop,elite)

	return gf, Nmut
end

function gmsEvoStep!(pop::Population{<:AbstractEvotype,<:VaryingTabularEnvironment,<:Vector{<:AbstractGenotype}}, ancestry)
	evolveEnvironment!(pop.env,pop.aGty)
	gf = replication!(pop, ancestry)
	Nmut = mutation!(pop)
	iGtySelected = selection!(pop)

	return gf, Nmut, ancestry[iGtySelected]
end

coord(gty::TabularGenotype, cardGty1::Integer) = gty.G[1] + cardGty1 * gty.G[2]

function evolve!(pop::AbstractPopulation,traj::TrajectoryData)
	for gen in 1:traj.NgenRelax
		gmsEvoStep!(pop)
	end

	ancestry = [ i for i in 1:length(pop.aGty) ]

	for gen in 1:traj.NgenSample
		aGtyPast = deepcopy(pop.aGty)
		envStatePast = pop.env.e[1]

		traj.growthFactor[gen], traj.mutationFactor[gen], iGtyPast = gmsEvoStep!(pop, ancestry)
		traj.envState[gen] = pop.env.e[1]
		traj.avePerformance[gen] = mean( [pop.aGty[i].aF[1] for i in 1:pop.pN[2]] )

		traj.jointProb[ coord(pop.aGty[1], pop.ety.G.Nv), traj.envState[gen], coord(aGtyPast[iGtyPast], pop.ety.G.Nv), envStatePast ] += 1
	end
end

function generateVaryingTabularSystemsITrajectories(
		repFactor, minRepCoef, mutCoef, λM, G::AbstractGraph, aFtb, repStrength, selStrength, transitionMatrix, Npop::Integer, NgenRelax::Integer, NgenSample::Integer
	)

	ety = mEvoTypes.VaryingTabularEvotypeI([Float64(repFactor)], Float64(minRepCoef), Float64(mutCoef), Float64(λM), G)
	env = mEvoTypes.VaryingTabularEnvironment(aFtb, [repStrength,selStrength], transitionMatrix)
	aGty = [ mEvoTypes.TabularGenotype( Int32[rand( 1:G.Nv ),0] ) for i in 1:Npop ];
	pop = mEvoTypes.init_Population( Int32(Npop), ety, env, aGty )

	#  assuming that you deal with a square grid
	traj = TrajectoryData( NgenRelax, NgenSample, 5G.Nv, length(aFtb) )
	evolve!(pop, traj)

	return traj
end

export generateVaryingTabularSystemsITrajectories
