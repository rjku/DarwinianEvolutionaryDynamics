
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
	NgenSample::Int32

	avePerformance::Array{Float64,1}
	growthFactor::Array{Float64,1}
	mutationFactor::Array{Float64,1}

	pop::Population
	pGty::Vector{TabularGenotype{Vector{Int32}}}
	fluxMatrix::Matrix{Int32}

	envState::Array{Int32,1}
end
	TrajectoryData(pop::Population,NgenRelax::Integer,NgenSample::Integer,dimG::Integer) = TrajectoryData( Int32(NgenRelax), Int32(NgenSample),
		Array{Float64}(undef,NgenSample), Array{Float64}(undef,NgenSample), Array{Float64}(undef,NgenSample),
		pop, Array{TabularGenotype{Vector{Int32}}}(undef,1), zeros(Int32,dimG,dimG), Array{Int32}(undef,NgenSample)
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

function gmsEvoStep!(pop::Population{<:AbstractEvotype,<:VaryingTabularEnvironment,<:Vector{<:AbstractGenotype}},elite::Bool=false)
	evolveEnvironment!(pop.env,pop.aGty)
	gf = replication!(pop)
	Nmut = mutation!(pop)
	selection!(pop,elite)

	return gf, Nmut #, ancestry
end

function evolve!(traj::TrajectoryData)
	for gen in 1:traj.NgenRelax
		gmsEvoStep!(traj.pop)
	end
	for gen in 1:traj.NgenSample
		# copy initial population

		traj.growthFactor[gen], traj.mutationFactor[gen] = gmsEvoStep!(traj.pop)
		traj.envState[gen] = traj.pop.env.e[1]
		traj.avePerformance[gen] = mean( [traj.pop.aGty[i].aF[1] for i in 1:traj.pop.pN[2]] )

		# feed J with one transition
	end
end

function generateVaryingTabularSystemsITrajectories(
		repFactor,minRepCoef,mutCoef,λM,G::AbstractGraph,
		aFtb,repStrength,selStrength,transitionMatrix,
		Npop::Integer,NgenRelax::Integer,NgenSample::Integer
	)

	ety = mEvoTypes.VaryingTabularEvotypeI([Float64(repFactor)],Float64(minRepCoef),Float64(mutCoef),Float64(λM),G)
	env = mEvoTypes.VaryingTabularEnvironment(aFtb,[repStrength,selStrength],transitionMatrix)
	aGty = [ mEvoTypes.TabularGenotype( Int32[rand( 1:G.Nv ),0] ) for i in 1:Npop ];

	#  assuming that you deal with a square grid
	traj = TrajectoryData( mEvoTypes.init_Population( Int32(Npop), ety, env, aGty ), NgenRelax, NgenSample, 5G.Nv )
	evolve!(traj)

	# traj.pGty[1] = traj.pop.aGty[rand(THREADRNG[threadid()],1:traj.pop.pN[2])]

	return traj
end

export generateVaryingTabularSystemsITrajectories
