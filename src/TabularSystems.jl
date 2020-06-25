
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

struct VarData <: AbstractEvolutionData
	NgenRelax::Int32
	NgenSample::Int32

	avePerformance::Array{Float64,1}
	growthFactor::Array{Float64,1}
	mutationFactor::Array{Float64,1}

	pGty::Vector{TabularGenotype{Vector{Int32}}}
	fluxMatrix::Matrix{Int32}
end
	VarData(NgenRelax::Integer,NgenSample::Integer,dimG::Integer) = VarData( Int32(NgenRelax), Int32(NgenSample), Array{Float64}(undef,Ngen), Array{Float64}(undef,Ngen), Array{Float64}(undef,Ngen),
 		Array{TabularGenotype{Vector{Int32}}}(undef,1), zeros(Int32,dimG,dimG) )

export EvoData, VarData


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
	@threads for i in 1:N fitness!(aGty[i],env) end
	return Population{Tety,Tenv,Vector{TG}}( Int32[N,N,length(aGty)],ety,env,aGty )
end

function fitness!(gty::AbstractGenotype,env::AbstractTabularEnvironment)
	gty.aF .= [ env.fTb[1][gty.G[1]], exp(env.fTb[1][gty.G[1]]*env.aSelCoef[1]), exp(env.fTb[1][gty.G[1]]*env.aSelCoef[2]) ]
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
	ρhat = 1.0 + ety.minRepCoef + ety.pRepFactor[1]*gty.aF[2]
	aMutProb .= ( aMutProb .* [ i == gty.G[2] ? ety.λM : 1.0 for i in eachindex(aMutProb) ] ) / ρhat

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
	if r <= ety.mutCoef    # <----- check for the presence of rhat
		if Nmut == 1
			aAdjV, aMutProb = neighbors(ety.G, gty.G[1])
		end
		i = draw( [ ety.mutCoef/(length(aAdjV)+1) for i in 0:length(aAdjV) ], r ) - 1
		gty.G[2] = typeof(gty.G[2])(i)
		Nmut += 1
	end

	return Nmut
end

function evolveEnvironment!(env::VaryingTabularEnvironment,aGty::Vector{<:AbstractGenotype})
	i = draw( env.W[:,env.e[1]], rand(THREADRNG[threadid()]) )

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

function evolve!(pop::AbstractPopulation,traj::VarData)
	for gen in 1:data.NgenRelax
		gmsEvoStep!(pop)
	end
	for gen in 1:data.NgenSample
		# copy initial population

		traj.avePerformance[gen] = mean( [pop.aGty[i].aF[1] for i in 1:pop.pN[2]] )
		traj.growthFactor[gen], traj.mutationFactor[gen] = gmsEvoStep!(pop)

		# feed J with one transition
	end
end

function generateVaryingTabularSystemsITrajectories(
		repFactor,minRepCoef,mutCoef,λM,G::AbstractGraph,
		aFtb,repStrength,selStrength,transitionMatrix,
		Npop::Integer,
		NgenRelax::Integer,
		NgenSample::Integer
	)

	ety = mEvoTypes.VaryingTabularEvotypeI([Float64(repFactor)],Float64(minRepCoef),Float64(mutCoef),Float64(λM),G)
	env = mEvoTypes.VaryingTabularEnvironment(aFtb,[repStrength,selStrength],transitionMatrix)
	aGty = [ mEvoTypes.TabularGenotype( [rand( 1:G.Nv),0] ) for i in 1:Npop ];
	pop = mEvoTypes.init_Population( Int32(Npop),ety,env,aGty );
	traj = VarData(NgenRelax,NgenSample,5G.Nv)  # <----- here you are assuming that you deal with a square grid

	evolve!(pop,traj)

	traj.pGty[1] = pop.aGty[rand(THREADRNG[threadid()],1:pop.pN[2])]
	return varData
end
