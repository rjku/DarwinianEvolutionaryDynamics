
module Types

abstract type AbstractEvotype end							# abstract evotypes are endowed with a vector "replicationType", and a "selectionType"
abstract type AbstractEnvironment end
abstract type AbstractGenotype end							# abstract genotypes are endowed with a vector "genome" and a fitness value array "aFitness"
	Base.length(gty::AbstractGenotype) = length(gty.genome)

	fitness(gty::AbstractGenotype) = gty.aFitness[1]

abstract type AbstractPopulation end
abstract type AbstractEvolutionData end

export AbstractEvotype, AbstractEnvironment, AbstractGenotype, AbstractPopulation, AbstractEvolutionData, fitness
# --------------------
# Environmental Traits

abstract type IsVaryingEnvironment end
struct StationaryEnvironment <: IsVaryingEnvironment end
struct MarkovianEnvironment <: IsVaryingEnvironment end				# characterized by envState::Vector{Int32}, and transMtx::Matrix{Float64}

IsVaryingEnvironment(::Type{<:AbstractEnvironment}) = StationaryEnvironment()

envState(env::Tenv) where Tenv <: AbstractEnvironment = _envState(IsVaryingEnvironment(Tenv), env)

	_envState(::IsVaryingEnvironment, env) = error("environmental state not defined")

	_envState(::MarkovianEnvironment, env) = env.envState[1]

export IsVaryingEnvironment, StationaryEnvironment, MarkovianEnvironment, envState
# -----------
# Replication

abstract type ReplicationType end
struct FitnessReplication <: ReplicationType
	repFactor::Float64
end
struct NeutralReplication <: ReplicationType
	repCoef::Int32
end
struct WithoutReplication <: ReplicationType end

# default behavior
ReplicationType(::Type) = NeutralReplication()

export ReplicationType, FitnessReplication, NeutralReplication, WithoutReplication
# --------
# Mutation

abstract type MutationType end
struct StandardMutation <: MutationType end
struct RescaledMutation <: MutationType end

# default behavior
MutationType(::Type) = StandardMutation()

export MutationType, StandardMutation, RescaledMutation
# ---------
# Selection

abstract type SelectionType end
struct FitnessSelection <: SelectionType end
struct NeutralSelection <: SelectionType end
struct ElitismSelection <: SelectionType end

# default behavior
SelectionType(::Type) = FitnessSelection()

export SelectionType, FitnessSelection, NeutralSelection, ElitismSelection
# --------------------
# Basic Concrete Types

struct Evotype{Trep<:ReplicationType,Tmut<:MutationType,Tsel<:SelectionType} <: AbstractEvotype
	repType::Trep
	mutType::Tmut
	selType::Tsel
end

struct Genotype{TG<:AbstractArray} <: AbstractGenotype
	genome::TG
	aFitness::Vector{Float64}
end
	Genotype(genome::AbstractArray) = Genotype(genome, [0.0, 0.0])

	Base.copy(gty::Genotype) = Genotype(deepcopy(gty.genome), deepcopy(gty.aFitness))

struct Population{Tety<:AbstractEvotype, Tenv<:AbstractEnvironment, TaG<:Vector{<:AbstractGenotype}} <: AbstractPopulation
	pN::Vector{Int32}		# population number: effective population, fixed population value
	ety::Tety
	env::Tenv
	aGty::TaG
end

export Evotypes, Genotype, Population

end  # module Types

# ===================
# computational: input--ideal-output | T relates to the representation ...
# abstract type atMetaGenotype end
# abstract type atCompEnv <: AbstractEnvironment end
#
# export AbstractEnvironment, atCompEnv
#
# abstract type atGtyMutEty <: AbstractEvotype end
# abstract type atPntMutEty <: AbstractEvotype end
#
# abstract type atSystemGty{Tmgty} <: AbstractGenotype end		# genotypes endowed with a pointer to a atMetaGenotype: pMetaGty
#
# export AbstractPopulation, AbstractEvotype, atGtyMutEty, atPntMutEty, atMetaGenotype, atIsingMetaGty, atChannelMetaGty
# export AbstractGenotype, atPhenotype, atSystemGty
#
# # type: population
# struct tEvoPop{Tety<:AbstractEvotype,Tenv<:AbstractEnvironment,Tamgty<:Vector{<:atMetaGenotype},TaG<:Vector{<:AbstractGenotype}} <: AbstractPopulation
# 	pN::Vector{Int32}		# population number: effective population, fixed population value, array size
# 	ety::Tety
# 	env::Tenv
# 	aMetaGty::Tamgty
# 	aGty::TaG
# end
#
# export tEvoPop
