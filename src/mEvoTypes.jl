
module mEvoTypes

using mUtils, Random
import Distributions.Categorical
import Statistics.mean

# ===================
# environmental types
abstract type AbstractEnvironment end

# computational: input--ideal-output | T relates to the representation ...
abstract type atCompEnv <: AbstractEnvironment end

export AbstractEnvironment, atCompEnv

# ===================
# population dynamics types
abstract type AbstractPopulation end

abstract type AbstractEvotype end
abstract type atMetaGenotype end
abstract type AbstractGenotype end							# abstract genotypes are endowed with a genome vector G

Base.length(gty::AbstractGenotype) = length(gty.G)

abstract type atPhenotype end

abstract type atGtyMutEty <: AbstractEvotype end
abstract type atPntMutEty <: AbstractEvotype end

abstract type atSystemGty{Tmgty} <: AbstractGenotype end		# genotypes endowed with a pointer to a atMetaGenotype: pMetaGty

abstract type AbstractEvolutionData end

export AbstractPopulation, AbstractEvotype, atGtyMutEty, atPntMutEty, atMetaGenotype, atIsingMetaGty, atChannelMetaGty
export AbstractGenotype, atPhenotype, atSystemGty


struct Population{Tety<:AbstractEvotype,Tenv<:AbstractEnvironment,TaG<:Vector{<:AbstractGenotype}} <: AbstractPopulation
	pN::Vector{Int32}		# population number: effective population, fixed population value, array size
	ety::Tety
	env::Tenv
	aGty::TaG
end

# type: population
struct tEvoPop{Tety<:AbstractEvotype,Tenv<:AbstractEnvironment,Tamgty<:Vector{<:atMetaGenotype},TaG<:Vector{<:AbstractGenotype}} <: AbstractPopulation
	pN::Vector{Int32}		# population number: effective population, fixed population value, array size
	ety::Tety
	env::Tenv
	aMetaGty::Tamgty
	aGty::TaG
end

export Population, tEvoPop

include("EvolvingSystems.jl")

include("TabularSystems.jl")

include("EDComponents.jl")

include("EDAnalysis.jl")

end
