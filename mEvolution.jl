
module mEvolution
using mAbstractTypes

export tCompEnv, tVecGty, tTupGty, tLivingPop

# ===================
# type: computational environment
struct tCompEnv{T} <: atCompEnv{T}
	idealInputOutput::Array{Array{T},1}
end

# ===================
# type: tuplonic genotype
struct tTupGty{T} <: at1dGty{T}
	X::Array{T,1}
	F::Array{Float64,1}
end

# constructor: undefined fitness
tTupGty{T}(X::Array{T,1}) where {T} = tTupGty{T}(X,Array{Float64}(undef,1))

# ===================
# type: vectorial genotype
mutable struct tVecGty{T} <: at1dGty{T}
	dX::Int32
	X::Array{T,1}
end

# ===================
# type: population
struct tLivingPop{T} <: atPopulation
	pN::Array{Int32,1}		# population number: effective population, fixed population value, array size
	ety::atEvotype
	env::atEnv
	aGty::Array{atGenotype{T},1}		# genotypes
end

# simplified constructor
tLivingPop{T}(N::Int32,ety::atEvotype,env::atEnv,aGty::Array{<:atGenotype{T},1}) where {T} =
	tLivingPop{T}( Int32[N,N,length(aGty)],ety,env,aGty )

end
