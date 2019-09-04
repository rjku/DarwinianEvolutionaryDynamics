
module mEvoTypes

# *******************
# ABSTRACT TYPES
# *******************

# ===================
# system types
abstract type atSystem end
abstract type atThermoSystem <: atSystem end

abstract type atMonteCarloPrm end

export atSystem, atThermoSystem, atMonteCarloPrm

# ===================
# environment types
abstract type atEnvironment end

# computational: input--ideal-output | T relates to the representation ...
abstract type atCompEnv{T} <: atEnvironment end

export atEnvironment, atCompEnv

# ===================
# population dynamics types
abstract type atPopulation end
abstract type atEvotype end
abstract type atGenotype end			# T relates to the representation of the genotypic variables
abstract type atPhenotype{T} end		# T relates to the representation of the phenotypic variables

abstract type at1dGty <: atGenotype end

export atPopulation, atEvotype, atGenotype, atPhenotype, at1dGty


# *******************
# CONCRETE TYPES
# *******************

# ===================
# type: simulation parameters
struct tDTMCprm <: atMonteCarloPrm
	Nsmpl::Int32							# number of Samples
	Nmcsps::Int32							# number of Monte Carlo Steps per Sample
end

export tDTMCprm

# ===================
# type: computational environment
# struct tCompEnv{T} <: atCompEnv{T}
# 	idealInputOutput::Array{Array{T},1}
# 	selFactor::Float64
# end
# to improve:
struct tCompEnv{T<:AbstractArray} <: atCompEnv{T}
	idealInputOutput::Array{T,1}
	selFactor::Float64
end

# ===================
# type: trivial environment
struct tTrivialEnv <: atEnvironment end

# ===================
# type: vectorial genotype
struct tVecGty{T<:AbstractVector} <: at1dGty
	pdX::Array{Int32,1}
	X::T
	pF::Array{Float64,1}
end

# constructor: undefined fitness
tVecGty{T}(X::AbstractVector) where {T<:AbstractVector} = tVecGty{T}([length(X)],X,[0.0])

# ===================
# type: evolutionary dynamics data
struct tEvoData
	Ngen::Int32
	growthFactor::Array{Float64,1}
	aveFitness::Array{Float64,1}
	mutationNumber::Array{Float64,1}
end

# simplified constructor
tEvoData(Ngen::Int32) = tEvoData(Ngen,Array{Float64}(undef,Ngen),Array{Float64}(undef,Ngen),Array{Float64}(undef,Ngen))

# ===================
# type: population
struct tLivingPop{T, Tevo<:atEvotype, Tenv<:atEnvironment, TaGty<:Array{<:atGenotype,1}} <: atPopulation
	pN::Array{Int32,1}		# population number: effective population, fixed population value, array size
	ety::Tevo
	env::Tenv
	aGty::TaGty		# genotypes

	repFactor::Float64
	mutFactor::Float64
	Xvar::T
end

export tCompEnv, tVecGty, tEvoData, tLivingPop

# ===================
# type: ising signal transduction
struct tIsingSigTransEty <: atEvotype
	L::Int32				# system size
	lJij::Int32				# interaction matrix number of entries

	β::Float64 				# inverse temperature [ critical inverse temperature ≃ 1/2.3 ≃ 0.43 ]
	he::Float64 			# global external field H

	L2::Int32;	halfL::Int32
	li::Int32;	li2::Int32			# input/redout region size set to one tenth of the system size
	ℍe::Float64

	# definition: useful nearest neighbour coordinate vectors
	jp::Array{Int32,1}; jm::Array{Int32,1}
	Jpi::Array{Int32,2}; Jmi::Array{Int32,2}; Jpj::Array{Int32,2}; Jmj::Array{Int32,2}
end

# constructor: ising signal transduction
tIsingSigTransEty(L::Int32, β::Real, he::Real) = tIsingSigTransEty(
	L, 2L^2, β, he, L^2, L÷2, L÷10, (L÷10)^2, he*β,
	Int32[i%L+1 for i in 1:L],
	Int32[(L-2+i)%L+1 for i in 1:L],
	Int32[i+2(j-1)*L for i in 1:L, j in 1:L],
	Int32[(L-2+i)%L+1+2(j-1)*L for i in 1:L, j in 1:L],
	Int32[i+(2j-1)*L for i in 1:L, j in 1:L],
	Int32[i+(2(L-2+j)%L+1)*L for i in 1:L, j in 1:L]
)

export tIsingSigTransEty

end
