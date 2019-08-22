
module mEvoTypes

# *******************
# ABSTRACT TYPES
# *******************

# ===================
# system types
abstract type atSystem end
abstract type atThermoSystem <: atSystem end

export atSystem, atThermoSystem

# ===================
# environment types
abstract type atEnv end
abstract type atCompEnv{T} <: atEnv end		# computational: input--ideal-output | T relates to the representation ...

export atEnv, atCompEnv

# ===================
# population dynamics types
abstract type atPopulation end
abstract type atEvotype end
abstract type atGenotype{T} end			# T relates to the representation of the genotypic variables
abstract type atPhenotype{T} end		# T relates to the representation of the phenotypic variables

abstract type at1dGty{T} <: atGenotype{T} end

export atPopulation, atEvotype, atGenotype, atPhenotype, at1dGty


# *******************
# CONCRETE TYPES
# *******************

# ===================
# type: computational environment
struct tCompEnv{T} <: atCompEnv{T}
	idealInputOutput::Array{Array{T},1}
end

# ===================
# type: trivial environment
struct tTrivialEnv{T} <: atCompEnv{T} end

# ===================
# type: vectorial genotype
struct tVecGty{T} <: at1dGty{T}
	X::Array{T,1}
	pF::Array{Float64,1}
	pdX::Array{Int32,1}
end

# constructor: undefined fitness
tVecGty{T}(X::Array{T,1}) where {T} = tVecGty{T}(X,Array{Float64}(undef,1),[length(X)])
# tVecGty{T}(X::Array{T,1}) where {T} = tVecGty{T}(X,[0.0],[length(X)])

# ===================
# type: evolutionary dynamics data
struct tEvoData
	Ngen::Int32
	growthFactor::Array{Float64,1}
	aveFitness::Array{Float64,1}
end

# simplified constructor
tEvoData(Ngen::Int32) = tEvoData(Ngen, Array{Float64}(undef,Ngen), Array{Float64}(undef,Ngen))

# ===================
# type: population
struct tLivingPop{T} <: atPopulation
	pN::Array{Int32,1}		# population number: effective population, fixed population value, array size
	ety::atEvotype
	env::atEnv
	aGty::Array{atGenotype{T},1}		# genotypes

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
