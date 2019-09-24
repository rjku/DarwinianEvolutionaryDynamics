
module mEvoTypes

# *******************
# ABSTRACT TYPES
# *******************

# ===================
# simulation types
abstract type atMonteCarloPrm end

export atMonteCarloPrm

# ===================
# environment types
abstract type atEnvironment end

# computational: input--ideal-output | T relates to the representation ...
abstract type atCompEnv <: atEnvironment end

export atEnvironment, atCompEnv

# ===================
# population dynamics types
abstract type atPopulation end

abstract type atEvotype end
abstract type atMetaGenotype end
abstract type atGenotype end

abstract type atPhenotype end		# T relates to the representation of the phenotypic variables

abstract type atIsingMetaGty <: atMetaGenotype end
abstract type atVecGty <: atGenotype end

export atPopulation, atEvotype, atIsingMetaGty, atMetaGenotype, atGenotype, atPhenotype, atVecGty


# *******************
# CONCRETE TYPES
# *******************

# ===================
# type: discrete time monte carlo parameters
struct tDTMCprm <: atMonteCarloPrm
	Nsmpl::Int32							# number of Samples
	Nmcsps::Int32							# number of Monte Carlo Steps per Sample
	Ntrials::Int32							# number of trials for fitness evaluation
end

export tDTMCprm

# ===================
# environment types
# ===================
# type: computational environment
struct tCompEnv{T<:AbstractArray} <: atCompEnv
	idealInputOutput::Array{T,1}
	selFactor::Float64
end

# ===================
# type: trivial environment
struct tTrivialEnv <: atEnvironment end

export tTrivialEnv

# ===================
# population dynamics types
# ===================
# type: evotype
struct tEty{Tx<:Number} <: atEvotype
	repRate::Float64
	mutRate::Float64
	ΔtOffset::Float64
	pRepFactor::Vector{Float64}
	pMutFactor::Vector{Float64}
	Xvar::Tx
end

# ===================
# type: vectorial genotype
struct tVecGty{TpMGty<:Vector{<:atMetaGenotype},Tx<:AbstractVector} <: atVecGty
	pMetaGty::TpMGty
	X::Tx
	pF::Array{Float64,1}
end

# constructor: undefined fitness
tVecGty(pMGty::TpMGty,X::TX) where {TpMGty<:Vector{<:atMetaGenotype},TX<:AbstractVector} =
	tVecGty{TpMGty,TX}(pMGty,X,[0.0])

# ===================
# type: population
struct tLivingPop{Tety<:atEvotype,Tenv<:atEnvironment,TaMGty<:Vector{<:atMetaGenotype},TaGty<:Vector{<:atGenotype}} <: atPopulation
	pN::Vector{Int32}		# population number: effective population, fixed population value, array size
	ety::Tety
	env::Tenv
	aMetaGty::TaMGty
	aGty::TaGty
end

# ===================
# type: evolutionary dynamics data
struct tEvoData
	Ngen::Int32
	growthFactor::Array{Float64,1}
	aveFitness::Array{Float64,1}
	mutationFactor::Array{Float64,1}

	pAveFt::Vector{Float64}
	aveFtinc::Float64
	aLivingPop::Vector{tLivingPop}
end

# initializer constructor
tEvoData(Ngen::Int32,aveFt::Float64,aveFtinc::Float64) =
	tEvoData(Ngen,Array{Float64}(undef,Ngen),Array{Float64}(undef,Ngen),Array{Float64}(undef,Ngen),[aveFt],aveFtinc,Array{tLivingPop}(undef,0))

export tCompEnv, tEty, tVecGty, tEvoData, tLivingPop

# ===================
# type: ising signal transduction
struct tIsingSigTransMGty{Tprm<:atMonteCarloPrm} <: atIsingMetaGty
	L::Int32				# system size
	dX::Int32				# interaction matrix number of entries

	β::Float64 				# inverse temperature [ critical inverse temperature ≃ 1/2.3 ≃ 0.43 ]
	he::Float64 			# global external field H

	L2::Int32;	halfL::Int32
	li::Int32;	li2::Int32			# input/redout region size set to one tenth of the system size
	βhe::Float64

	# definition: useful nearest neighbour coordinate vectors
	jp::Array{Int32,1}; jm::Array{Int32,1}
	Jpi::Array{Int32,2}; Jmi::Array{Int32,2}; Jpj::Array{Int32,2}; Jmj::Array{Int32,2}

	prms::Tprm
end

# constructor: ising signal transduction
tIsingSigTransMGty(L::Int32, β::Real, he::Real, prms::Tprm) where {Tprm<:atMonteCarloPrm} =
	tIsingSigTransMGty{Tprm}(
	L, 2L^2, β, he, L^2, L÷2, L÷10+1, (L÷10+1)^2, β*he,
	Int32[i%L+1 for i in 1:L],									# i -> i+1  with periodic boundary conditions
	Int32[(L-2+i)%L+1 for i in 1:L],							# i -> i-1  with periodic boundary conditions
	Int32[i+2(j-1)*L for i in 1:L, j in 1:L],
	Int32[(L-2+i)%L+1+2(j-1)*L for i in 1:L, j in 1:L],
	Int32[i+(2j-1)*L for i in 1:L, j in 1:L],
	Int32[i+(2(L-2+j)%L+1)*L for i in 1:L, j in 1:L],
	prms
)

export tIsingSigTransMGty

# *******************
# TRIVIAL STUFF
# *******************

struct tTrivialEty <: atEvotype end

export tTrivialEty

end
