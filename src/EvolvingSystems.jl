
abstract type atIsingMetaGty <: atMetaGenotype end
abstract type atChannelMetaGty <: atMetaGenotype end

abstract type atMonteCarloPrm end

export atIsingMetaGty, atChannelMetaGty, atMonteCarloPrm


# ==============
#   SIMULATION
# ==============

# type: discrete time monte carlo parameters
struct tDTMCprm <: atMonteCarloPrm
	Nsmpl::Int32							# number of Samples
	Nmcsps::Int32							# number of Monte Carlo Steps per Sample
	Ntrials::Int32							# number of trials for fitness evaluation
end

export tDTMCprm

# =================
#   ENVIRONMENTAL
# =================

# type: computational environment
struct tDctEnv{T<:Vector{<:Dict{<:AbstractArray,<:Real}}} <: AbstractEnvironment
	fTb::T
	aSelCoef::Vector{Float64}
end

tDctEnv(aGty,fFn,aSelCoef::Vector{Float64}) = tDctEnv([ Dict( e => fFn(e) for e in aGty ) ],aSelCoef)
tDctEnv(aGty,afFn::AbstractArray,aSelCoef::Vector{Float64}) = tDctEnv([ Dict( e => fFn(e) for e in aGty ) for fFn in afFn ],aSelCoef)

# type: computational environment
struct tCompEnv{T<:Vector{<:AbstractVector}} <: atCompEnv
	IOidl::T
	aSelCoef::Vector{Float64}
end

# type: trivial environment
struct tTrivialEnv <: AbstractEnvironment end

export tCompEnv, tTrivialEnv

# ==============
#   POPULATION
# ==============

# -------------------
#   GENOTYPES
# -------------------

# type: genotype with additive genotypic variations
struct tAddGty{Tmgty<:atMetaGenotype,Tpmgty<:Vector{Tmgty},Tg<:Number,TG<:AbstractArray{Tg}} <: atSystemGty{Tmgty}
	pMetaGty::Tpmgty
	G::TG
	Δg::Tg
	aF::Vector{Float64}		# Fitness Indicators: [ loss, replication rate function, selection function ]
end

# constructor: undefined fitness
tAddGty(pMetaGty::Vector{Tmgty},G::Vector{Tg},Δg::Tg) where {Tmgty<:atMetaGenotype,Tg<:Number} =
	tAddGty{Tmgty,Vector{Tmgty},Tg}(pMetaGty,G,Δg,[0.0,0.0,0.0])

# type: genotype with multiplicative genotypic variations
struct tMltGty{Tmgty<:atMetaGenotype,Tpmgty<:Vector{Tmgty},Tg<:Number,TG<:AbstractArray{Tg}} <: atSystemGty{Tmgty}
	pMetaGty::Tpmgty
	G::TG
	δg::Tg
	aF::Vector{Float64}		# Fitness Indicators: [ loss, replication rate function, selection function ]
end

# constructor: undefined fitness
tMltGty(pMetaGty::Vector{Tmgty},G::Vector{Tg},δg::Tg) where {Tmgty<:atMetaGenotype,Tg<:Number} =
	tMltGty{Tmgty,Vector{Tmgty},Tg}(pMetaGty,G,δg,[0.0,0.0,0.0])

# type: genotype with alphabetical/discrete (not unbounded) genotypic variables
struct tAlphaGty{Tmgty<:atMetaGenotype,Tpmgty<:Vector{Tmgty},TG<:AbstractArray} <: atSystemGty{Tmgty}
	pMetaGty::Tpmgty
	G::TG
	g::TG					# Symbols/Genotypic Values array
	aF::Vector{Float64}		# Fitness Indicators: [ loss, replication rate function, selection function ]
end

# constructor: undefined fitness
tAlphaGty(pMetaGty::Vector{Tmgty},G::TG,g::TG) where {Tmgty<:atMetaGenotype,TG<:AbstractArray} =
	tAlphaGty{Tmgty,Vector{Tmgty},TG}(pMetaGty,G,g,[0.0,0.0,0.0])

Base.copy(gty::tAlphaGty) = tAlphaGty( gty.pMetaGty, deepcopy(gty.G), gty.g, deepcopy(gty.aF) )

# type: genotype with binary (not unbounded) genotypic variables
struct tBinaryGty{Tmgty<:atMetaGenotype,Tpmgty<:Vector{Tmgty},TG<:AbstractArray} <: atSystemGty{Tmgty}
	pMetaGty::Tpmgty
	G::TG
	aF::Vector{Float64}		# Fitness Indicators: [ loss, replication rate function, selection function ]
end

# constructor: undefined fitness
tBinaryGty(pMetaGty::Vector{Tmgty},G::TG) where {Tmgty<:atMetaGenotype,TG<:AbstractArray} =
	tBinaryGty{Tmgty,Vector{Tmgty},TG}(pMetaGty,G,[0.0,0.0,0.0])

Base.copy(gty::tBinaryGty) = tBinaryGty( gty.pMetaGty, deepcopy(gty.G), deepcopy(gty.aF) )

# type: genotype with continuous bounded genotypic variables
struct tCntGty{Tmgty<:atMetaGenotype,Tpmgty<:Vector{Tmgty},TG<:AbstractArray} <: atSystemGty{Tmgty}
	pMetaGty::Tpmgty
	G::TG
	gbounds::TG
	aF::Vector{Float64}		# Fitness Indicators: [ loss, replication rate function, selection function ]
end

# constructor: undefined fitness
tCntGty(pMetaGty::Vector{Tmgty},G::TG,gbounds::TG) where {Tmgty<:atMetaGenotype,TG<:AbstractArray} =
	tCntGty{Tmgty,Vector{Tmgty},TG}(pMetaGty,G,gbounds,[0.0,0.0,0.0])

Base.copy(gty::tCntGty) = tCntGty( gty.pMetaGty, deepcopy(gty.G), gty.gbounds, deepcopy(gty.aF))

# type: genotype with continuous bounded genotypic variables and genetic assimilation mechanism
struct tCntGenAssGty{Tmgty<:atMetaGenotype,Tpmgty<:Vector{Tmgty},TG<:AbstractArray} <: atSystemGty{Tmgty}
	pMetaGty::Tpmgty
	G::TG
	gbounds::TG
	aPmut::Vector{Float64}
	aF::Vector{Float64}		# Fitness Indicators: [ loss, replication rate function, selection function ]
end

# constructor: undefined fitness
tCntGenAssGty(pMetaGty::Vector{Tmgty},G::TG,gbounds::TG,aPmut::Vector{Float64}) where {Tmgty<:atMetaGenotype,TG<:AbstractArray} =
	tCntGenAssGty{Tmgty,Vector{Tmgty},TG}(pMetaGty,G,gbounds,aPmut,[0.0,0.0,0.0])

function Base.copy!(gtyDst::Tgty,gtySrc::Tgty) where {Tgty<:AbstractGenotype}
	gtyDst.G .= gtySrc.G
	gtyDst.aF .= gtySrc.aF
end

export tAddGty, tMltGty, tAlphaGty, tCntGty

# -------------------
#   EVOTYPES
# -------------------

# type: discretized time evotype
struct tDscdtEty <: atGtyMutEty
	repRate::Float64
	mutRate::Float64
	ΔtOffset::Float64
	pRepFactor::Vector{Float64}
	pMutFactor::Vector{Float64}
end

# costructor. tDscdtEty
tDscdtEty(repRate::Float64,mutRate::Float64,ΔtOffset::Float64,gty::AbstractGenotype) =
	tDscdtEty( repRate,mutRate,ΔtOffset,[repRate/(2gty.pMetaGty[1].dG*mutRate+ΔtOffset)],[mutRate/(2gty.pMetaGty[1].dG*mutRate+ΔtOffset)] )

# costructor. tDscdtEty
tDscdtEty(repRate::Float64,mutRate::Float64,ΔtOffset::Float64,gty::tAlphaGty) =
	tDscdtEty( repRate,mutRate,ΔtOffset,[repRate/(length(gty.g)*gty.pMetaGty[1].dG*mutRate+ΔtOffset)],
	[mutRate/(length(gty.g)*gty.pMetaGty[1].dG*mutRate+ΔtOffset)] )

# type: evotype
struct tPntMutEty <: atPntMutEty
	pRepFactor::Vector{Float64}
	pPntMutFactor::Vector{Float64}
	aSize2PMF::Dict{Int32,Vector{Float64}}			# genotypic dimension => [ binomial PMF of mutation ]
	minRepCoef::Float64
end

# costructor. tPointMutationEty
tPntMutEty(repFactor::Float64,pntMutFactor::Float64,adG::Vector{<:Integer},NmutMax::Integer) = tPntMutEty([repFactor], [pntMutFactor],
	Dict( dG => [ binomial(BigInt(dG),n)*(pntMutFactor^n)*(1-pntMutFactor)^(dG-n) for n in 1:NmutMax ] for dG in adG), 0.0)

# costructor. tPointMutationEty
tPntMutEty(repFactor::Float64,pntMutFactor::Float64,adG::Vector{<:Integer},NmutMax::Integer,minRepCoef::Float64) = tPntMutEty([repFactor],[pntMutFactor],
	Dict( dG => [ binomial(BigInt(dG),n)*(pntMutFactor^n)*(1-pntMutFactor)^(dG-n) for n in 1:NmutMax ] for dG in adG), minRepCoef)

export tDscdtEty, tPntMutEty

# -------------------
#   METAGENOTYPES
# -------------------

# type: ising signal transduction
struct tIsingSigTransMetaGty{Tprm<:atMonteCarloPrm} <: atIsingMetaGty
	L::Int32				# system size
	dG::Int32				# interaction matrix number of entries

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
tIsingSigTransMetaGty(L::Int32, β::Real, he::Real, prms::Tprm) where {Tprm<:atMonteCarloPrm} =
	tIsingSigTransMetaGty{Tprm}(
	L, 2L^2, β, he, L^2, L÷2, L÷10+1, (L÷10+1)^2, β*he,
	Int32[i%L+1 for i in 1:L],									# i -> i+1  with periodic boundary conditions
	Int32[(L-2+i)%L+1 for i in 1:L],							# i -> i-1  with periodic boundary conditions
	Int32[i+2(j-1)*L for i in 1:L, j in 1:L],
	Int32[(L-2+i)%L+1+2(j-1)*L for i in 1:L, j in 1:L],
	Int32[i+(2j-1)*L for i in 1:L, j in 1:L],
	Int32[i+(2(L-2+j)%L+1)*L for i in 1:L, j in 1:L],
	prms
)

# type: disordered channel metagenotype <: atChannelMetaGty
struct tDisChnMetaGty <: atChannelMetaGty
	L::Int32					# system size
	L2::Int32;	L2mL::Int32
	dG::Int32;	Nvbl::Int32		# genotype length and viable transitions
	p::Float64 					# viability probability
	V::Vector{Vector{Int32}}	# [ t(e), o(e) ] vector
	kout::Float64				# output rate constant
end

export tIsingSigTransMetaGty, tDisChnMetaGty

# -------------------
#   POPULATION and DATA
# -------------------

# type: evolutionary dynamics data
struct tEvoData <: AbstractEvolutionData
	Ngen::Int32
	avePerformance::Array{Float64,1}
	growthFactor::Array{Float64,1}
	mutationFactor::Array{Float64,1}

	pAveFt::Vector{Float64}
	aveFtinc::Float64
	pMinF::Vector{Float64}

	aGen::Vector{Int32}
	aEvoPop::Vector{tEvoPop}
end

# initializer constructors
tEvoData(Ngen::Int32,aveFt::Float64,aveFtinc::Float64,minF::Float64) = tEvoData(
	Ngen, Array{Float64}(undef,Ngen), Array{Float64}(undef,Ngen), Array{Float64}(undef,Ngen),
	[aveFt], aveFtinc, [minF],
	Array{Int32}(undef,0), Array{tEvoPop}(undef,0)
)

tEvoData(Ngen::Int32,aveFtinc::Float64,minF::Float64) = tEvoData(
	Ngen, Array{Float64}(undef,Ngen), Array{Float64}(undef,Ngen), Array{Float64}(undef,Ngen),
	[0.0], aveFtinc, [minF],
	Array{Int32}(undef,0), Array{tEvoPop}(undef,0)
)

# type: flux pattern
struct tFluxPattern
	Nq::Int32
	V::Vector{Vector{Int32}}	# [ t(e), o(e) ] vector
	jave::Vector{Vector{Float64}}
	jcov::Vector{Array{Float64,2}}
	jcor::Array{Float64,2}
end

function tFluxPattern(env::atCompEnv,gty::atSystemGty{<:atChannelMetaGty})
	Nq = gty.pMetaGty[1].dG + 2gty.pMetaGty[1].L
	tFluxPattern( Nq,
		[ Vector{Int32}(undef, Nq), Vector{Int32}(undef, Nq) ],
		[ Vector{Float64}(undef, Nq) for i in eachindex(env.IOidl) ],
		[ Array{Float64}(undef, Nq, Nq) for i in eachindex(env.IOidl) ],
		Array{Float64}(undef, Nq, Nq)
	)
end

struct tCrntPattern
	Nq::Int32
	aV::Vector{Vector{Vector{Int32}}}	# [ t(e), o(e) ] vectors for each i/o
	Jave::Vector{Vector{Float64}}
end

function tCrntPattern(env::atCompEnv,gty::atSystemGty{<:atChannelMetaGty})
	Nq = Int32(gty.pMetaGty[1].dG/2 + 2gty.pMetaGty[1].L)
	tCrntPattern( Nq,
		[ [ Vector{Int32}(undef, Nq), Vector{Int32}(undef, Nq) ] for i in eachindex(env.IOidl) ],
		[ Vector{Float64}(undef, Nq) for i in eachindex(env.IOidl) ]
	)
end

# type: statistics: average, covariance, and correlations
struct tStat{Tave<:Vector{<:Number},Tcov<:Array{<:Number}}
	ave::Tave
	cov::Tcov
	cor::Tcov
end

tStat(gty::AbstractGenotype) = tStat(
	Array{typeof(gty.G[1])}(undef,length(gty)),
	Array{typeof(gty.G[1])}(undef,length(gty),length(gty)),
	Array{typeof(gty.G[1])}(undef,length(gty),length(gty))
)

tStat(Nv::Integer) = tStat( Array{Float64}(undef,Nv), Array{Float64}(undef,Nv,Nv), Array{Float64}(undef,Nv,Nv) )

export EvoData, tEvoData, tFluxPattern, tCrntPattern, tStat

# ==========
#  METHODS

function mutation!(gty::tAddGty,i::Integer)
	gty.G[ i % length(gty) + 1 ] += i % 2 == 0 ? gty.Δg : -gty.Δg
end

function mutation!(gty::tMltGty,i::Integer)
	gty.G[ i % length(gty) + 1 ] *= i % 2 == 0 ? gty.δg : 1.0/gty.δg
end

function mutation!(gty::tAlphaGty,i::Integer)
	gty.G[i] = rand(THREADRNG[threadid()], filter( e -> e != gty.G[i], gty.g ))
end

function mutation!(gty::tCntGty,i::Integer)
	gty.G[i] = gty.gbounds[1] + rand(THREADRNG[threadid()])*(gty.gbounds[2] - gty.gbounds[1])
end
