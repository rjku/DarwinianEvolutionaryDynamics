
module mEvoFunc
using mEvoTypes, Base.Threads, Random, Statistics, SparseArrays, Distributions, ForwardDiff, DelimitedFiles, Dates, Distances, ProgressMeter
import Future

# threads random generators initialization
const THREADRNG = let m = MersenneTwister(1)
            [m; accumulate(Future.randjump, fill(big(10)^20, nthreads()-1), init=m)]
        end;

const FITNESSOFFSET, FITNESSTHRESHOLD, BASALFITNESS, MAXSIZE = .001, .99, 3.0, 12

# initializer constructor for discrete time evolution. Same MetaGenotype for all.
function initLivingPop( N::Int32,ety::Tety,env::Tenv,aMGty::Vector{Tmgty},aGty::Vector{Tgty} ) where {
		Tety<:atEvotype,Tenv<:atEnvironment,Tmgty<:atMetaGenotype,Tgty<:atGenotype }
	N <= length(aGty) || throw(DimensionMismatch("Inconsistent Dimensions: N > length(aGty)"))
	@threads for i in 1:N fitness!(env,aGty[i]) end
	return tLivingPop{Tety,Tenv,Vector{Tmgty},Vector{Tgty}}( Int32[N,N,length(aGty)],ety,env,aMGty,aGty )
end

# function. changing rep and mut -factors in tDscdtEty
function set_tDscdtEtyFactors(ety::tDscdtEty,gty::atGenotype)
	ety.pRepFactor[1] = ety.repRate/(2gty.pMetaGty[1].dG*ety.mutRate+ety.ΔtOffset)
	ety.pMutFactor[1] = ety.mutRate/(2gty.pMetaGty[1].dG*ety.mutRate+ety.ΔtOffset)
end

# function. changing rep and mut -factors in tDscdtEty
function set_tDscdtEtyFactors(ety::tDscdtEty,gty::tAlphaGty)
	ety.pRepFactor[1] = ety.repRate/(gty.pdg[1]*gty.pMetaGty[1].dG*ety.mutRate+ety.ΔtOffset)
	ety.pMutFactor[1] = ety.mutRate/(gty.pdg[1]*gty.pMetaGty[1].dG*ety.mutRate+ety.ΔtOffset)
end

export initLivingPop


# *********************************
# | PROPER EVOLUTIONARY FUNCTIONS \
# *********************************

# function. population replication!
function replication!(pop::tLivingPop)
	G = zeros(Int32,pop.pN[2])
	Kr = Vector{Float64}(undef,nthreads())
	ipKr = Vector{Int32}(undef,nthreads())

	@threads for i in 1:pop.pN[2]
		Kr[threadid()] = pop.ety.pRepFactor[1]*pop.aGty[i].pF[1]
		ipKr[threadid()] = trunc(Int32,Kr[threadid()])
		G[i] = rand(THREADRNG[threadid()]) < Kr[threadid()] - ipKr[threadid()] ? ipKr[threadid()] + 1 : ipKr[threadid()]
	end

	for i in 1:pop.pN[2]
		for inew in 1:G[i]
			pop.pN[1] += 1
			if pop.pN[1] <= pop.pN[3]
				pop.aGty[pop.pN[1]] = deepcopy(pop.aGty[i]) 		# <- check whether is necessary to deepcopy
			else
				pop.pN[3] += 1
				push!(pop.aGty,pop.aGty[i])
			end
		end
	end

	return log(pop.pN[1]/pop.pN[2])
end

function blindMutation!(gty::tAddGty,ety::atGtyMutEty)
	Pmut = ety.pMutFactor[1]/(1+ety.pRepFactor[1]*gty[i].pF[1])
	CPmut = 2gty.pMetaGty[1].dG*Pmut
	CPmut <= 1 || throw("cumulative probability of mutation exceeds 1")

	CDFmut::Float64 = Pmut; ig::Int32 = 0; r::Float64 = rand(THREADRNG[threadid()])
	if r <= CPmut
		while CDFmut < r
			CDFmut += Pmut
			ig += 1
		end
		gty.G[ ig % gty.pMetaGty[1].dG + 1 ] += ig % 2 == 0 ? gty.Δg : -gty.Δg
		return Int32(1)
	else
		return Int32(0)
	end
end

function blindMutation!(gty::tMltGty,ety::atGtyMutEty)
	Pmut = ety.pMutFactor[1]/(1+ety.pRepFactor[1]*gty[i].pF[1])
	CPmut = 2gty.pMetaGty[1].dG*Pmut
	CPmut <= 1 || throw("cumulative probability of mutation exceeds 1")

	CDFmut::Float64 = Pmut; ig::Int32 = 0; r::Float64 = rand(THREADRNG[threadid()])
	if r <= CPmut
		while CDFmut < r
			CDFmut += Pmut
			ig += 1
		end
		gty.G[ ig % gty.pMetaGty[1].dG + 1 ] *= ig % 2 == 0 ? gty.δg : 1.0/gty.δg
		return Int32(1)
	else
		return Int32(0)
	end
end

function blindMutation!(gty::tAlphaGty,ety::atGtyMutEty)
	Pmut = ety.pMutFactor[1]/(1+ety.pRepFactor[1]*gty[i].pF[1])
	CPmut = gty.pMetaGty[1].dG*gty.pdg[1]*Pmut
	CPmut <= 1 || throw("cumulative probability of mutation exceeds 1")

	CDFmut::Float64 = Pmut; ig::Int32 = 0; r::Float64 = rand(THREADRNG[threadid()])
	if r <= CPmut
		while CDFmut < r
			CDFmut += Pmut
			ig += 1
		end
		gty.G[ ig % gty.pMetaGty[1].dG + 1 ] = gty.g[ ig % gty.pdg[1] + 1 ]
		return Int32(1)
	else
		return Int32(0)
	end
end

function blindMutation!(gty::tCntGty,ety::atPntMutEty)
	K = 1.0/(1.0 + ety.pRepFactor[1]*gty.pF[1])
	CDFmut = 1.0-K*(1.0-(1.0-ety.pPntMutFactor[1])^gty.pMetaGty[1].dG)		# probability of no mutation
	Nmut::Int32 = 0;	rNmut = rand(THREADRNG[threadid()])

	while CDFmut < rNmut && Nmut < length(ety.aSize2PDF[gty.pMetaGty[1].dG])
		CDFmut += K*ety.aSize2PDF[gty.pMetaGty[1].dG][Nmut+=1]
	end

	aMut = Vector{Int32}(undef, Nmut)
	for nmut in 1:Nmut
		aMut[nmut] = rand(THREADRNG[threadid()], filter( e -> !(e in aMut), 1:gty.pMetaGty[1].dG ))
	end

	for i in aMut
		gty.G[i] = gty.gbounds[1] + rand(THREADRNG[threadid()])*(gty.gbounds[2] - gty.gbounds[1])
	end

	return Nmut
end

# function. effective mutation: blindly mutate the population's genotype according to the effective dynamical mutation rate
function effMutation!(pop::tLivingPop)
	Nmutations = zeros(Int32,nthreads())
	@threads for i in 1:pop.pN[1]
		Nmutations[threadid()] += blindMutation!(pop.aGty[i],pop.ety)
		fitness!(pop.env,pop.aGty[i])
	end
	return sum(Nmutations)/pop.pN[2] 	# normalized number of mutations
end

# function. effective selection: pruning of the population
function effSelection!(pop::tLivingPop, ubermode::Bool)
	popGtyRef::Array{atGenotype,1} = copy(pop.aGty)

	if ubermode
		# survival of the fittest
		survivedGty = sortperm([pop.aGty[i].pF[1] for i in 1:pop.pN[1]],rev=true)
		# sort!(pop.aGty, by= x -> x.pF[1], rev=true)
	else
		# selection with replacement of individuals
		survivedGty = rand(THREADRNG[threadid()],1:pop.pN[1],pop.pN[2])
	end
	for i in 1:pop.pN[2]
		pop.aGty[i] = popGtyRef[survivedGty[i]]
	end

	# population size renormalization
	pop.pN[1] = pop.pN[2]
end

# function. condition for genotypic upgrade
function upgradeCondition(gty::atSystemGty{<:atIsingMetaGty})
	return gty.pF[1] - (gty.pMetaGty[1].halfL - BASALFITNESS) > FITNESSTHRESHOLD && gty.pMetaGty[1].L < MAXSIZE
end

# function. new genotypic variables choice. additive genotype case
function newg!(gty::tAddGty{<:atIsingMetaGty},Gref::Vector)
	# duplication of previous line genotype + some randomness
	for k in 0:1, j in 1:gty.pMetaGty[1].L, i in 1:2
		gty.G[i*(gty.pMetaGty[1].halfL+1)+(j+k*(gty.pMetaGty[1].L+2)-1)*(gty.pMetaGty[1].L+2)] =
			Gref[i*gty.pMetaGty[1].halfL+(j+k*(gty.pMetaGty[1].L)-1)*gty.pMetaGty[1].L] + rand(THREADRNG[threadid()],-1:1)*gty.Δg
	end
	for u in 0:1, k in 0:1, j in 1:2, i in 1:gty.pMetaGty[1].halfL
		gty.G[i+k*(gty.pMetaGty[1].halfL+1)+(j+gty.pMetaGty[1].L+u*(gty.pMetaGty[1].L+2)-1)*(gty.pMetaGty[1].L+2)] =
			Gref[i+k*gty.pMetaGty[1].halfL+(j+gty.pMetaGty[1].L*(u+1)-3)*gty.pMetaGty[1].L] + rand(THREADRNG[threadid()],-1:1)*gty.Δg
	end
end

# function. new genotypic variables choice. alphabetic genotype case
function newg!(gty::tAlphaGty{<:atIsingMetaGty},Gref::Vector)
	# duplication of previous line genotype + some randomness
	for k in 0:1, j in 1:gty.pMetaGty[1].L, i in 1:2
		gty.G[i*(gty.pMetaGty[1].halfL+1)+(j+k*(gty.pMetaGty[1].L+2)-1)*(gty.pMetaGty[1].L+2)] = rand(THREADRNG[threadid()],gty.g)
	end
	for u in 0:1, k in 0:1, j in 1:2, i in 1:gty.pMetaGty[1].halfL
		gty.G[i+k*(gty.pMetaGty[1].halfL+1)+(j+gty.pMetaGty[1].L+u*(gty.pMetaGty[1].L+2)-1)*(gty.pMetaGty[1].L+2)] = rand(THREADRNG[threadid()],gty.g)
	end
	# there are a couple of connection missing in this routine
end

# function. new genotypic variables choice. alphabetic genotype case
function newg!(gty::tCntGty{<:atIsingMetaGty},Gref::Vector)
	# duplication of previous line genotype + some randomness
	for k in 0:1, j in 1:gty.pMetaGty[1].L, i in 1:2
		gty.G[i*(gty.pMetaGty[1].halfL+1)+(j+k*(gty.pMetaGty[1].L+2)-1)*(gty.pMetaGty[1].L+2)] =
			Gref[i*gty.pMetaGty[1].halfL+(j+k*(gty.pMetaGty[1].L)-1)*gty.pMetaGty[1].L]
	end

	for u in 0:1, k in 0:1, j in 1:2, i in 1:gty.pMetaGty[1].halfL
		gty.G[i+k*(gty.pMetaGty[1].halfL+1)+(j+gty.pMetaGty[1].L+u*(gty.pMetaGty[1].L+2)-1)*(gty.pMetaGty[1].L+2)] =
			Gref[i+k*gty.pMetaGty[1].halfL+(j+gty.pMetaGty[1].L*(u+1)-3)*gty.pMetaGty[1].L]
	end

	gty.G[ ( gty.pMetaGty[1].halfL+1 )*( 2gty.pMetaGty[1].L+1 ) ] = gty.gbounds[1] + rand(THREADRNG[threadid()])*(gty.gbounds[2] - gty.gbounds[1])
	gty.G[ ( gty.pMetaGty[1].halfL+1 )*( 2gty.pMetaGty[1].L+3 ) ] = gty.gbounds[1] + rand(THREADRNG[threadid()])*(gty.gbounds[2] - gty.gbounds[1])
	gty.G[ ( gty.pMetaGty[1].L+2 )*( 2gty.pMetaGty[1].L+3 ) ] = gty.gbounds[1] + rand(THREADRNG[threadid()])*(gty.gbounds[2] - gty.gbounds[1])
	gty.G[ ( gty.pMetaGty[1].L+2 )*( 2gty.pMetaGty[1].L+4 ) ] = gty.gbounds[1] + rand(THREADRNG[threadid()])*(gty.gbounds[2] - gty.gbounds[1])
	# there are a couple of connection missing in this routine
end

# function. upgrade genotypic variables
function upgradeGtyG!(gty::atSystemGty{<:atIsingMetaGty})
	Gref = copy(gty.G)
	append!( gty.G, ones(Float64, 8gty.pMetaGty[1].L+8) )

	for u in 0:1, j in 1:gty.pMetaGty[1].L, k in 0:1, i in 1:gty.pMetaGty[1].halfL
		gty.G[ i + k*(gty.pMetaGty[1].halfL + 1) + (j + u*(gty.pMetaGty[1].L + 2) - 1)*(gty.pMetaGty[1].L+2) ] =
			Gref[ i + k*gty.pMetaGty[1].halfL + (j + u*gty.pMetaGty[1].L - 1)*gty.pMetaGty[1].L ]
	end
	newg!(gty,Gref)
end

# function. upgrade metagenotype
function upgradeMetaGty!(ety::atEvotype,aMetaGty::Vector{<:tIsingSigTransMGty},gty::atSystemGty{<:tIsingSigTransMGty})
	foundMetaGty = false
	for metaGty in aMetaGty
		if gty.pMetaGty[1].L+2 == metaGty.L
			foundMetaGty = true
			gty.pMetaGty[1] = metaGty
		end
	end
	if !foundMetaGty
		push!( aMetaGty, tIsingSigTransMGty(gty.pMetaGty[1].L+Int32(2),gty.pMetaGty[1].β,gty.pMetaGty[1].he,gty.pMetaGty[1].prms) )
		gty.pMetaGty[1] = aMetaGty[end]
		# set_tDscdtEtyFactors(ety,gty)
	end
end

# function. evolutionary upgrade
function evoUpgrade!(pop::tLivingPop)
	for i in 1:pop.pN[2]
		if upgradeCondition(pop.aGty[i])
			upgradeGtyG!(pop.aGty[i])
			upgradeMetaGty!(pop.ety,pop.aMetaGty,pop.aGty[i])
			fitness!(pop.env,pop.aGty[i])
		end
	end
end

# function: set next evolutionary achievement
function setNextAch!(evo::tEvoData,gen::Integer)
	while evo.aveFitness[gen] >= 1.0 - 10^(- evo.pAveFt[1])
		evo.pAveFt[1] += evo.aveFtinc
	end
end

# function: genetic evolution a la Giardina--Kurchan--Peliti with genetic upgrade
function evolutionGKP!(pop::tLivingPop,evo::tEvoData;ubermode::Bool=false)

	@showprogress 1 "Evolutionary Dynamics Status: " for gen in 1:evo.Ngen
		evo.aveFitness[gen] = mean( [pop.aGty[i].pF[1] for i in 1:pop.pN[2]] )
		if evo.aveFitness[gen] >= 1.0 - 10^(- evo.pAveFt[1])
			push!(evo.aLivingPop,deepcopy(pop))
			push!(evo.aGen,gen)
			setNextAch!(evo,gen)
			println("\nEvolutionary achievement at generation $gen: ⟨f⟩ = ", evo.aveFitness[gen] )
		end

		evo.growthFactor[gen] = replication!(pop)
		evo.mutationFactor[gen] = effMutation!(pop)
		effSelection!(pop,ubermode)
	end
end

# function: genetic evolution a la Giardina--Kurchan--Peliti with genetic upgrade
function evolutionGKPup!(pop::tLivingPop,evo::tEvoData;ubermode::Bool=false)

	push!(evo.aLivingPop,deepcopy(pop))
	push!(evo.aGen,1)

	@showprogress 1 "Evolutionary Dynamics Status: " for gen in 1:evo.Ngen
		evo.aveFitness[gen] = mean( [pop.aGty[i].pF[1] for i in 1:pop.pN[2]] )
		if evo.aveFitness[gen] >= evo.pAveFt[1]
			push!(evo.aLivingPop,deepcopy(pop))
			push!(evo.aGen,gen)
			evo.pAveFt[1] += evo.aveFtinc
			println("\nEvolutionary achievement at generation $gen: ⟨f⟩ = ", evo.aveFitness[gen] )
		end

		evo.growthFactor[gen] = replication!(pop)
		evo.mutationFactor[gen] = effMutation!(pop)
		effSelection!(pop,ubermode)
		evoUpgrade!(pop)
	end
end

export replication!, effMutation!, effSelection!, evoUpgrade!, upgradeGtyG!, setNextAch!, evolutionGKP!, evolutionGKPup!

# =============
# |  SYSTEMS  \
# =============

function fitness!(env::atEnvironment,gty::atGenotype)
	gty.pF[1] = fitness(env,gty)
end

export fitness, fitness!

# ***********
# |  ISING  \
# ***********

# function: flipping --- for metropolis
function flipping!(istMGty::tIsingSigTransMGty, βJij::Array{<:Real,1}, βhi::Real, n::Array{<:Integer,2})
	# Monte Carlo step evaluated using Glauber transition rates:
	# next possibly transitioning spin (coordinates)
	rng = THREADRNG[threadid()]
	i, j = rand(rng,1:istMGty.L), rand(rng,1:istMGty.L)
	# transitioning?!
	if i > istMGty.li || j > istMGty.li
		if rand(rng) < ( 1. - n[i,j]*tanh(
				βJij[istMGty.Jpi[i,j]]*n[istMGty.jp[i],j] + βJij[istMGty.Jmi[i,j]]*n[istMGty.jm[i],j] +
				βJij[istMGty.Jpj[i,j]]*n[i,istMGty.jp[j]] + βJij[istMGty.Jmj[i,j]]*n[i,istMGty.jm[j]] +
				istMGty.βhe ))/2
			n[i,j] = - n[i,j]
		end
	# else
	# 	if rand(rng) < ( 1. - n[i,j]*tanh(
	# 			βJij[istMGty.Jpi[i,j]]*n[istMGty.jp[i],j] + βJij[istMGty.Jmi[i,j]]*n[istMGty.jm[i],j] +
	# 			βJij[istMGty.Jpj[i,j]]*n[i,istMGty.jp[j]] + βJij[istMGty.Jmj[i,j]]*n[i,istMGty.jm[j]] +
	# 			istMGty.βhe + βhi ))/2
	# 		n[i,j] = - n[i,j]
	# 	end
	end
end

# MonteCarlo simulation function for ising signal transduction: readout magnetization evaluation
# 	( tIsingST istMGty, interaction matrix Jij, input field h ) → readout magnetization m
function metropolis(istMGty::tIsingSigTransMGty,Jij::Array{T,1},hi::T)::Real where {T<:Real}
	βJij::Array{Float64,1} = Jij*istMGty.β;		βhi::Float64 = hi*istMGty.β

	# the initial state vector
	n = Vector{Array{Int8,2}}(undef,nthreads())
	for t in 1:nthreads()
		n[t] = Int8(istMGty.he != 0.0 ? sign(istMGty.he) : 1).*ones(Int8, istMGty.L, istMGty.L)
		n[t][1:istMGty.li,1:istMGty.li] .= sign(hi)
	end

	# definition: time-averaged readout magnetization
	mro = zeros(Float64, nthreads())

	NsmplsThreaded = istMGty.prms.Nsmpl ÷ nthreads()
	# istMGty.prms.Nsmpl % nthreads() == 0 || throw("mismatch between nthreads and Nsmpl. Choose: Nsmpl = multiple of $nthreads()")
	for t in 1:nthreads()
		for is in 1:NsmplsThreaded
			for imcs in 1:istMGty.prms.Nmcsps, ilp in 1:istMGty.L2
				flipping!(istMGty,βJij,βhi,n[threadid()])
			end
			# evaluation: time-averaged readout magnetization
			mro[threadid()] += sum(n[threadid()][istMGty.halfL+1:istMGty.halfL+istMGty.li,istMGty.halfL+1:istMGty.halfL+istMGty.li])
		end
	end
	# evaluation: time-averaged readout magnetization
	return sum(mro)/(NsmplsThreaded*nthreads()*istMGty.li2)
end

# MonteCarlo simulation function for ising signal transduction: time-averaged spin config evaluation
# 	( tIsingST istMGty, state n, interaction matrix Jij, input field h ) → readout magnetization m
function metropolis!(istMGty::tIsingSigTransMGty,Jij::Array{<:Real,1},hi::Real,an::Vector{Array{Int8,2}},prms::tDTMCprm)
	βJij::Array{Float64,1} = Jij.*istMGty.β;		βhi::Float64 = hi.*istMGty.β

	n = Vector{Array{Int8,2}}(undef,nthreads())
	@threads for t in 1:nthreads()
		n[t] = Int8(istMGty.he != 0.0 ? sign(istMGty.he) : 1).*ones(Int8, istMGty.L, istMGty.L)
		n[t][1:istMGty.li,1:istMGty.li] .= sign(hi)
	end

	# @showprogress 1 "Monte Carlo Dynamics Status: "
	NsmplsThreaded = prms.Nsmpl ÷ nthreads()
	prms.Nsmpl % nthreads() == 0 || throw("mismatch between nthreads and Nsmpl. Choose: Nsmpl = multiple of $nthreads()")
	@threads for t in 1:nthreads()
		for is in 1:NsmplsThreaded
			for imcs in 1:prms.Nmcsps, ilp in 1:istMGty.L2
				flipping!(istMGty,βJij,βhi,n[threadid()])
			end
			# an[ is + (threadid() - 1)*NsmplsThreaded ] = deepcopy( n[threadid()] )
			an[ is + (t - 1)*NsmplsThreaded ] = deepcopy( n[threadid()] )	# <- when no threading is used
		end
	end
end

# function: fitness for ising signal transduction
function fitness(env::tCompEnv{<:Vector{<:Vector{<:Real}}},gty::atSystemGty{<:atIsingMetaGty})::Float64
	fValues = zeros(Float64,gty.pMetaGty[1].prms.Ntrials)
	for t in 1:gty.pMetaGty[1].prms.Ntrials
		d2::Float64 = 0.0
		for io in env.IOidl
			d2 += ( metropolis(gty.pMetaGty[1],broadcast(x->10^(x),gty.G),io[1]) - io[2] )^2
		end
		fValues[t] = exp(-sqrt(d2)/env.selFactor)
	end
	return minimum(fValues) + (gty.pMetaGty[1].halfL - BASALFITNESS)
end

export metropolis, metropolis!

# **************
# |  CHANNELS  \
# **************

# constructor: disordered channel
function tDisChnMGty(L::Int32, p::Float64, kout::Float64)
	L2mL = L*(L-1)

	viability = rand(THREADRNG[threadid()],Bernoulli(p),2L2mL);	Nvbl = sum(viability)
	V = [ Vector{Int32}(undef, 2Nvbl), Vector{Int32}(undef, 2Nvbl) ]

	g = 0
	for (i,v) in enumerate(viability)
		if i <= L2mL && Bool(v)
			g+=1
			V[1][g] = i + (i-1)÷(L-1) + 1;	V[2][g] = i + (i-1)÷(L-1);
			V[1][g+Nvbl] = V[2][g];			V[2][g+Nvbl] = V[1][g];
		elseif i > L2mL && Bool(v)
			g+=1
			V[1][g] = i - L2mL + L;		V[2][g] = i - L2mL;
			V[1][g+Nvbl] = V[2][g];		V[2][g+Nvbl] = V[1][g];
		end
	end

	mEvoTypes.tDisChnMGty(L, L^2, L2mL, 2Nvbl, Nvbl, p, V, kout)
end

function getW(gty::atSystemGty{<:atChannelMetaGty})
	W = sparse(gty.pMetaGty[1].V..., broadcast(x->10^(x),gty.G), gty.pMetaGty[1].L2+1, gty.pMetaGty[1].L2+1)
	W[end,gty.pMetaGty[1].L2mL+1:gty.pMetaGty[1].L2] .= gty.pMetaGty[1].kout
	for i in 1:gty.pMetaGty[1].L2
		W[i,i] = -sum(W[:,i])
	end
	return W
end

# function: response of channel using forward differentiation
function responseFD(gty::atSystemGty{<:atChannelMetaGty},W::AbstractMatrix,Input::Vector{<:Real} )
	length(Input) == gty.pMetaGty[1].L || throw( "inconsistent dimensions between input currents and system" )
	W[1:gty.pMetaGty[1].L,end] = Input
	W[end,end] = -sum(W[1:end-1,end])

	Λ(q::Vector) = det( Array(W .* (1. .+ sparse(fill(gty.pMetaGty[1].L2+1,gty.pMetaGty[1].L), gty.pMetaGty[1].L2mL+1:gty.pMetaGty[1].L2,
		map(e -> - 1. + exp(e), q[1:end-1]),gty.pMetaGty[1].L2+1,gty.pMetaGty[1].L2+1)) ) - q[end] .* I )
	dΛ(q::Vector) = ForwardDiff.gradient(Λ,q)

	q0 = zeros(Float64,gty.pMetaGty[1].L+1)
	y = dΛ(q0)

	return [ -y[i]/y[end] for i in 1:gty.pMetaGty[1].L ]
end

# function: constructor of normalized stochastic matrix for disordered channel metagenotype
# The last column is kept null as it is reserved for input transitions. The last row encodes the normalization.
function getWnrmd(gty::atSystemGty{<:atChannelMetaGty})
	W = sparse(gty.pMetaGty[1].V..., broadcast(x->10^(x),gty.G), gty.pMetaGty[1].L2+1, gty.pMetaGty[1].L2+1)
	W[end,gty.pMetaGty[1].L2mL+1:gty.pMetaGty[1].L2] .= gty.pMetaGty[1].kout
	for i in 1:gty.pMetaGty[1].L2
		W[i,i] = -sum(W[:,i])
	end
	W[end,:] = ones(Float64,gty.pMetaGty[1].L2+1)
	return W
end

# function: response of channel using linear solver
function response(gty::atSystemGty{<:atChannelMetaGty},W::AbstractMatrix,Input::Vector{<:Real})
	length(Input) == gty.pMetaGty[1].L || throw( "inconsistent dimensions between input currents and system" )
	W[1:gty.pMetaGty[1].L,end] = Input
	b = zeros(Float64,gty.pMetaGty[1].L2+1);	b[end] = 1.0

	return gty.pMetaGty[1].kout .* ( ( W \ b )[gty.pMetaGty[1].L2mL+1:gty.pMetaGty[1].L2] )
end

# function: response of channel using linear solver -- compact version
function response(gty::atSystemGty{<:atChannelMetaGty},Input::Vector{<:Real})
	W = getWnrmd(gty)
	response(gty, W, Input)
end

function fitness(istEnv::tCompEnv{<:Vector{<:Vector{<:Vector{<:Real}}}},gty::atSystemGty{<:atChannelMetaGty})::Float64
	W = getWnrmd(gty)

	d2::Float64 = 0.0
	for io in istEnv.IOidl
		d2 += sum( (response(gty, W, io[1]) - io[2]).^2 )
	end
	return exp(-d2/istEnv.selFactor)
end

export tDisChnMGty, response


# *************************
# | STATISTICS FUNCTIONS  \
# *************************

function getGStat!(pop::tLivingPop,GAve::Vector{Float64},GCov::Array{Float64,2},GCor::Array{Float64,2})
	mapreduce( x -> x ==  pop.aGty[1].pMetaGty[1].dG, &, [ pop.aGty[i].pMetaGty[1].dG for i in 2:pop.pN[2] ] ) ||
		throw(DimensionMismatch("inconsistent dimensions"))

	for x in 1:pop.aGty[1].pMetaGty[1].dG
		GAve[x] = mean( [pop.aGty[i].G[x] for i in 1:pop.pN[2]] )
	end

	for y in 1:pop.aGty[1].pMetaGty[1].dG, x in 1:pop.aGty[1].pMetaGty[1].dG
		GCov[x,y] = myCov( [ pop.aGty[i].G[x] for i in 1:pop.pN[2] ], GAve[x], [ pop.aGty[i].G[y] for i in 1:pop.pN[2] ], GAve[y] )
	end

	for y in 1:pop.aGty[1].pMetaGty[1].dG, x in 1:pop.aGty[1].pMetaGty[1].dG
		GCor[x,y] = GCov[x,y]/sqrt(GCov[x,x]*GCov[y,y])
	end
end

export getGStat!

# **********************************
# | STATISTICS FUNCTIONS --- ISING \
# **********************************

function myCov(X::AbstractArray,XAve::Number,Y::AbstractArray,YAve::Number,)
	N = length(X)
	N == length(Y) || throw(DimensionMismatch("inconsistent dimensions"))
	cov = 0.0
	for (x, y) in zip(X,Y)
		cov += (x - XAve)*(y - YAve)
	end
	return cov/(N-1)
end

# function: showing the spin config of the ising signal transduction system
function getSpinStat!(env::tCompEnv,gty::atSystemGty{<:atIsingMetaGty},aSpinAve::Vector{Array{Float64,2}},
		aSpinCov::Vector{Array{Float64,2}},sCor::Array{Float64,2},prms::tDTMCprm)
	aan = Vector{Vector{Array{Int8,2}}}(undef,length(env.IOidl))
	aSpinCorFisherz = Vector{Array{Float64,2}}(undef,length(env.IOidl))

	for i in eachindex(env.IOidl)
		aan[i] = [ Array{Int8}(undef, gty.pMetaGty[1].L, gty.pMetaGty[1].L) for ismpl in 1:prms.Nsmpl ]
		aSpinCorFisherz[i] = Array{Float64}(undef, gty.pMetaGty[1].L2, gty.pMetaGty[1].L2)

		metropolis!( gty.pMetaGty[1], broadcast(x->10^(x),gty.G), env.IOidl[i][1], aan[i], prms )

		for s in 1:gty.pMetaGty[1].L2
			aSpinAve[i][s] = mean([aan[i][t][s] for t in 1:prms.Nsmpl])
		end

		for sj in 1:gty.pMetaGty[1].L2, si in 1:gty.pMetaGty[1].L2
			aSpinCov[i][si,sj] = myCov( [aan[i][t][si] for t in 1:prms.Nsmpl],aSpinAve[i][si],[aan[i][t][sj] for t in 1:prms.Nsmpl],aSpinAve[i][sj] )
		end

		for sj in 1:gty.pMetaGty[1].L2, si in 1:gty.pMetaGty[1].L2
			# aSpinCorFisherz[i][si,sj] = map( r -> log( (1 + r)/(1 - r) )/2, aSpinCov[i][si,sj]/sqrt(aSpinCov[i][sj,sj]*aSpinCov[i][si,si]) )
			aSpinCorFisherz[i][si,sj] = aSpinCov[i][si,sj]/sqrt(aSpinCov[i][sj,sj]*aSpinCov[i][si,si])
		end
	end

	for sj in 1:gty.pMetaGty[1].L2, si in 1:gty.pMetaGty[1].L2
		# sCor[si,sj] = tanh( mean( [ aSpinCorFisherz[i][si,sj] for i in 1:length(env.IOidl)] ) )
		sCor[si,sj] = mean( [ aSpinCorFisherz[i][si,sj] for i in eachindex(env.IOidl)] )
	end
end

function showG!(gty::tCntGty{<:atIsingMetaGty},GMat::Array{Float64,2})
	ii::Int32 = 0
	for j in 1:2gty.pMetaGty[1].L, i in 1:2gty.pMetaGty[1].L
		#	column is odd => ( row is even => G ; otherwise spin place ) ; otherwise ( row is odd => G ; otherwise empty space )
		GMat[i,j] = j%2==1 ? ( i%2==0 ? gty.G[ii+=1] : gty.gbounds[1] - 10 ) : ( i%2==1 ? gty.G[ii+=1] : gty.gbounds[2] + 10 )
	end
end

function showJij!(gty::atSystemGty{<:atIsingMetaGty},JijMat::Array{Float64,2})
	ii::Int32 = 0
	for j in 1:2gty.pMetaGty[1].L, i in 1:2gty.pMetaGty[1].L
		JijMat[i,j] = j%2==1 ? ( i%2==0 ? 10^(gty.G[ii+=1]) : -1 ) : ( i%2==1 ? 10^(gty.G[ii+=1]) : 10^6+1 )
	end
end

function getJijStat!(pop::tLivingPop,JijAve::Vector{Float64},JijCov::Array{Float64,2},JijCor::Array{Float64,2})
	mapreduce( x -> x ==  pop.aGty[1].pMetaGty[1].dG, &, [ pop.aGty[i].pMetaGty[1].dG for i in 2:pop.pN[2] ] ) ||
		throw(DimensionMismatch("inconsistent dimensions"))

	for x in 1:pop.aGty[1].pMetaGty[1].dG
		JijAve[x] = mean( [ 10^pop.aGty[i].G[x] for i in 1:pop.pN[2]] )
	end

	for y in 1:pop.aGty[1].pMetaGty[1].dG, x in 1:pop.aGty[1].pMetaGty[1].dG
		JijCov[x,y] = myCov( [ 10^pop.aGty[i].G[x] for i in 1:pop.pN[2] ], JijAve[x], [ 10^pop.aGty[i].G[y] for i in 1:pop.pN[2] ], JijAve[y] )
	end

	for y in 1:pop.aGty[1].pMetaGty[1].dG, x in 1:pop.aGty[1].pMetaGty[1].dG
		JijCor[x,y] = JijCov[x,y]/sqrt(JijCov[x,x]*JijCov[y,y])
	end
end

function getJij!(gty::atSystemGty{<:atIsingMetaGty},JijMat::Array{Float64,2})
	JijMat .= 0
	for x in 1:gty.pMetaGty[1].L, y in 1:gty.pMetaGty[1].L
		JijMat[ x+(y-1)*gty.pMetaGty[1].L, gty.pMetaGty[1].jp[x]+(y-1)*gty.pMetaGty[1].L ] = 10.0^gty.G[ gty.pMetaGty[1].Jpi[x,y] ]
		JijMat[ x+(y-1)*gty.pMetaGty[1].L, gty.pMetaGty[1].jm[x]+(y-1)*gty.pMetaGty[1].L ] = 10.0^gty.G[ gty.pMetaGty[1].Jmi[x,y] ]
		JijMat[ x+(y-1)*gty.pMetaGty[1].L, x+(gty.pMetaGty[1].jp[y]-1)*gty.pMetaGty[1].L ] = 10.0^gty.G[ gty.pMetaGty[1].Jpj[x,y] ]
		JijMat[ x+(y-1)*gty.pMetaGty[1].L, x+(gty.pMetaGty[1].jm[y]-1)*gty.pMetaGty[1].L ] = 10.0^gty.G[ gty.pMetaGty[1].Jmj[x,y] ]
	end
end

export getSpinStat!, showJij!, getJijStat!

# ************************************
# | STATISTICS FUNCTIONS --- CHANNEL \
# ************************************

# function getW(gty::atSystemGty{<:atChannelMetaGty})
# 	W = sparse(gty.pMetaGty[1].V..., broadcast(x->10^(x),gty.G), gty.pMetaGty[1].L2+1, gty.pMetaGty[1].L2+1)
# 	W[end,gty.pMetaGty[1].L2mL+1:gty.pMetaGty[1].L2] .= gty.pMetaGty[1].kout
# 	for i in 1:gty.pMetaGty[1].L2
# 		W[i,i] = -sum(W[:,i])
# 	end
# 	return W
# end
#
# # function: response of channel using forward differentiation
# function responseFD(gty::atSystemGty{<:atChannelMetaGty},W::AbstractMatrix,Input::Vector{<:Real} )
# 	length(Input) == gty.pMetaGty[1].L || throw( "inconsistent dimensions between input currents and system" )
# 	W[1:gty.pMetaGty[1].L,end] = Input
# 	W[end,end] = -sum(W[1:end-1,end])
#
# 	Λ(q::Vector) = det( Array(W .* (1. .+ sparse( fill(gty.pMetaGty[1].L2+1,gty.pMetaGty[1].L), gty.pMetaGty[1].L2mL+1:gty.pMetaGty[1].L2,
# 		map(e -> exp(e) - 1.0, q[1:end-1]),gty.pMetaGty[1].L2+1,gty.pMetaGty[1].L2+1)) ) - q[end] .* I )
# 	dΛ(q::Vector) = ForwardDiff.gradient(Λ,q)
#
# 	q0 = zeros(Float64,gty.pMetaGty[1].L+1)
# 	y = dΛ(q0)
#
# 	return [ -y[i]/y[end] for i in 1:gty.pMetaGty[1].L ]
# end
#
# # function: showing the spin config of the ising signal transduction system
# function getCrntStat!(env::tCompEnv,gty::atSystemGty{<:atChannelMetaGty},crntAve::tCrntPattern,
# 		aCrntCov::Vector{Vector{Float64}},crntCor::Vector{Float64})
#
# 	# aCrntAve = Vector{Vector{Array{Int8,2}}}(undef,length(env.IOidl))
# 	# aCrntCorFisherz = Vector{Array{Float64,2}}(undef,length(env.IOidl))
#
# 	W = getW(gty)
#
# 	Vout = [ fill(gty.pMetaGty[1].L2+1,gty.pMetaGty[1].L), collect(gty.pMetaGty[1].L2mL+1:gty.pMetaGty[1].L2) ]
# 	Vin = [ collect(1:gty.pMetaGty[1].L), fill(gty.pMetaGty[1].L2+1,gty.pMetaGty[1].L) ]
#
# 	crntAve.V[1] = copy( vcat( gty.pMetaGty[1].V[1][gty.pMetaGty[1].Nvbl+1:gty.pMetaGty[1].dG],
# 		fill(gty.pMetaGty[1].L2+1,gty.pMetaGty[1].L), collect(1:gty.pMetaGty[1].L) ) )
# 	crntAve.V[2] = copy( vcat( gty.pMetaGty[1].V[2][gty.pMetaGty[1].Nvbl+1:gty.pMetaGty[1].dG],
# 		collect(gty.pMetaGty[1].L2mL+1:gty.pMetaGty[1].L2), fill(gty.pMetaGty[1].L2+1,gty.pMetaGty[1].L) ) )
#
# 	q0 = zeros(Float64, gty.pMetaGty[1].Nvbl + 2gty.pMetaGty[1].L + 1 )
# 	for (i, io) in enumerate(env.IOidl)
# 		W[1:gty.pMetaGty[1].L,end] = io[1]
# 		W[end,end] = -sum(W[1:end-1,end])
#
# 		Λ(q::Vector) = det( Array(W .*
# 			( 1.0 .+ sparse( crntAve.V[1][1:gty.pMetaGty[1].Nvbl], crntAve.V[2][1:gty.pMetaGty[1].Nvbl],
# 				map(e -> exp(e) - 1.0, q[1:gty.pMetaGty[1].Nvbl]), gty.pMetaGty[1].L2+1, gty.pMetaGty[1].L2+1 ) ) .*
# 			( 1.0 .+ sparse( crntAve.V[1][gty.pMetaGty[1].Nvbl+1:gty.pMetaGty[1].dG], crntAve.V[2][gty.pMetaGty[1].Nvbl+1:gty.pMetaGty[1].dG],
# 				map(e -> exp(-e) - 1.0, q[1:gty.pMetaGty[1].Nvbl]), gty.pMetaGty[1].L2+1, gty.pMetaGty[1].L2+1) ) .*
# 			( 1.0 .+ sparse( Vout..., map(e -> exp(e) - 1.0, q[gty.pMetaGty[1].Nvbl+1:aDisChnGty[1].pMetaGty[1].Nvbl+gty.pMetaGty[1].L]),
# 					gty.pMetaGty[1].L2+1, gty.pMetaGty[1].L2+1) ) .*
# 			( 1.0 .+ sparse( Vin..., map(e -> exp(e) - 1.0, q[gty.pMetaGty[1].Nvbl+gty.pMetaGty[1].L+1:aDisChnGty[1].pMetaGty[1].Nvbl+2gty.pMetaGty[1].L]),
# 					gty.pMetaGty[1].L2+1, gty.pMetaGty[1].L2+1) )
# 			) - q[end] .* I )
# 		dΛ(q::Vector) = ForwardDiff.gradient(Λ,q)
#
# 		y = dΛ(q0)
# 		[ -y[i]/y[end] for i in 1:gty.pMetaGty[1].L ]
# 	end
#
# 	# for i in eachindex(env.IOidl)
# 	# 	aan[i] = [ Array{Int8}(undef, gty.pMetaGty[1].L, gty.pMetaGty[1].L) for ismpl in 1:prms.Nsmpl ]
# 	# 	aSpinCorFisherz[i] = Array{Float64}(undef, gty.pMetaGty[1].L2, gty.pMetaGty[1].L2)
# 	#
# 	# 	metropolis!( gty.pMetaGty[1], broadcast(x->10^(x),gty.G), env.IOidl[i][1], aan[i], prms )
# 	#
# 	# 	for s in 1:gty.pMetaGty[1].L2
# 	# 		aSpinAve[i][s] = mean([aan[i][t][s] for t in 1:prms.Nsmpl])
# 	# 	end
# 	#
# 	# 	for sj in 1:gty.pMetaGty[1].L2, si in 1:gty.pMetaGty[1].L2
# 	# 		aSpinCov[i][si,sj] = myCov( [aan[i][t][si] for t in 1:prms.Nsmpl],aSpinAve[i][si],[aan[i][t][sj] for t in 1:prms.Nsmpl],aSpinAve[i][sj] )
# 	# 	end
# 	#
# 	# 	for sj in 1:gty.pMetaGty[1].L2, si in 1:gty.pMetaGty[1].L2
# 	# 		# aSpinCorFisherz[i][si,sj] = map( r -> log( (1 + r)/(1 - r) )/2, aSpinCov[i][si,sj]/sqrt(aSpinCov[i][sj,sj]*aSpinCov[i][si,si]) )
# 	# 		aSpinCorFisherz[i][si,sj] = aSpinCov[i][si,sj]/sqrt(aSpinCov[i][sj,sj]*aSpinCov[i][si,si])
# 	# 	end
# 	# end
# 	#
# 	# for sj in 1:gty.pMetaGty[1].L2, si in 1:gty.pMetaGty[1].L2
# 	# 	# sCor[si,sj] = tanh( mean( [ aSpinCorFisherz[i][si,sj] for i in 1:length(env.IOidl)] ) )
# 	# 	sCor[si,sj] = mean( [ aSpinCorFisherz[i][si,sj] for i in eachindex(env.IOidl)] )
# 	# end
# end

# *******************
# I/O FUNCTIONS
# *******************

# function: saving the Discrete Time Monte Carlo parameters
function write_DTMCprm(prms::tDTMCprm,fileTag::String)
	open( fileTag * "_DTMCprm" * ".dat", "w" ) do f
		print(f, prms.Nsmpl, "\t")
		print(f, prms.Nmcsps, "\t")
		print(f, prms.Ntrials, "\t")
	end
end

# function: saving the MetaGenotype of Ising systmes
function write_MetaGty(metaGty::atIsingMetaGty,fileTag::String)
	open( fileTag * "_MetaGty" * ".dat", "w" ) do f
		print(f, metaGty.β, "\t")
		print(f, metaGty.he, "\t")
	end
end

# function: saving the genotypes
function write_aGty(aGty::Vector{<:atGenotype},Npop::Int32,fileTag::String)
	length(aGty) >= Npop || throw(DimensionMismatch("population size exceeds Npop variable"))
	open( fileTag * "_a" * split(string(typeof(aGty[1])),"{")[1] * ".dat", "w" ) do f
		for i in 1:Npop
			print(f, aGty[i].pMetaGty[1].dG, "\t")
			for x in aGty[i].G
				print(f, x, "\t")
			end
			print(f, aGty[i].pF[1], "\n")
		end
	end
end

function write_aGty(aGty::Vector{tAlphaGty},Npop::Int32,fileTag::String)
	length(aGty) >= Npop || throw(DimensionMismatch("population size exceeds Npop variable"))
	open( fileTag * "_a" * split(string(typeof(aGty[1])),"{")[1] * ".dat", "w" ) do f
		for i in 1:Npop
			print(f, aGty[i].pMetaGty[1].dG, "\t")
			for x in aGty[i].G
				print(f, x, "\t")
			end
			print(f, aGty[i].pF[1], "\n")
		end
	end
	open( fileTag * "_a" * split(string(typeof(aGty[1])),"{")[1] * ".dat", "w" ) do f
		print(f, aGty[i].pdg[1], "\t")
		for x in aGty[1].g
			print(f, x, "\t")
		end
	end
end

# function: saving the ising population
function write_tLivingPop(pop::tLivingPop{<:atEvotype,<:atEnvironment,<:Vector{<:atIsingMetaGty},<:Vector{<:atGenotype}},fileTag::String)
	write_DTMCprm(pop.aGty[1].pMetaGty[1].prms,fileTag)
	write_MetaGty(pop.aGty[1].pMetaGty[1],fileTag)
	write_aGty(pop.aGty,pop.pN[2],fileTag)
end

# function: saving the evolutionary data
function write_tEvoData(aData::Vector{tEvoData},fileTag::String)
	open( "data/" * fileTag * "_evoData" * ".dat", "w" ) do f
		for dataBatch in aData, t in 1:dataBatch.Ngen print(f,dataBatch.aveFitness[t],"\t") end, print(f,"\n")
		for dataBatch in aData, t in 1:dataBatch.Ngen print(f,dataBatch.growthFactor[t],"\t") end, print(f,"\n")
		for dataBatch in aData, t in 1:dataBatch.Ngen print(f,dataBatch.mutationFactor[t],"\t") end, print(f,"\n")
	end

	for dataBatch in aData
		for (gen, pop)  in zip(dataBatch.aGen, dataBatch.aLivingPop)
			write_tLivingPop( pop, "data/" * fileTag * "_g$gen")
		end
	end
end

# to do: saving evolution parameters

export write_DTMCprm, write_MetaGty, write_aGty, write_tLivingPop, write_tEvoData

# read metagenotype
function read_MetaGty(L::Integer,prms::atMonteCarloPrm,fileTag::String)::atIsingMetaGty
	metaGtyMat = readdlm( "data/" * fileTag * "_MetaGty" * ".dat" )
	return tIsingSigTransMGty(L,metaGtyMat[1],metaGtyMat[2],prms)
end

# read population and spits metagenotype and genotype
function read_aIsingSigTransAlphaGty(prms::atMonteCarloPrm,fileTag::String)
	gtyMat = readdlm( "data/" * fileTag * "_atAlphaGty" * ".dat" )
	gMat = readdlm( "data/" * fileTag * "_atAlphaGty_g" * ".dat" )
	Npop = size(gtyMat)[1]

	sysSizes = Int32[]
	L = zeros(Int32,Npop)

	for i in 1:Npop
		L[i] = Int32( sqrt(gtyMat[i,1]/2) )
		if !( L[i] in sysSizes )
			push!(sysSizes,L[i])
		end
	end

	aMetaGty = [ read_MetaGty(L,prms,fileTag) for L in sysSizes ]
	aGty = [ tAlphaGty( [aMetaGty[collect(1:length(sysSizes))[sysSizes .== L[i]][1]]], gtyMat[i,1],
		gtyMat[i,2:1+Int32(gtyMat[i,1])], Int32[gMat[1]], gMat[2:end], Float64[gtyMat[i,end]] ) for i in 1:Npop ]
	return aMetaGty, aGty, Npop
end

export read_aIsingSigTransGty

# *******************
# TRIVIAL
# *******************

function fitness!(trivialEty::tTrivialEty,trivialEnv::tTrivialEnv,gty::atGenotype)
	gty.pF[1]=1/(euclidean(gty.G,ones(Float64,gty.pMetaGty[1].dG))+FITNESSOFFSET)
end

export fitness!

end
