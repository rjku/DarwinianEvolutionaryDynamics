
module mEvoFunc

using Base.Threads, SparseArrays
using LinearAlgebra, Random, Statistics, Distributions, ForwardDiff, Distances
using mEvoTypes, mUtils
using DelimitedFiles, Dates, ProgressMeter
import Future

# threads random generators initialization
const THREADRNG = let m = MersenneTwister(1)
            [m; accumulate(Future.randjump, fill(big(10)^20, nthreads()-1), init=m)]
        end;

const FITNESSOFFSET, FITNESSTHRESHOLD, BASALFITNESS, MAXSIZE = 1.0, .99, 3.0, 8

function initEvoPop( N::Int32,ety::Tety,env::Tenv,aMetaGty::Vector{Tmgty},aGty::Vector{TG} ) where {
		Tety<:atEvotype,Tenv<:atEnvironment,Tmgty<:atMetaGenotype,TG<:atGenotype }
	N <= length(aGty) || throw(DimensionMismatch("Inconsistent Dimensions: N > length(aGty)"))
	@threads for i in 1:N fitness!(aGty[i],env) end
	return tEvoPop{Tety,Tenv,Vector{Tmgty},Vector{TG}}( Int32[N,N,length(aGty)],ety,env,aMetaGty,aGty )
end

# function. changing rep and mut -factors in tDscdtEty
function set_tDscdtEtyFactors(ety::tDscdtEty,gty::atGenotype)
	ety.pRepFactor[1] = ety.repRate/(2gty.pMetaGty[1].dG*ety.mutRate+ety.ΔtOffset)
	ety.pMutFactor[1] = ety.mutRate/(2gty.pMetaGty[1].dG*ety.mutRate+ety.ΔtOffset)
end

# function. changing rep and mut -factors in tDscdtEty
function set_tDscdtEtyFactors(ety::tDscdtEty,gty::tAlphaGty)
	ety.pRepFactor[1] = ety.repRate/(length(gty.g)*gty.pMetaGty[1].dG*ety.mutRate+ety.ΔtOffset)
	ety.pMutFactor[1] = ety.mutRate/(length(gty.g)*gty.pMetaGty[1].dG*ety.mutRate+ety.ΔtOffset)
end

export initEvoPop


# *********************************
# | PROPER EVOLUTIONARY FUNCTIONS \
# *********************************

# function. population replication optimized for populations of individuals (one niches)
function replicationOne!(pop::tEvoPop)
	Kr = pop.ety.pRepFactor[1]*pop.aGty[1].aF[2]							# effective replication factor
	ipKr = trunc(Int32,Kr)													# floored effective replication factor
	G::Int32 = rand(THREADRNG[threadid()]) < Kr - ipKr ? ipKr + 1 : ipKr	# fluctuating growth coefficient

	for inew in 1:G															# offspring production
		pop.pN[1] += 1
		if pop.pN[1] <= length(pop.aGty)
			pop.aGty[pop.pN[1]] = copy(pop.aGty[1])
		else
			push!(pop.aGty,copy(pop.aGty[1]))
		end
	end

	return log(pop.pN[1])
end

# function. generic population replication
function replication!(pop::tEvoPop)
	Kr = Vector{Float64}(undef,nthreads())
	ipKr = Vector{Int32}(undef,nthreads())
	G = zeros(Int32,pop.pN[2])

	@threads for i in 1:pop.pN[2]
		Kr[threadid()] = pop.ety.pRepFactor[1]*pop.aGty[i].aF[2]
		ipKr[threadid()] = trunc(Int32,Kr[threadid()])
		G[i] = rand(THREADRNG[threadid()]) < Kr[threadid()] - ipKr[threadid()] ? ipKr[threadid()] + 1 : ipKr[threadid()]
	end

	for i in 1:pop.pN[2]
		for inew in 1:G[i]
			pop.pN[1] += 1
			if pop.pN[1] <= length(pop.aGty)
				pop.aGty[pop.pN[1]] = copy(pop.aGty[i])
			else
				push!(pop.aGty,copy(pop.aGty[i]))
			end
		end
	end

	return log(pop.pN[1]/pop.pN[2])
end

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

function mutation!(gty::atGenotype,ety::atGtyMutEty)
	Pmut = ety.pMutFactor[1]/(1+ety.pRepFactor[1]*gty[i].aF[2])
	CPmut = 2length(gty)*Pmut
	CPmut <= 1 || throw("cumulative probability of mutation exceeds 1")

	CDFmut::Float64 = Pmut; i::Int32 = 0; r::Float64 = rand(THREADRNG[threadid()])
	if r <= CPmut
		while CDFmut < r
			CDFmut += Pmut
			i += 1
		end
		mutation!(gty,i)

		return Int32(1)
	else
		return Int32(0)
	end
end

function mutation!(gty::atGenotype,ety::atPntMutEty)
	K = 1.0/(1.0 + ety.pRepFactor[1]*gty.aF[2])
	CDFmut = 1.0-K*(1.0-(1.0-ety.pPntMutFactor[1])^length(gty))		# probability of no mutation
	Nmut::Int32 = 0;	rNmut = rand(THREADRNG[threadid()])

	# determining the number of mutations
	while CDFmut < rNmut && Nmut < length(ety.aSize2PMF[length(gty)])
		CDFmut += K*ety.aSize2PMF[length(gty)][Nmut+=1]
	end

	# determining the mutating genes
	aMut = Vector{Int32}(undef, Nmut)
	for nmut in 1:Nmut
		aMut[nmut] = rand(THREADRNG[threadid()], filter( e -> !(e in aMut), 1:length(gty) ))
	end

	# determining the mutations
	for i in aMut
		mutation!(gty,i)
	end

	return Nmut
end

# function. effective mutation: blindly mutate the population's genotype according to the effective dynamical mutation rate
function mutationOne!(pop::tEvoPop)
	aNmut = Vector{Int32}(undef,pop.pN[1])
	for i in 1:pop.pN[1]
		aNmut[i] = mutation!(pop.aGty[i],pop.ety)
		if aNmut[i] > 0
			fitness!(pop.aGty[i],pop.env)
		end
	end
	return sum(aNmut)
end

# function. effective mutation: blindly mutate the population's genotype according to the effective dynamical mutation rate
function mutation!(pop::tEvoPop)
	aNmut = Vector{Int32}(undef,pop.pN[1])
	@threads for i in 1:pop.pN[1]
		aNmut[i] = mutation!(pop.aGty[i],pop.ety)
		if aNmut[i] > 0
			fitness!(pop.aGty[i],pop.env)
		end
	end
	return sum(aNmut)/pop.pN[2] 	# normalized number of mutations
end

function selectionOne(aSelFncVals::Vector{<:Real})
	i, CF = 1, aSelFncVals[1]
	r = sum(aSelFncVals)*rand(THREADRNG[threadid()])

	while CF < r
		CF += aSelFncVals[i += 1]
	end

	return i
end

# function. selection of the fittest within the niche population
function selectionOne!(pop::tEvoPop,elite::Bool=false)
	if elite
		fmax = maximum([ gty.aF[3] for gty in pop.aGty[1:pop.pN[1]] ])
		iSelected = findfirst( f -> f == fmax, [ gty.aF[2] for gty in pop.aGty[1:pop.pN[1]] ] )
	else
		iSelected = selectionOne([ pop.aGty[i].aF[3] for i in 1:pop.pN[1] ])
	end

	pop.aGty[1], pop.pN[1] = pop.aGty[iSelected], pop.pN[2]
end

function selection!(aSelFncVals::Vector{<:Real},aiSelected::Vector{<:Integer})
	Ftot = sum(aSelFncVals)

	@threads for j in eachindex(aiSelected)
		i, CF = 1, aSelFncVals[1]
		r = Ftot*rand(THREADRNG[threadid()])

		while CF < r
			CF += aSelFncVals[i += 1]
		end

		aiSelected[j] = i
	end
end

# function. effective selection: pruning of the population
function selection!(pop::tEvoPop,elite::Bool=false)
	popGtyRef::Array{atGenotype,1} = copy(pop.aGty)

	if elite
		# survival of the fittest
		aiSelected = sortperm([ pop.aGty[i].aF[3] for i in 1:pop.pN[1] ],rev=true)
		# sort!(pop.aGty, by= x -> x.aF[2], rev=true)
	else
		# selection with replacement of individuals based on selection function values
		aiSelected = zeros(Int32,pop.pN[2])
		selection!([ pop.aGty[i].aF[3] for i in 1:pop.pN[1] ], aiSelected)
		# aiSelected = rand(THREADRNG[threadid()],1:pop.pN[1],pop.pN[2])
	end
	for i in 1:pop.pN[2]
		pop.aGty[i] = popGtyRef[aiSelected[i]]
	end

	# population size renormalization
	pop.pN[1] = pop.pN[2]
end

# function. condition for genotypic upgrade
function upgradeCondition(gty::atSystemGty{<:atIsingMetaGty},maxSize::Integer)
	return gty.aF[2] - (gty.pMetaGty[1].halfL - BASALFITNESS) > FITNESSTHRESHOLD && gty.pMetaGty[1].L < maxSize
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
function upgradeMetaGty!(ety::atEvotype,aMetaGty::Vector{<:tIsingSigTransMetaGty},gty::atSystemGty{<:tIsingSigTransMetaGty},evo::tEvoData)
	foundMetaGty = false

	im = 2
	while !foundMetaGty && im <= length(aMetaGty)
		if gty.pMetaGty[1].L+2 == aMetaGty[im].L
			foundMetaGty = true
			gty.pMetaGty[1] = aMetaGty[im]
		end
	end

	# for metaGty in aMetaGty
	# 	if gty.pMetaGty[1].L+2 == metaGty.L
	# 		foundMetaGty = true
	# 		gty.pMetaGty[1] = metaGty
	# 		break
	# 	end
	# end

	if !foundMetaGty
		push!( aMetaGty, tIsingSigTransMetaGty(gty.pMetaGty[1].L+Int32(2),gty.pMetaGty[1].β,gty.pMetaGty[1].he,gty.pMetaGty[1].prms) )
		gty.pMetaGty[1] = aMetaGty[end]
		evo.pAveFt = 0.0
		evo.pMinF[1] -= 1.0
		# set_tDscdtEtyFactors(ety,gty)
	end
end

# function. evolutionary upgrade
function evoUpgrade!(pop::tEvoPop,evo::tEvoData,maxSize::Integer)
	for i in 1:pop.pN[2]
		if upgradeCondition(pop.aGty[i],maxSize)
			upgradeGtyG!(pop.aGty[i])
			upgradeMetaGty!(pop.ety,pop.aMetaGty,pop.aGty[i],evo)
			fitness!(pop.aGty[i],pop.env)
		end
	end
end

# function: set next evolutionary achievement
function setNextAch!(evo::tEvoData,gen::Integer)
	while evo.avePerformance[gen] <= evo.pMinF[1] + 10^( -evo.pAveFt[1] )
		evo.pAveFt[1] += evo.aveFtinc
	end
end

# function: population evolutionary step: growth, mutation, and selection
# returns: growth factor and number of mutations
function gmsEvoStep!(pop::tEvoPop,elite::Bool=false)
	gf = replication!(pop)
	Nmut = mutation!(pop)
	selection!(pop,elite)

	return gf, Nmut
end

function gmsEvoStep!(pop::tEvoPop,aNichePop::Vector{<:tEvoPop},elite::Bool=false)
	aGf = zeros(Float64,pop.pN[2])			# growth factor array
	aNmut = zeros(Float64,pop.pN[2])		# N mutations array

	# this may be processed in parallel, but for now we do so already downstream
	for (i,nichePop) in enumerate(aNichePop)
		aGf[i], aNmut[i] = gmsEvoStep!(nichePop,elite)
	end

	# sampling among the individuals in each niche to populate the population
	aiSmpl = rand(1:aNichePop[1].pN[2],pop.pN[2])
	pop.aGty[1:pop.pN[2]] .= [ aNichePop[i].aGty[aiSmpl[i]] for i in 1:pop.pN[2] ]

	return sum(aGf)/pop.pN[2], sum(aNmut)/pop.pN[2]
end

# function: individual evolutionary step: growth, mutation, and selection
# returns: growth factor and number of mutations
function gmsOneEvoStep!(pop::tEvoPop,elite::Bool=false)
	growth = replicationOne!(pop)
	Nmut = mutationOne!(pop)
	selectionOne!(pop,elite)

	return growth, Nmut
end

function gmsNicOneEvoStep!(aNicPop::Vector{<:tEvoPop},elite::Bool=false)
	aGf = zeros(Float64,length(aNicPop))		# growth factor array
	aNmut = zeros(Float64,length(aNicPop))		# N mutations array

	@threads for i in eachindex(aNicPop)
		aGf[i], aNmut[i] = gmsOneEvoStep!(aNicPop[i],elite)
	end

	return sum(aGf)/length(aNicPop), sum(aNmut)/length(aNicPop)
end

function metropolisEvoStep!(gty::atGenotype,env::atEnvironment)
	Nmut = 0
	gtyClone = copy(gty)
	for i in 1:length(gty)
		mutation!(gtyClone,rand(THREADRNG[threadid()],1:length(gty)))
		fitness!(gtyClone,env)
		if rand(THREADRNG[threadid()]) < min(1,gtyClone.aF[3]/gty.aF[3])
			copy!(gty,gtyClone)
			Nmut += 1
		else
			copy!(gtyClone,gty)
		end
	end
	return Nmut
end

function metropolisEvoStep!(pop::tEvoPop)
	aNmut = Vector{Int32}(undef,pop.pN[2])
	@threads for i in 1:pop.pN[2]
		aNmut[i] = metropolisEvoStep!(pop.aGty[i],pop.env)
	end
	return sum(aNmut)/pop.pN[2]
end

# Fave + ( 1. - ( Fave % 1 ) ) % ( 1/3 ) + .2

# function: genetic evolution with selection within population niches
function gmsNicED!(pop::tEvoPop,aNichePop::Vector{<:tEvoPop},evo::tEvoData;elite::Bool=false)
	@showprogress 1 "Evolutionary Dynamics Status: " for gen in 1:evo.Ngen
		evo.avePerformance[gen] = mean( [pop.aGty[i].aF[1] for i in 1:pop.pN[2]] )
		if evo.avePerformance[gen] <= evo.pMinF[1] + 10^(- evo.pAveFt[1])
			push!(evo.aEvoPop,deepcopy(pop))		# record the population
			push!(evo.aGen,gen)						# record the generation
			setNextAch!(evo,gen)
			println("\nEvolutionary achievement at generation $gen: ⟨f⟩ = ", evo.avePerformance[gen] )
		end
		evo.growthFactor[gen], evo.mutationFactor[gen] = gmsEvoStep!(pop, aNichePop, elite)
	end
end

# function: genetic evolution with selection within individual niches
function gmsNicOneED!(pop::tEvoPop,evo::tEvoData;elite::Bool=false)
	aNicPop = [ tEvoPop( [Int32(1), Int32(1)], pop.ety, pop.env, pop.aMetaGty, [pop.aGty[i]] ) for i in 1:pop.pN[2] ]

	@showprogress 1 "Evolutionary Dynamics Status: " for gen in 1:evo.Ngen
		evo.avePerformance[gen] = mean( [aNicPop[i].aGty[1].aF[1] for i in 1:pop.pN[2]] )
		if evo.avePerformance[gen] <= evo.pMinF[1] + 10^(- evo.pAveFt[1])
			pop.aGty[1:pop.pN[2]] .= [ aNicPop[i].aGty[1] for i in 1:pop.pN[2] ]
			push!(evo.aEvoPop,deepcopy(pop))
			push!(evo.aGen,gen)
			setNextAch!(evo,gen)
			println("\nEvolutionary achievement at generation $gen: ⟨f⟩ = ", evo.avePerformance[gen] )
		end
		evo.growthFactor[gen], evo.mutationFactor[gen] = gmsNicOneEvoStep!(aNicPop,elite)
	end

	for i in 1:pop.pN[2]
		pop.aGty[i] = aNicPop[i].aGty[1]
	end
end

# function: genetic evolution a la Giardina--Kurchan--Peliti with genetic upgrade
function gmsPopED!(pop::tEvoPop,evo::tEvoData;elite::Bool=false)
	@showprogress 1 "Evolutionary Dynamics Status: " for gen in 1:evo.Ngen
		evo.avePerformance[gen] = mean( [pop.aGty[i].aF[1] for i in 1:pop.pN[2]] )
		if evo.avePerformance[gen] <= evo.pMinF[1] + 10^(- evo.pAveFt[1])
			push!(evo.aEvoPop,deepcopy(pop))
			push!(evo.aGen,gen)
			setNextAch!(evo,gen)
			println("\nEvolutionary achievement at generation $gen: ⟨f⟩ = ", evo.avePerformance[gen] )
		end

		evo.growthFactor[gen], evo.mutationFactor[gen] = gmsEvoStep!(pop,elite)
	end

	return pop, evo
end

# function: genetic evolution a la Giardina--Kurchan--Peliti with genetic upgrade
function gmsPopEDup!(pop::tEvoPop,evo::tEvoData,maxSize::Integer;elite::Bool=false)
	@showprogress 1 "Evolutionary Dynamics Status: " for gen in 1:evo.Ngen
		evo.avePerformance[gen] = mean( [pop.aGty[i].aF[1] for i in 1:pop.pN[2]] )
		if evo.avePerformance[gen] <= evo.pMinF[1] + 10^( -evo.pAveFt[1] )
			push!(evo.aEvoPop,deepcopy(pop))
			push!(evo.aGen,gen)
			setNextAch!(evo,gen)
			println("\nEvolutionary achievement at generation $gen: ⟨f⟩ = ", evo.avePerformance[gen] )
		end

		evo.growthFactor[gen], evo.mutationFactor[gen] = gmsEvoStep!(pop,elite)
		evoUpgrade!(pop,evo,maxSize)
	end
end

function metropolisED!(pop::tEvoPop,evo::tEvoData)
	@showprogress 1 "Evolutionary Dynamics Status: " for gen in 1:evo.Ngen
		evo.avePerformance[gen] = mean( [pop.aGty[i].aF[1] for i in 1:pop.pN[2]] )
		if evo.avePerformance[gen] <= evo.pMinF[1] + 10^(- evo.pAveFt[1])
			push!(evo.aEvoPop,deepcopy(pop))		# record the population
			push!(evo.aGen,gen)						# record the generation
			setNextAch!(evo,gen)
			println("\nEvolutionary achievement at generation $gen: ⟨f⟩ = ", evo.avePerformance[gen] )
		end
		evo.mutationFactor[gen] = metropolisEvoStep!(pop)
	end
end

function generateEvoGty!(pop::tEvoPop,Ngen::Integer;elite::Bool=false)
	for gen in 1:Ngen
		gmsEvoStep!(pop,elite)
	end
	return [ pop.aGty[rand(1:pop.pN[2])] ]
end

export replication!, mutation!, selection!, replicationOne!, mutationOne!, selectionOne!
export gmsEvoStep!, gmsOneEvoStep!, gmsNicOneEvoStep!
export gmsNicED!, gmsNicOneED!, gmsPopED!, gmsPopEDup!, generateEvoGty!
export evoUpgrade!, upgradeGtyG!, setNextAch!

# =============
# |  SYSTEMS  \
# =============

function fitness!(gty::atGenotype,env::atEnvironment)
	gty.aF .= fitness(gty,env)
end

export fitness, fitness!

# ***********
# |  ISING  \
# ***********

# function: flipping --- for metropolis
function flipping!(istMetaGty::tIsingSigTransMetaGty, βJij::Array{<:Real,1}, βhi::Real, n::Array{<:Integer,2})
	# Monte Carlo step evaluated using Glauber transition rates:
	# next possibly transitioning spin (coordinates)
	rng = THREADRNG[threadid()]
	i, j = rand(rng,1:istMetaGty.L), rand(rng,1:istMetaGty.L)
	# transitioning?!
	if i > istMetaGty.li || j > istMetaGty.li
		if rand(rng) < ( 1.0 - n[i,j]*tanh(
				βJij[istMetaGty.Jpi[i,j]]*n[istMetaGty.jp[i],j] + βJij[istMetaGty.Jmi[i,j]]*n[istMetaGty.jm[i],j] +
				βJij[istMetaGty.Jpj[i,j]]*n[i,istMetaGty.jp[j]] + βJij[istMetaGty.Jmj[i,j]]*n[i,istMetaGty.jm[j]] +
				istMetaGty.βhe ))/2.0
			n[i,j] = - n[i,j]
		end
	# else
	# 	if rand(rng) < ( 1. - n[i,j]*tanh(
	# 			βJij[istMetaGty.Jpi[i,j]]*n[istMetaGty.jp[i],j] + βJij[istMetaGty.Jmi[i,j]]*n[istMetaGty.jm[i],j] +
	# 			βJij[istMetaGty.Jpj[i,j]]*n[i,istMetaGty.jp[j]] + βJij[istMetaGty.Jmj[i,j]]*n[i,istMetaGty.jm[j]] +
	# 			istMetaGty.βhe + βhi ))/2
	# 		n[i,j] = - n[i,j]
	# 	end
	end
end

# MonteCarlo simulation function for ising signal transduction: readout magnetization evaluation
function metropolis(istMetaGty::tIsingSigTransMetaGty,Jij::Array{T,1},hi::T)::Real where {T<:Real}
	βJij::Array{Float64,1} = Jij*istMetaGty.β;		βhi::Float64 = hi*istMetaGty.β

	# the initial state vector
	n = Vector{Array{Int8,2}}(undef,nthreads())
	for t in 1:nthreads()
		n[t] = Int8(istMetaGty.he != 0.0 ? sign(istMetaGty.he) : 1).*ones(Int8, istMetaGty.L, istMetaGty.L)
		n[t][1:istMetaGty.li,1:istMetaGty.li] .= sign(hi)
	end

	# definition: time-averaged readout magnetization
	mro = zeros(Float64, nthreads())

	NsmplsThreaded = istMetaGty.prms.Nsmpl ÷ nthreads()
	# istMetaGty.prms.Nsmpl % nthreads() == 0 || throw("mismatch between nthreads and Nsmpl. Choose: Nsmpl = multiple of $nthreads()")
	for t in 1:nthreads()
		for is in 1:NsmplsThreaded
			for imcs in 1:istMetaGty.prms.Nmcsps, ilp in 1:istMetaGty.L2
				flipping!(istMetaGty,βJij,βhi,n[threadid()])
			end
			# evaluation: time-averaged readout magnetization
			mro[threadid()] += sum(n[threadid()][istMetaGty.halfL+1:istMetaGty.halfL+istMetaGty.li,istMetaGty.halfL+1:istMetaGty.halfL+istMetaGty.li])
		end
	end
	# evaluation: time-averaged readout magnetization
	return sum(mro)/(NsmplsThreaded*nthreads()*istMetaGty.li2)
end

# MonteCarlo simulation function for ising signal transduction: time-averaged spin config evaluation
# 	( tIsingST istMetaGty, state n, interaction matrix Jij, input field h ) → readout magnetization m
function metropolis!(istMetaGty::tIsingSigTransMetaGty,Jij::Array{<:Real,1},hi::Real,an::Vector{Array{Int8,2}},prms::tDTMCprm)
	βJij::Array{Float64,1} = Jij.*istMetaGty.β;		βhi::Float64 = hi.*istMetaGty.β

	n = Vector{Array{Int8,2}}(undef,nthreads())
	@threads for t in 1:nthreads()
		n[t] = Int8(istMetaGty.he != 0.0 ? sign(istMetaGty.he) : 1).*ones(Int8, istMetaGty.L, istMetaGty.L)
		n[t][1:istMetaGty.li,1:istMetaGty.li] .= sign(hi)
	end

	# @showprogress 1 "Monte Carlo Dynamics Status: "
	NsmplsThreaded = prms.Nsmpl ÷ nthreads()
	prms.Nsmpl % nthreads() == 0 || throw("mismatch between nthreads and Nsmpl. Choose: Nsmpl = multiple of $nthreads()")
	@threads for t in 1:nthreads()
		for is in 1:NsmplsThreaded
			for imcs in 1:prms.Nmcsps, ilp in 1:istMetaGty.L2
				flipping!(istMetaGty,βJij,βhi,n[threadid()])
			end
			# an[ is + (threadid() - 1)*NsmplsThreaded ] = deepcopy( n[threadid()] )
			an[ is + (t - 1)*NsmplsThreaded ] = deepcopy( n[threadid()] )	# <- when no threading is used
		end
	end
end

# function: fitness for ising signal transduction
function fitness(gty::atSystemGty{<:atIsingMetaGty},env::tCompEnv{<:Vector{<:Vector{<:Real}}})
	ℓVals = zeros(Float64,gty.pMetaGty[1].prms.Ntrials)
	for t in 1:gty.pMetaGty[1].prms.Ntrials
		for io in env.IOidl
			ℓVals[t] += ( metropolis(gty.pMetaGty[1],broadcast(x->10.0^(x),gty.G),io[1]) - io[2] )^2
		end
	end
	d2 = maximum(ℓVals)
	return d2, exp(-d2*env.aSelCoef[1]) + (gty.pMetaGty[1].halfL - BASALFITNESS), exp(-d2*env.aSelCoef[2]) + (gty.pMetaGty[1].halfL - BASALFITNESS)
end

export metropolis, metropolis!

# **************
# |  CHANNELS  \
# **************

# constructor: disordered channel
function tDisChnMetaGty(L::Int32, p::Float64, kout::Float64)
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

	mEvoTypes.tDisChnMetaGty(L, L^2, L2mL, 2Nvbl, Nvbl, p, V, kout)
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

	Λ(q::Vector) = det( Array(W .* ( 1.0 .+ sparse( fill(gty.pMetaGty[1].L2+1,gty.pMetaGty[1].L),
		gty.pMetaGty[1].L2mL+1:gty.pMetaGty[1].L2,
		map(e -> exp(e) - 1.0, q[1:end-1]),gty.pMetaGty[1].L2+1,gty.pMetaGty[1].L2+1 ) ) ) - q[end] .* I )
	dΛ(q::Vector) = ForwardDiff.gradient(Λ,q)

	q0 = zeros(Float64,gty.pMetaGty[1].L+1);	y = dΛ(q0);

	if y[end] != 0.0
		return [ -y[i]/y[end] for i in 1:gty.pMetaGty[1].L ]
	else
		q0[end] = 10^-14;	y = dΛ(q0);
		if y[end] != 0.0
			return [ -y[i]/y[end] for i in 1:gty.pMetaGty[1].L ]
		else
			throw( "a₁(0) = 0. Irreducible stochastic matrix ... or numerical error" )
		end
	end
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

function fitness(gty::atSystemGty{<:atChannelMetaGty},env::tCompEnv{<:Vector{<:Vector{<:Vector{<:Real}}}})
	W = getWnrmd(gty)

	d2::Float64 = 0.0
	for io in env.IOidl
		d2 += sum( (response(gty, W, io[1]) - io[2]).^2 )
	end

	return d2, exp(-d2*env.aSelCoef[1]), exp(-d2*env.aSelCoef[2])		# loss, replication rate function, and selection function
end

export tDisChnMetaGty, response


# *************************
# | STATISTICS FUNCTIONS  \
# *************************

function getGStat!(pop::tEvoPop,Gstat::tStat)
	mapreduce( x -> x ==  pop.aGty[1].pMetaGty[1].dG, &, [ pop.aGty[i].pMetaGty[1].dG for i in 2:pop.pN[2] ] ) ||
		throw(DimensionMismatch("inconsistent dimensions"))

	for x in 1:pop.aGty[1].pMetaGty[1].dG
		Gstat.ave[x] = mean( [pop.aGty[i].G[x] for i in 1:pop.pN[2]] )
	end

	for y in 1:pop.aGty[1].pMetaGty[1].dG, x in 1:pop.aGty[1].pMetaGty[1].dG
		Gstat.cov[x,y] = myCov( [ pop.aGty[i].G[x] for i in 1:pop.pN[2] ], Gstat.ave[x], [ pop.aGty[i].G[y] for i in 1:pop.pN[2] ], Gstat.ave[y] )
	end

	for y in 1:pop.aGty[1].pMetaGty[1].dG, x in 1:pop.aGty[1].pMetaGty[1].dG
		Gstat.cor[x,y] = Gstat.cov[x,y]/sqrt(Gstat.cov[x,x]*Gstat.cov[y,y])
	end
end

export getGStat!

# **********************************
# | STATISTICS FUNCTIONS --- ISING \
# **********************************

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

function getJijStat!(pop::tEvoPop,JijAve::Vector{Float64},JijCov::Array{Float64,2},JijCor::Array{Float64,2})
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

# fills fluxPtrn type with the flux patterns of the genotype in its differnt environments
function getFluxStat!(env::tCompEnv,gty::atSystemGty{<:atChannelMetaGty},fluxPtrn::tFluxPattern)
	length(fluxPtrn.jave) == length(env.IOidl) || throw( "inconsistent dimensions between _fluxPtrn.j_ currents and _env.IOidl_" )
	aCrntCorFisherz = Vector{Array{Float64,2}}(undef,length(env.IOidl))

	W = getW(gty)
	Nq = gty.pMetaGty[1].dG + 2gty.pMetaGty[1].L
	q0 = zeros(Float64, Nq + 1)
	for (i, io) in enumerate(env.IOidl)
		W[1:gty.pMetaGty[1].L,end] = io[1]
		W[end,end] = -sum(W[1:end-1,end])

		Λ(q::Vector) = det( Array(W .* ( 1.0 .+
			sparse( gty.pMetaGty[1].V[1][1:gty.pMetaGty[1].Nvbl], gty.pMetaGty[1].V[2][1:gty.pMetaGty[1].Nvbl],
				map(e -> exp(e) - 1.0, q[1:gty.pMetaGty[1].Nvbl]), gty.pMetaGty[1].L2+1, gty.pMetaGty[1].L2+1 ) .+
			sparse( gty.pMetaGty[1].V[2][1:gty.pMetaGty[1].Nvbl], gty.pMetaGty[1].V[1][1:gty.pMetaGty[1].Nvbl],
				map(e -> exp(e) - 1.0, q[gty.pMetaGty[1].Nvbl+1:gty.pMetaGty[1].dG]), gty.pMetaGty[1].L2+1,
				gty.pMetaGty[1].L2+1) .+
			sparse( fill(gty.pMetaGty[1].L2+1,gty.pMetaGty[1].L), collect(gty.pMetaGty[1].L2mL+1:gty.pMetaGty[1].L2),
				map(e -> exp(e) - 1.0, q[gty.pMetaGty[1].dG+1:gty.pMetaGty[1].dG+gty.pMetaGty[1].L]),
				gty.pMetaGty[1].L2+1, gty.pMetaGty[1].L2+1) .+
			sparse( collect(1:gty.pMetaGty[1].L), fill(gty.pMetaGty[1].L2+1,gty.pMetaGty[1].L),
				map(e -> exp(e) - 1.0, q[gty.pMetaGty[1].dG+gty.pMetaGty[1].L+1:Nq]), gty.pMetaGty[1].L2+1,
				gty.pMetaGty[1].L2+1) ) ) -
			q[end] .* I )
		dΛ(q::Vector) = ForwardDiff.gradient(Λ,q)
		ddΛ(q::Vector) = ForwardDiff.hessian(Λ,q)

		jave = dΛ(q0)
		if jave[end] == 0.0
			q0[end] = 10^-14
			jave = dΛ(q0)
		end
		if jave[end] == 0.0
			throw( "a₁(0) = 0. Irreducible stochastic matrix ... or numerical error" )
		end

		jcov = ddΛ(q0)
		chopArray!(jcov)

		fluxPtrn.jave[i] .= [ -jave[j]/jave[end] for j in 1:Nq ]
		fluxPtrn.jcov[i] .= [ -( jcov[ji,jj] + jcov[ji,end]*fluxPtrn.jave[i][jj] + jcov[jj,end]*fluxPtrn.jave[i][ji] + jcov[end,end]*fluxPtrn.jave[i][ji]*fluxPtrn.jave[i][jj] )/jave[end] for ji in 1:Nq, jj in 1:Nq ]
		aCrntCorFisherz[i] = broadcast( r -> log( (1 + r)/(1 - r) )/2, [ fluxPtrn.jcov[i][ji,jj] / sqrt( fluxPtrn.jcov[i][ji,ji] * fluxPtrn.jcov[i][jj,jj] ) for ji in 1:Nq, jj in 1:Nq ] )

		q0[end] = 0.0
	end

	fluxPtrn.jcor .= [ tanh( mean( [ aCrntCorFisherz[i][ji,jj] for i in eachindex(env.IOidl) ] ) ) for ji in 1:Nq, jj in 1:Nq ]

	fluxPtrn.V[1] .= vcat( gty.pMetaGty[1].V[1], fill(gty.pMetaGty[1].L2+1,gty.pMetaGty[1].L), collect(1:gty.pMetaGty[1].L) )
	fluxPtrn.V[2] .= vcat( gty.pMetaGty[1].V[2], collect(gty.pMetaGty[1].L2mL+1:gty.pMetaGty[1].L2), fill(gty.pMetaGty[1].L2+2,gty.pMetaGty[1].L) )

	return 0
end

function getCrntPtrn!(gty::atSystemGty{<:atChannelMetaGty},fluxPtrn::tFluxPattern,crntPtrn::tCrntPattern)
	Ncrnt = Int32(gty.pMetaGty[1].dG/2)
	for (io,jave) in enumerate(fluxPtrn.jave)
		for i in 1:Ncrnt
			if jave[i] > jave[i + Ncrnt]
				crntPtrn.aV[io][1][i] = fluxPtrn.V[1][i]
				crntPtrn.aV[io][2][i] = fluxPtrn.V[2][i]
				crntPtrn.Jave[io][i] = jave[i] - jave[i + Ncrnt]
			else
				crntPtrn.aV[io][1][i] = fluxPtrn.V[2][i]
				crntPtrn.aV[io][2][i] = fluxPtrn.V[1][i]
				crntPtrn.Jave[io][i] = jave[i + Ncrnt] - jave[i]
			end
		end
		crntPtrn.aV[io][1][Ncrnt+1:Ncrnt+2gty.pMetaGty[1].L] .= fluxPtrn.V[1][gty.pMetaGty[1].dG+1:gty.pMetaGty[1].dG+2gty.pMetaGty[1].L]
		crntPtrn.aV[io][2][Ncrnt+1:Ncrnt+2gty.pMetaGty[1].L] .= fluxPtrn.V[2][gty.pMetaGty[1].dG+1:gty.pMetaGty[1].dG+2gty.pMetaGty[1].L]
		crntPtrn.Jave[io][Ncrnt+1:Ncrnt+2gty.pMetaGty[1].L] = jave[gty.pMetaGty[1].dG+1:gty.pMetaGty[1].dG+2gty.pMetaGty[1].L]
	end
end

function response!(gty::atSystemGty{<:atChannelMetaGty},W::AbstractMatrix,Input::Vector{<:Real},p)
	length(Input) == gty.pMetaGty[1].L || throw( "inconsistent dimensions between input currents and system" )
	W[1:gty.pMetaGty[1].L,end] = Input
	b = zeros(Float64,gty.pMetaGty[1].L2+1);	b[end] = 1.0

	p .= ( W \ b )[1:end-1]
end

function collectResponses!(pop::tEvoPop{<:atEvotype,<:atEnvironment,<:Vector{<:atChannelMetaGty},<:Vector{<:atGenotype}},aResp::AbstractArray)
	i = 0
	for gty in pop.aGty[1:pop.pN[2]]
		W = getWnrmd(gty)
		for io in pop.env.IOidl
			response!(gty,W,io[1],view(aResp,i+=1,:))
		end
	end
end

function getRStat!(pop::tEvoPop{<:atEvotype,<:atEnvironment,<:Vector{<:atChannelMetaGty},<:Vector{<:atGenotype}},Rstat::tStat)

	aResp = Array{Float64}(undef, pop.pN[2]*length(pop.env.IOidl), pop.aMetaGty[1].L2)
	collectResponses!(pop,aResp)

	for x in 1:pop.aMetaGty[1].L2
		Rstat.ave[x] = mean( aResp[:,x] )
	end

	for y in 1:pop.aMetaGty[1].L2, x in 1:pop.aMetaGty[1].L2
		Rstat.cov[x,y] = myCov( aResp[:,x], Rstat.ave[x], aResp[:,y], Rstat.ave[y] )
	end

	for y in 1:pop.aMetaGty[1].L2, x in 1:pop.aMetaGty[1].L2
		Rstat.cor[x,y] = Rstat.cov[x,y]/sqrt(Rstat.cov[x,x]*Rstat.cov[y,y])
	end

	return aResp
end

export getFluxStat!, response!, collectResponses!, getRStat!

# *********************************
# | ROBUSTNESS ANALYSIS FUNCTIONS \
# *********************************

sensitivity(ℓm::Real,ℓ0::Real) = (ℓm-ℓ0)/ℓ0

# function: distribution of sensitivity for a given genome in a given environment
function testRbst!(gty::atGenotype,env::atEnvironment,aRbst::Vector{Float64})
	length(gty) == length(aRbst) || throw(DimensionMismatch("sensitivity vector's length must be the same as gty's"))

	gtyClone = copy(gty)
	for i in eachindex(aRbst)
		copy!(gtyClone,gty)
		mutation!(gtyClone,i)
		fitness!(gtyClone,env)
		aRbst[i] = sensitivity(gtyClone.aF[1],gty.aF[1])
	end
end

# function: distribution of sensitivity for the population
function testRbst!(pop::tEvoPop,aRbst::Vector{Vector{Float64}})
	pop.pN[2] == length(aRbst) || throw(DimensionMismatch("sensitivity array's length must be the same as pop.pN[2]"))

	# @showprogress 1 "Robustness Test Status: "
	@threads for i in eachindex(aRbst)
		testRbst!(pop.aGty[i],pop.env,aRbst[i])
	end
end

function testRbst!(pop::tEvoPop,aRbst::Vector{Vector{Float64}},prfBins::Vector{Float64})
	length(prfBins)-1 == length(aRbst) || throw(DimensionMismatch("sensitivity array's length must be the same as prfBins-1"))

	# @showprogress 1 "Robustness Test Status: "
	@threads for gty in pop.aGty[1:pop.pN[2]]
		iBin = get_binIndex(gty.aF[1],prfBins)
		if iBin > 0
			vRbst = Vector{Float64}(undef,length(gty))
			testRbst!(gty,pop.env,vRbst)
			append!(aRbst[iBin],vRbst)
		end
	end
end

# fucntion. array of vector of equal length --> array with their Hamming distances
function hamming(a::Array{<:AbstractArray})
	aH, iaH = Vector{Float64}(undef, Int64(( length(a) * (length(a) - 1) ) / 2) ), 0

	for i in eachindex(a), j in i+1:length(a)
		aH[iaH+=1] = Distances.hamming(a[i], a[j])
	end

	return aH
end

# function. population ⊗ performance bin --> [ array of hamming distances for each bin ]
function hamming(pop::tEvoPop,prfBins::Vector{Float64})
	aaH = Vector{Vector{Int64}}(undef, length(prfBins) - 1)

	aBinIndex = [ get_binIndex(pop.aGty[i].aF[1],prfBins) for i in 1:pop.pN[2] ]
	@threads for iBin in eachindex(aaH)
		aGtyIndex = findall( e -> e == iBin, aBinIndex )
		aaH[iBin] = hamming( [ pop.aGty[i].G for i in aGtyIndex ] )
	end

	return aaH
end

# # move to Utils
# function get_binIndex(val::T,aBins::Vector{<:T}) where T
# 	for i in 1:length(aBins)-2
# 		if val >= aBins[i] && val < aBins[i+1]
# 			return i
# 		end
# 	end
# 	if val >= aBins[end-1] && val <= aBins[end]
# 		return length(aBins)-1
# 	else
# 		return 0
# 	end
# end

export sensitivity, testRbst!

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
			print(f, aGty[i].aF[2], "\n")
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
			print(f, aGty[i].aF[2], "\n")
		end
	end
	open( fileTag * "_a" * split(string(typeof(aGty[1])),"{")[1] * ".dat", "w" ) do f
		print(f, length(aGty[i].g), "\t")
		for x in aGty[1].g
			print(f, x, "\t")
		end
	end
end

# function: saving the ising population
function write_tEvoPop(pop::tEvoPop{<:atEvotype,<:atEnvironment,<:Vector{<:atIsingMetaGty},<:Vector{<:atGenotype}},fileTag::String)
	write_DTMCprm(pop.aGty[1].pMetaGty[1].prms,fileTag)
	write_MetaGty(pop.aGty[1].pMetaGty[1],fileTag)
	write_aGty(pop.aGty,pop.pN[2],fileTag)
end

# function: saving the evolutionary data
function write_tEvoData(aData::Vector{tEvoData},fileName::String)
	open( fileName, "w" ) do f
		for dataBatch in aData, t in 1:dataBatch.Ngen print(f,dataBatch.avePerformance[t],"\t") end, print(f,"\n")
		for dataBatch in aData, t in 1:dataBatch.Ngen print(f,dataBatch.growthFactor[t],"\t") end, print(f,"\n")
		for dataBatch in aData, t in 1:dataBatch.Ngen print(f,dataBatch.mutationFactor[t],"\t") end, print(f,"\n")
	end

	# for dataBatch in aData
	# 	for (gen, pop)  in zip(dataBatch.aGen, dataBatch.aEvoPop)
	# 		write_tEvoPop( pop, fileTag * "_g$gen")
	# 	end
	# end
end

# to do: saving evolution parameters

export write_DTMCprm, write_MetaGty, write_aGty, write_tEvoPop, write_tEvoData

function read_tEvoData(fileName::String)
	dataMat = readdlm( fileName )
	return dataMat[1,:], dataMat[2,:], dataMat[3,:]
end

# read metagenotype
function read_MetaGty(L::Integer,prms::atMonteCarloPrm,fileTag::String)::atIsingMetaGty
	metaGtyMat = readdlm( "data/" * fileTag * "_MetaGty" * ".dat" )
	return tIsingSigTransMetaGty(L,metaGtyMat[1],metaGtyMat[2],prms)
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

export read_aIsingSigTransGty, read_tEvoData

# *******************
# TRIVIAL
# *******************

function fitness!(gty::atGenotype,trivialEnv::tTrivialEnv)
	gty.aF[2]=1/(euclidean(gty.G,ones(Float64,gty.pMetaGty[1].dG))+FITNESSOFFSET)
end

export fitness!

end

# to do:
# * relax Int32 constraint from functions whenever possible;
# * generation function randomly generating its initial population
# * robustness anaalysis from cluster
# * everything from jupyter
# * check the representation of integers and Int32
