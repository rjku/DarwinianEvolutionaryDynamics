
module Methods

# using ..EvolutionaryDynamics
using ..Types

using Base.Threads, Random
import Future

# threads random generators initialization
const THREADRNG = let m = MersenneTwister(1)
            [m; accumulate(Future.randjump, fill(big(10)^20, nthreads()-1), init=m)]
        end;

import Distributions.Categorical

#---------
# Ancestry

ancestry(ety::AbstractEvotype, Npop::Integer) = _ancestry(ety.repType, Npop)

function _ancestry(repType::NeutralReplication, Npop)
	ancestryVec = [ i for i in 1:((repType.repCoef + 1) * Npop) ]
	iAncestry = Npop
	for i in 1:Npop
		for inew in 1:repType.repCoef
			iAncestry += 1
			ancestryVec[iAncestry] = i
		end
	end

	return ancestryVec
end # function _ancestry NeutralReplication

_ancestry(::FitnessReplication, Npop) = [ i for i in 1:Npop ]

_ancestry(::WithoutReplication, Npop) = [ i for i in 1:Npop ]

export ancestry
#------------
# Replication

replication!(pop::AbstractPopulation) = _replication!(pop.ety.repType, pop)

function _replication!(repType::FitnessReplication, pop::AbstractPopulation)
	Kr = Vector{Float64}(undef,	nthreads())
	ipKr = Vector{Int32}(undef,	nthreads())
	G = zeros(Int32, pop.pN[2])

	@threads for i in 1:pop.pN[2]
		Kr[threadid()] = repType.repFactor * fitness(pop.aGty[i])
		ipKr[threadid()] = trunc(Int32, Kr[threadid()])
		# G[i] = rand(THREADRNG[threadid()]) < Kr[threadid()] - ipKr[threadid()] ? ipKr[threadid()] + 1 : ipKr[threadid()]
		G[i] = rand() < Kr[threadid()] - ipKr[threadid()] ? ipKr[threadid()] + 1 : ipKr[threadid()]
	end

	for i in 1:pop.pN[2]
		for inew in 1:G[i]
			pop.pN[1] += 1
			if pop.pN[1] <= length(pop.aGty)
				pop.aGty[pop.pN[1]] = copy(pop.aGty[i])
			else
				push!(pop.aGty, copy(pop.aGty[i]))
			end
		end
	end

	return log(pop.pN[1] / pop.pN[2])
end

function _replication!(repType::NeutralReplication, pop)
	for i in 1:pop.pN[2]
		for inew in 1:repType.repCoef
			pop.pN[1] += 1
			if pop.pN[1] <= length(pop.aGty)
				pop.aGty[pop.pN[1]] = copy(pop.aGty[i])
			else
				push!(pop.aGty,copy(pop.aGty[i]))
			end
		end
	end

	return log(pop.pN[1] / pop.pN[2])
end # function _replication! NeutralReplication

_replication!(repType::WithoutReplication, pop) = 0.0


replication!(pop::AbstractPopulation, ancestryVec::Vector{<:Integer}) = _replication!(pop.ety.repType, pop, ancestryVec)

_replication!(repType::ReplicationType, pop, ancestryVec) = _replication!(repType, pop)

function _replication!(repType::FitnessReplication, pop, ancestryVec)
	Kr = Vector{Float64}(undef,	nthreads())
	ipKr = Vector{Int32}(undef,	nthreads())
	G = zeros(Int32, pop.pN[2])

	@threads for i in 1:pop.pN[2]
		Kr[threadid()] = repType.repFactor * fitness(pop.aGty[i])
		ipKr[threadid()] = trunc(Int32, Kr[threadid()])
		# G[i] = rand(THREADRNG[threadid()]) < Kr[threadid()] - ipKr[threadid()] ? ipKr[threadid()] + 1 : ipKr[threadid()]
		G[i] = rand() < Kr[threadid()] - ipKr[threadid()] ? ipKr[threadid()] + 1 : ipKr[threadid()]
	end

	for i in 1:pop.pN[2]
		for inew in 1:G[i]
			pop.pN[1] += 1
			if pop.pN[1] <= length(pop.aGty)
				pop.aGty[pop.pN[1]] = copy(pop.aGty[i])
				ancestryVec[pop.pN[1]] = i
			else
				push!(pop.aGty, copy(pop.aGty[i]))
				push!(ancestryVec, i)
			end
		end
	end

	return log(pop.pN[1] / pop.pN[2])
end

export replication!
#---------
# Mutation

function mutation!(pop::AbstractPopulation)
	aNmut = Vector{Int32}(undef,pop.pN[1])

	@threads for i in 1:pop.pN[1]
		aNmut[i] = mutation!(pop.aGty[i], pop.ety)

		if aNmut[i] > 0
			fitness!(pop.aGty[i], pop.env)
		end
	end
	return sum(aNmut) / pop.pN[2] 	# normalized number of mutations
end # function mutation!

mutation!(gty::AbstractGenotype, ety::AbstractEvotype) = _mutation!(ety.mutType, gty, ety)

_mutation!(::MutationType, gty::AbstractGenotype, ety::AbstractEvotype) = 0

export mutation!
#----------
# Selection

selection!(pop::AbstractPopulation) = _selection!(pop.ety.selType, pop)

function _selection!(::NeutralSelection, pop::AbstractPopulation)
	aGtyRef::Vector{AbstractGenotype} = copy(pop.aGty)

	# aiSelected = rand(THREADRNG[threadid()], 1:pop.pN[1], pop.pN[2])
	aiSelected = rand(1:pop.pN[1], pop.pN[2])

	for i in 1:pop.pN[2]
		pop.aGty[i] = deepcopy(aGtyRef[aiSelected[i]])
	end

	# population size renormalization
	pop.pN[1] = pop.pN[2]

	return aiSelected[1]
end # function _selection! NeutralSelection

function _selection!(::FitnessSelection, pop::AbstractPopulation)
	aGtyRef::Vector{AbstractGenotype} = copy(pop.aGty)

	aSelFncVals = [ fitness(pop.aGty[i]) for i in 1:pop.pN[1] ] / sum([ fitness(pop.aGty[i]) for i in 1:pop.pN[1] ])
	# aiSelected = rand(THREADRNG[threadid()], Categorical(aSelFncVals), pop.pN[2])
	aiSelected = rand(Categorical(aSelFncVals), pop.pN[2])

	for i in 1:pop.pN[2]
		pop.aGty[i] = deepcopy(aGtyRef[aiSelected[i]])
	end

	# population size renormalization
	pop.pN[1] = pop.pN[2]

	return aiSelected[1]
end # function _selection! FitnessSelection

function _selection!(::ElitismSelection, pop::AbstractPopulation)
	aGtyRef::Vector{AbstractGenotype} = copy(pop.aGty)

 	aiSelected = sortperm([ fitness(pop.aGty[i]) for i in 1:pop.pN[1] ], rev=true)

	for i in 1:pop.pN[2]
		pop.aGty[i] = aGtyRef[aiSelected[i]]
	end

	# population size renormalization
	pop.pN[1] = pop.pN[2]

	return aiSelected[1]	# the fittest of fittests
end # function _selection! ElitismSelection

export selection!
#--------
# Fitness

fitness!(gty::AbstractGenotype, env::Tenv) where Tenv <: AbstractEnvironment = _fitness!(IsVaryingEnvironment(Tenv), gty, env)

function fitness!(pop::AbstractPopulation)
	for i in 1:pop.pN[2]
		fitness!(pop.aGty[i], pop.env)
	end
end

_fitness!(::IsVaryingEnvironment, gty::AbstractGenotype, env::AbstractEnvironment) = error("fitness function undefined")

export fitness!
#------------
# Environment

evolveEnvironment!(env::Tenv) where Tenv <: AbstractEnvironment = _evolveEnvironment!(IsVaryingEnvironment(Tenv), env)

_evolveEnvironment!(::IsVaryingEnvironment, env) = error("environmental evolution not defined")

function _evolveEnvironment!(::MarkovianEnvironment, env)
	d = Categorical( env.transMtx[:, envState(env)] )
	# i = rand(THREADRNG[threadid()], d)
	i = rand(d)

	if envState(env) != i
		env.envState[1] = i
		return true
	else
		return false
	end
end

export evolveEnvironment!
#----------
# utilities

function draw(probVec::Vector{<:Real}, r::Real)
	if r <= sum(probVec)
		cumProb, i = probVec[1], 1
		while cumProb < r
			cumProb += probVec[i += 1]
		end
		return i
	else
		return 0
	end
end

export draw

end  # module Methods

# function selection!(pop::AbstractPopulation)
# 	aGtyRef::Array{AbstractGenotype,1} = copy(pop.aGty)
# 	# here you are copying the references of the genotypes
# 	# it is fine, because when you select the new generation you are not changing the genotypes themeselves, but simply picking new references
#
# 	# selection with replacement of individuals based on selection function values
# 	aSelFncVals = [ pop.aGty[i].aF[3] for i in 1:pop.pN[1] ] / sum([ pop.aGty[i].aF[3] for i in 1:pop.pN[1] ])
# 	aiSelected = rand(THREADRNG[threadid()], Categorical(aSelFncVals), pop.pN[2])
#
# 	# if you wanted elite selection
# 	# aiSelected = sortperm([ pop.aGty[i].aF[3] for i in 1:pop.pN[1] ],rev=true)
#
# 	for i in 1:pop.pN[2]
# 		pop.aGty[i] = aGtyRef[aiSelected[i]]
# 	end
#
# 	# population size renormalization
# 	pop.pN[1] = pop.pN[2]
#
# 	return aiSelected[1]
# end
