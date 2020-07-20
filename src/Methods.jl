
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

# generic population replication
function replication!(pop::AbstractPopulation)
	Kr = Vector{Float64}(undef,nthreads())
	ipKr = Vector{Int32}(undef,nthreads())
	G = zeros(Int32,pop.pN[2])

	# @threads
	for i in 1:pop.pN[2]
		Kr[threadid()] = pop.ety.minRepCoef + pop.ety.pRepFactor[1]*pop.aGty[i].aFitness[2]
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

function replication!(pop::AbstractPopulation, ancestry)
	Kr = Vector{Float64}(undef,nthreads())
	ipKr = Vector{Int32}(undef,nthreads())
	G = zeros(Int32,pop.pN[2])

	# @threads
	for i in 1:pop.pN[2]
		Kr[threadid()] = pop.ety.minRepCoef + pop.ety.pRepFactor[1]*pop.aGty[i].aFitness[2]
		ipKr[threadid()] = trunc(Int32,Kr[threadid()])
		G[i] = rand(THREADRNG[threadid()]) < Kr[threadid()] - ipKr[threadid()] ? ipKr[threadid()] + 1 : ipKr[threadid()]
	end

	for i in 1:pop.pN[2]
		for inew in 1:G[i]
			pop.pN[1] += 1
			if pop.pN[1] <= length(pop.aGty)
				pop.aGty[pop.pN[1]] = copy(pop.aGty[i])
				ancestry[pop.pN[1]] = i
			else
				push!(pop.aGty,copy(pop.aGty[i]))
				push!(ancestry,i)
			end
		end
	end

	return log(pop.pN[1]/pop.pN[2])
end


function mutation!(pop::AbstractPopulation)
	aNmut = Vector{Int32}(undef,pop.pN[1])
	# @threads
	for i in 1:pop.pN[1]
		aNmut[i] = mutation!(pop.aGty[i], pop.ety, pop.env)
		# aNmut[i] = EvolutionaryDynamics.mutation!(pop.aGty[i], pop.ety, pop.env)
	end
	return sum(aNmut)/pop.pN[2] 	# normalized number of mutations
end


selection!(pop::AbstractPopulation) = _selection(pop.ety.selType, pop)

function _selection!(::FitnessSelection, pop::AbstractPopulation)
	aGtyRef::Vector{AbstractGenotype} = copy(pop.aGty)
	# Here, you are copying the references of the genotypes. It is fine, because when you select the new generation,
	# you are not changing the genotypes themeselves, but simply picking new references.

	aSelFncVals = [ fitness(pop.aGty[i]) for i in 1:pop.pN[1] ] / sum([ fitness(pop.aGty[i]) for i in 1:pop.pN[1] ])
	aiSelected = rand(THREADRNG[threadid()], Categorical(aSelFncVals), pop.pN[2])

	for i in 1:pop.pN[2]
		pop.aGty[i] = aGtyRef[aiSelected[i]]
	end

	# population size renormalization
	pop.pN[1] = pop.pN[2]

	return aiSelected[1]
end

evolveEnvironment!(env::Tenv) where Tenv <: AbstractEnvironment = _evolveEnvironment!(IsVaryingEnvironment(Tenv), env)

_evolveEnvironment!(::IsVaryingEnvironment, env) = error("environmental evolution not defined")

function _evolveEnvironment!(::MarkovianEnvironment, env)
	d = Categorical( env.transMtx[:, envState(env)] )
	i = rand(THREADRNG[threadid()], d)

	if envState(env) != i
		env.envState[1] = i
		return true
	else
		return false
	end
end

export replication!, mutation!, selection!, _selection!

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
