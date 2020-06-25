
# generic population replication
function replication!(pop::AbstractPopulation)
	Kr = Vector{Float64}(undef,nthreads())
	ipKr = Vector{Int32}(undef,nthreads())
	G = zeros(Int32,pop.pN[2])

	# @threads
	for i in 1:pop.pN[2]
		Kr[threadid()] = pop.ety.minRepCoef + pop.ety.pRepFactor[1]*pop.aGty[i].aF[2]
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

function mutation!(pop::AbstractPopulation)
	aNmut = Vector{Int32}(undef,pop.pN[1])
	# @threads
	for i in 1:pop.pN[1]
		aNmut[i] = mutation!(pop.aGty[i],pop.ety,pop.env)
	end
	return sum(aNmut)/pop.pN[2] 	# normalized number of mutations
end


import Distributions.Categorical

function selection!(pop::AbstractPopulation,elite::Bool=false)
	popGtyRef::Array{AbstractGenotype,1} = copy(pop.aGty)

	if elite
		# survival of the fittest
		aiSelected = sortperm([ pop.aGty[i].aF[3] for i in 1:pop.pN[1] ],rev=true)
		# sort!(pop.aGty, by= x -> x.aF[2], rev=true)
	else
		# selection with replacement of individuals based on selection function values
		aSelFncVals = [ pop.aGty[i].aF[3] for i in 1:pop.pN[1] ] / sum([ pop.aGty[i].aF[3] for i in 1:pop.pN[1] ])
		d = Categorical(aSelFncVals)
		aiSelected = rand(THREADRNG[threadid()],d,pop.pN[2])
	end

	for i in 1:pop.pN[2]
		pop.aGty[i] = popGtyRef[aiSelected[i]]
	end

	# population size renormalization
	pop.pN[1] = pop.pN[2]
end
