
using Random, Distances, Plots, Base.Threads
import Future

const NGEN, NPOP, DGTY = 1000, 1000, 2
const REPRATE, MUTRATE = 10., 0.5
const FITNESSOFFSET, DELTATOFFSET =  1.0, 1.0

# fitness function
fitness(gty::Array{Float64,1}, Dgty::Int64)::Float64 = 1/(euclidean(gty,ones(Float64,Dgty))+FITNESSOFFSET)

# replication function: ( population of genotypes, initial population size, dimension genotypic space ) → replicated population
function replication!(popGty::Array{Array{Float64,1},1}, fitPopGty::Array{Float64,1}, pNpopGty::Array{Int64,1}, repFactor::Float64, R::Vector{MersenneTwister})
	G = zeros(Int64,pNpopGty[2])
	@threads for i in 1:pNpopGty[2]
		Kr::Float64 = repFactor*fitPopGty[i]
		ipKr::Int8 = trunc(Int32,Kr)

		G[i] = rand(R[threadid()]) < Kr - ipKr ? ipKr + 1 : ipKr
		# G::Int64 = rand() < Kr - ipKr ? ipKr + 1 : ipKr
	end

	for i in 1:pNpopGty[2]
		for inew in 1:G[i]
			pNpopGty[1] += 1
			if pNpopGty[1] <= pNpopGty[3]
				popGty[pNpopGty[1]] = deepcopy(popGty[i])
				fitPopGty[pNpopGty[1]] = fitPopGty[i]
			else
				pNpopGty[3] += 1
				push!(popGty,popGty[i])
				push!(fitPopGty,fitPopGty[i])
			end
		end
	end
end

# mutation function: ( genotype, fitness, dimension genotypic space ) → mutated genotype
function effMutation!(popGty::Array{Array{Float64,1},1}, fitPopGty::Array{Float64,1}, pNpopGty::Array{Int64,1}, Dgty::Int64, repFactor::Float64, mutProb::Float64, R::Vector{MersenneTwister})
	Δx::Float32 = .01

	@threads for i in 1:pNpopGty[1]
		r1::Float32 = rand(R[threadid()])
		pCum::Float32 = 0.
		xvar::Int16 = -1

		effMutProb::Float64 = mutProb/(1+repFactor*fitPopGty[i])

		while pCum < r1 && xvar < 2*Dgty
			xvar += 1
			pCum += effMutProb
		end
		if xvar < 2*Dgty
			popGty[i][xvar%Dgty+1] += xvar < Dgty ? Δx : -Δx
			fitPopGty[i] = fitness(popGty[i],Dgty)
		end
	end
end

# effective selection function: ( genotype populations, fitness, population size ) → selected population
function effSelection!(popGty::Array{Array{Float64,1},1}, fitPopGty::Array{Float64,1}, pNpopGty::Array{Int64,1}; ubermode::Bool=false)
	growthFactor::Float64 = log(pNpopGty[1]/pNpopGty[2])
	popGtyRef, fitPopGtyRef = copy(popGty), copy(fitPopGty)

	if ubermode
		# survival of the fittest
		survivedGty = sortperm(fitPopGty,rev=true)
	else
		# selection with replacement of individuals
		survivedGty = rand(collect(1:pNpopGty[1]),pNpopGty[2])
	end
	for i in 1:pNpopGty[2]
		popGty[i] = popGtyRef[survivedGty[i]]
		fitPopGty[i] = fitPopGtyRef[survivedGty[i]]
	end

	# population size renormalization
	pNpopGty[1] = pNpopGty[2]

	return growthFactor
end

# ==============================

# main function
function main()#::Int8
	R = let m = MersenneTwister(1)
	        [m; accumulate(Future.randjump, fill(big(10)^20, nthreads()-1), init=m)]
	    end;

	# variable Npopulation as a "pointer" : effective population, fixed population value, array size
	pNpop = Int64[NPOP,NPOP,3*NPOP]

	# ∀ individual → genotypic variable
	# the factor 2 accounts for the maximum population growth in one step
	𝕏 = [ zeros(Float64,DGTY) for i in 1:pNpop[3] ]

	# ∀ individual → fitness value
	𝔽 = Float64[ i <= pNpop[1] ? fitness(𝕏[i],DGTY) : 0.0 for i in 1:pNpop[3] ]

	# growth factors and average fitness
	μ = zeros(Float64,NGEN)
	ϕ = zeros(Float64,NGEN)

	# println("\n", 𝕏[1:pNpop[1]])
	# println("\n", 𝕏, " is a: ", typeof(𝕏))
	# println(𝔽, " is a: ", typeof(𝔽))


	# evolution time scale, mutation probability, and replication factor
	Δt::Float64 = 1/(2*DGTY*MUTRATE+DELTATOFFSET)
	κm::Float64 = MUTRATE*Δt
	κr::Float64 = REPRATE*Δt

	for gen in 1:NGEN
		replication!(𝕏,𝔽,pNpop,κr,R)
		effMutation!(𝕏,𝔽,pNpop,DGTY,κr,κm,R)
		μ[gen] = effSelection!(𝕏,𝔽,pNpop,ubermode=true)
		ϕ[gen] = sum(𝔽[1:pNpop[2]])/pNpop[2]
	end

	# println(𝕏[1:pNpop[1]])

	pyplot() # Switch to using the PyPlot.jl backend
	return plot(collect(1:NGEN),ϕ)
	# return plot(collect(1:NGEN),μ)
end

# main()
