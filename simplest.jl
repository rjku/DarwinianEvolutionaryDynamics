
using Random, Distances, Plots

const NGEN, NPOP, DGTY, REPRATE, MUTRATE  = 100000, 10000, 5, 1.0, 0.1

# fitness function
fitness(gty::Array{Float64,1}, Dgty::Int64)::Float64 = 1/(euclidean(gty,ones(Float64,Dgty))+.1)

# replication function: ( population of genotypes, initial population size, dimension genotypic space ) → replicated population
function replication!(popGty::Array{Array{Float64,1},1}, fitPopGty::Array{Float64,1}, pNpopGty::Array{Int64,1}, repFactor::Float64)
	for i in 1:pNpopGty[2]
		Kr::Float64 = repFactor*fitPopGty[i]
		ipKr::Int8 = trunc(Int32,Kr)

		G::Int8 = rand() < Kr - ipKr ? ipKr + 1 : ipKr

		for inew in 1:G
			pNpopGty[1] += 1
			if pNpopGty[1] <= length(popGty)
				popGty[pNpopGty[1]] = deepcopy(popGty[i])
			else
				push!(popGty,popGty[i])
				push!(fitPopGty,fitPopGty[i])
			end
		end
	end
end

# mutation function: ( genotype, dimension genotypic space ) → mutated genotype
function effMutation!(popGty::Array{Array{Float64,1},1}, fitPopGty::Array{Float64,1}, pNpopGty::Array{Int64,1}, Dgty::Int64, repFactor::Float64, mutProb::Float64)
	Δx::Float32 = .01

	for i in 1:pNpopGty[1]
		r1::Float32 = rand()
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

function effSelection!(popGty::Array{Array{Float64,1},1}, pNpopGty::Array{Int64,1})

	growthFactor::Float64 = log(pNpopGty[1]/pNpopGty[2])
	popGtyRef = copy(popGty)

	# selection with replacement of individuals
	survivedGty = rand(collect(1:pNpopGty[1]),pNpopGty[2])
	for i in 1:pNpopGty[2]
		popGty[i] = popGtyRef[survivedGty[i]]
	end

	# population size renormalization
	pNpopGty[1] = pNpopGty[2]

	return growthFactor
end

# main function
function main()#::Int8
	# variable Npopulation as a "pointer"
	pNpop = Int64[NPOP,NPOP]

	# ∀ individual → genotypic variable
	# the factor 2 accounts for the maximum population growth in one step
	𝕏 = [ zeros(Float64,DGTY) for i in 1:2*pNpop[1] ]

	# ∀ individual → fitness value
	𝔽 = Float64[ i <= pNpop[1] ? fitness(𝕏[i],DGTY) : 0.0 for i in 1:2*pNpop[1] ]

	# growth factors
	μ = zeros(Float64,NGEN)

	println("\n", 𝕏[1:pNpop[1]])
	# println("\n", 𝕏, " is a: ", typeof(𝕏))
	# println(𝔽, " is a: ", typeof(𝔽))


	# evolution time scale, mutation probability, and replication factor
	Δt::Float64 = 1/(2*DGTY*MUTRATE+.01)
	κm::Float64 = MUTRATE*Δt
	κr::Float64 = REPRATE*Δt

	for gen in 1:NGEN
		replication!(𝕏,𝔽,pNpop,κr)
		effMutation!(𝕏,𝔽,pNpop,DGTY,κr,κm)
		μ[gen] = effSelection!(𝕏,pNpop)
	end

	println(𝕏[1:pNpop[1]])

	pl = plot(collect(1:NGEN),μ)

	return μ
end

mu = main()

pyplot() # Switch to using the PyPlot.jl backend
plot(collect(1:NGEN),mu)
