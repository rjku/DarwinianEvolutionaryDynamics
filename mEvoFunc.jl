
module mEvoFunc
using mEvoTypes, Random, Base.Threads, Distances
import Future


# *******************
# PROPER EVOLUTIONARY FUNCTIONS
# *******************

# replication function: ( living population, replication factor, Marsenne Twister ) → replicated population
function replication!(pop::tLivingPop,R::Vector{MersenneTwister})
	G = zeros(Int32,pop.pN[2])
	@threads for i in 1:pop.pN[2]
		Kr::Float64 = pop.repFactor*pop.aGty[i].pF[1]
		ipKr::Int32 = trunc(Int32,Kr)

		G[i] = rand(R[threadid()]) < Kr - ipKr ? ipKr + 1 : ipKr
		# G::Int64 = rand() < Kr - ipKr ? ipKr + 1 : ipKr
	end

	for i in 1:pop.pN[2]
		for inew in 1:G[i]
			pop.pN[1] += 1
			if pop.pN[1] <= pop.pN[3]
				pop.aGty[pop.pN[1]] = deepcopy(pop.aGty[i])
			else
				pop.pN[3] += 1
				push!(pop.aGty,pop.aGty[i])
			end
		end
	end
end

# effective selection function: ( population, MT ) → selected population
function effSelection!(pop::tLivingPop; ubermode::Bool=false)
	growthFactor::Float64 = log(pop.pN[1]/pop.pN[2])
	popGtyRef::Array{atGenotype,1} = deepcopy(pop.aGty)

	if ubermode
		# survival of the fittest
		survivedGty = sortperm([pop.aGty[i].pF[1] for i in 1:pop.pN[1]],rev=true)
	else
		# selection with replacement of individuals
		survivedGty = rand(collect(1:pop.pN[1]),pop.pN[2])
	end
	for i in 1:pop.pN[2]
		pop.aGty[i] = popGtyRef[survivedGty[i]]		# error?
	end

	# population size renormalization
	pop.pN[1] = pop.pN[2]

	return growthFactor
end

# effective mutation function: ( population, variation, fitness!Function acting on genotypes, MT )
function effMutation!(pop::tLivingPop, R::Vector{MersenneTwister})
	@threads for i in 1:pop.pN[1]
		r1::Float64 = rand(R[threadid()])
		cumProb::Float64 = 0.
		xvar::Int32 = -1

		effMutProb::Float64 = pop.mutFactor/(1+pop.repFactor*pop.aGty[i].pF[1])

		while cumProb < r1 && xvar < 2pop.aGty[i].pdX[1]
			xvar += 1
			cumProb += effMutProb
		end
		if xvar < 2pop.aGty[i].pdX[1]
			pop.aGty[i].X[xvar%pop.aGty[i].pdX[1]+1] += xvar < pop.aGty[i].pdX[1] ? pop.Xvar : -pop.Xvar
			fitness!(pop.ety,pop.env,pop.aGty[i])
		end
	end
end


# evolution function: ( population, evolutionary dynamics data )
function evolution!(pop::tLivingPop, evo::tEvoData)#::Int8
	R = let m = MersenneTwister(1)
	        [m; accumulate(Future.randjump, fill(big(10)^20, nthreads()-1), init=m)]
	    end;

	for gen in 1:evo.Ngen
		replication!(pop,R)
		effMutation!(pop,R)
		evo.growthFactor[gen] = effSelection!(pop,ubermode=true)
		evo.aveFitness[gen] = sum([pop.aGty[i].pF[1] for i in 1:pop.pN[2]])/pop.pN[2]
	end
end

export replication!, effMutation!, effSelection!, evolution!


# *******************
# ISING
# *******************

# ===================
# MonteCarlo simulation function for ising signal transduction: readout magnetization evaluation
# 	( tIsingST isingST, interaction matrix Jij, input field h ) → readout magnetization m
function metropolis(isingST::tIsingSigTransEty, Jij::Array{T,1}, hi::Real)::Real where {T<:Real}
	Nmcs::Int32 =  50isingST.L2		# number of Monte Carlo steps
	Nsamplings::Int16 = 20				# number of samplings points
	𝕁::Array{Float64,1} = Jij*isingST.β;	ℍi::Float64 = hi*isingST.β

	# the initial state vector
	n = (isingST.he != 0.0 ? sign(isingST.he) : 1)*ones(Int8, isingST.L, isingST.L)

	# definition: time-averaged readout magnetization
	mro::Float64 = 0.0

	for is in 1:Nsamplings
		for imcs in 1:Nmcs
			# Monte Carlo step evaluated using Glauber transition rates:
			# next possibly transitioning spin (coordinates)
			i, j = rand(collect(1:isingST.L)), rand(collect(1:isingST.L))
			# transitioning?!
			if i <= isingST.li && j <= isingST.li
				if rand() < ( 1. - n[i,j]*tanh(
						𝕁[isingST.Jpi[i,j]]*n[isingST.jp[i],j] + 𝕁[isingST.Jmi[i,j]]*n[isingST.jm[i],j] +
						𝕁[isingST.Jpj[i,j]]*n[i,isingST.jp[j]] + 𝕁[isingST.Jmj[i,j]]*n[i,isingST.jm[j]] +
						isingST.ℍe + ℍi ))/2
					n[i,j] = - n[i,j]
				end
			else
				if rand() < ( 1. - n[i,j]*tanh(
						𝕁[isingST.Jpi[i,j]]*n[isingST.jp[i],j] + 𝕁[isingST.Jmi[i,j]]*n[isingST.jm[i],j] +
						𝕁[isingST.Jpj[i,j]]*n[i,isingST.jp[j]] + 𝕁[isingST.Jmj[i,j]]*n[i,isingST.jm[j]] +
						isingST.ℍe ))/2
					n[i,j] = - n[i,j]
				end
			end
		end

		# evaluation: time-averaged readout magnetization
		mro += sum(n[isingST.halfL+1:isingST.halfL+isingST.li,isingST.halfL+1:isingST.halfL+isingST.li])
	end

	# evaluation: time-averaged readout magnetization
	mro /= Nsamplings * isingST.li2

	return mro
end

# ===================
# MonteCarlo simulation function for ising signal transduction: readout magnetization evaluation
# 	( tIsingST isingST, state n, interaction matrix Jij, input field h ) → readout magnetization m
function metropolis(isingST::tIsingSigTransEty, n::Array{T1,2}, Jij::Array{T2,1}, hi::Real)::Real where
		{T1<:Integer,T2<:Real}
	Nmcs::Int32 =  50isingST.L2		# number of Monte Carlo steps
	Nsamplings::Int16 = 20				# number of samplings points
	𝕁::Array{Float64,1} = Jij*isingST.β;	ℍi::Float64 = hi*isingST.β

	# definition: time-averaged readout magnetization
	mro::Float64 = 0.0

	for is in 1:Nsamplings
		for imcs in 1:Nmcs
			# Monte Carlo step evaluated using Glauber transition rates:
			# next possibly transitioning spin (coordinates)
			i, j = rand(collect(1:isingST.L)), rand(collect(1:isingST.L))
			# transitioning?!
			if i <= isingST.li && j <= isingST.li
				if rand() < ( 1. - n[i,j]*tanh(
						𝕁[isingST.Jpi[i,j]]*n[isingST.jp[i],j] + 𝕁[isingST.Jmi[i,j]]*n[isingST.jm[i],j] +
						𝕁[isingST.Jpj[i,j]]*n[i,isingST.jp[j]] + 𝕁[isingST.Jmj[i,j]]*n[i,isingST.jm[j]] +
						isingST.ℍe + ℍi ))/2
					n[i,j] = - n[i,j]
				end
			else
				if rand() < ( 1. - n[i,j]*tanh(
						𝕁[isingST.Jpi[i,j]]*n[isingST.jp[i],j] + 𝕁[isingST.Jmi[i,j]]*n[isingST.jm[i],j] +
						𝕁[isingST.Jpj[i,j]]*n[i,isingST.jp[j]] + 𝕁[isingST.Jmj[i,j]]*n[i,isingST.jm[j]] +
						isingST.ℍe ))/2
					n[i,j] = - n[i,j]
				end
			end
		end

		# evaluation: time-averaged readout magnetization
		mro += sum(n[isingST.halfL+1:isingST.halfL+isingST.li,isingST.halfL+1:isingST.halfL+isingST.li])
	end

	# evaluation: time-averaged readout magnetization
	mro /= Nsamplings * isingST.li2

	return mro
end

# ===================
# MonteCarlo simulation function for ising signal transduction: time-averaged spin config evaluation
# 	( tIsingST isingST, state n, interaction matrix Jij, input field h ) → readout magnetization m
function metropolis(isingST::tIsingSigTransEty, n::Array{T1,2}, Jij::Array{T2,1}, hi::Real, aves::Array{T3,2}) where
		{T1<:Integer,T2<:Real,T3<:Real}
	Nmcs::Int32 =  50isingST.L2		# number of Monte Carlo steps
	Nsamplings::Int16 = 20				# number of samplings points
	𝕁::Array{Float64,1} = Jij*isingST.β;	ℍi::Float64 = hi*isingST.β

	# initialization: time-averaged spin config
	aves = zeros(eltype(aves), isingST.L, isingST.L)

	for is in 1:Nsamplings
		for imcs in 1:Nmcs
			# Monte Carlo step evaluated using Glauber transition rates:
			# next possibly transitioning spin (coordinates)
			i, j = rand(collect(1:isingST.L)), rand(collect(1:isingST.L))
			# transitioning?!
			if i <= isingST.li && j <= isingST.li
				if rand() < ( 1. - n[i,j]*tanh(
						𝕁[isingST.Jpi[i,j]]*n[isingST.jp[i],j] + 𝕁[isingST.Jmi[i,j]]*n[isingST.jm[i],j] +
						𝕁[isingST.Jpj[i,j]]*n[i,isingST.jp[j]] + 𝕁[isingST.Jmj[i,j]]*n[i,isingST.jm[j]] +
						isingST.ℍe + ℍi ))/2
					n[i,j] = - n[i,j]
				end
			else
				if rand() < ( 1. - n[i,j]*tanh(
						𝕁[isingST.Jpi[i,j]]*n[isingST.jp[i],j] + 𝕁[isingST.Jmi[i,j]]*n[isingST.jm[i],j] +
						𝕁[isingST.Jpj[i,j]]*n[i,isingST.jp[j]] + 𝕁[isingST.Jmj[i,j]]*n[i,isingST.jm[j]] +
						isingST.ℍe ))/2
					n[i,j] = - n[i,j]
				end
			end
		end

		# evaluation: time-averaged spin config
		aves += n
	end

	# evaluation: time-averaged spin config
	aves /= Nsamplings

	matshow(aves,cmap="Greys_r"); gcf()
end

# fitness function for ising signal transduction
# 	( evotype isingST, genotype gty, environment isingSTenv )
function fitness(isingST::tIsingSigTransEty, gty::at1dGty{<:Real}, isingSTenv::tCompEnv{<:Real})::Float64
	# where {Tgty,Tenv<:Real}
	invf::Float64 = 0.0
	for iio in isingSTenv.idealInputOutput
		invf += ( metropolis(isingST, gty.X, iio[1]) - iio[2] )^2
	end
	return 1/(sqrt(invf)+0.1)	# fitness offset = 0.1
end

function fitness!(isingST::tIsingSigTransEty, isingSTenv::tCompEnv{<:Real}, gty::at1dGty{<:Real})
	# where {Tgty,Tenv<:Real}
	invf::Float64 = 0.0
	for iio in isingSTenv.idealInputOutput
		invf += ( metropolis(isingST, gty.X, iio[1]) - iio[2] )^2
	end
	gty.pF[1] = 1/(sqrt(invf)+0.1)	# fitness offset = 0.1
end

# constructor of tLivingPop based on ising fitness
# function tLivingPop{T}( N::Int32,ety::tIsingSigTransEty,env::atEnv,aGty::Array{<:atGenotype{T},1},
# 	repFactor::Float64,mutFactor::Float64,Xvar::T,ΔtOffset::Float64 ) where {T}
# 	for i in 1:N fitness!(ety,env,aGty[i]) end
# 	tLivingPop{T}( Int32[N,N,length(aGty)],ety,env,aGty )
# end

# simplified constructor for discrete time evolution
function tLivingPop{T}( N::Int32,ety::atEvotype,env::atEnv,aGty::Array{<:atGenotype{T},1},
		repFactor::Float64,mutFactor::Float64,Xvar::T,ΔtOffset::Float64 ) where {T}
	for i in 1:N fitness!(ety,env,aGty[i]) end

	tLivingPop{T}( Int32[N,N,length(aGty)],ety,env,aGty,
	repFactor/(2maximum([aGty[i].pdX[1] for i in 1:N])*mutFactor+ΔtOffset),
	mutFactor/(2maximum([aGty[i].pdX[1] for i in 1:N])*mutFactor+ΔtOffset),
	Xvar )
end

export metropolis, fitness, fitness!, tLivingPop


# *******************
# TRIVIAL
# *******************

struct tTrivialEty <: atEvotype end
struct tTrivialEnv <: atEnv end

function fitness!(trivialEty::tTrivialEty,trivialEnv::tTrivialEnv,gty::at1dGty{<:Real})
	gty.pF[1]=1/(euclidean(gty.X,ones(Float64,gty.pdX[1]))+1.)
end

export tTrivialEty, tTrivialEnv, fitness!

end
