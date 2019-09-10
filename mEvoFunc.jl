
module mEvoFunc
using mEvoTypes, Random, Base.Threads, DelimitedFiles, Dates, Distances, ProgressMeter
import Future

const FITNESSOFFSET = .001

# simplified constructor for discrete time evolution
function tLivingPop{T, Tevo, Tenv, TaGty}( N::Int32,ety::Tevo,env::Tenv,aGty::Array{<:atGenotype,1},
		repFactor::Float64,mutFactor::Float64,Xvar::T,ΔtOffset::Float64,prm::tDTMCprm ) where {T,
		Tevo<:atEvotype, Tenv<:atEnvironment, TaGty<:Array{<:atGenotype,1}}
	for i in 1:N fitness!(ety,env,aGty[i],prm) end

	tLivingPop{T, Tevo, Tenv, TaGty}( Int32[N,N,length(aGty)],ety,env,aGty,
	repFactor/(2maximum([aGty[i].pdX[1] for i in 1:N])*mutFactor+ΔtOffset),
	mutFactor/(2maximum([aGty[i].pdX[1] for i in 1:N])*mutFactor+ΔtOffset),
	Xvar )
end

# function: saving the genotypes for 1dGty's
function write_aGty(pop::tLivingPop{<:Number,<:atEvotype,<:atEnvironment,<:Array{<:at1dGty,1}})
	open( "population_" * string(now()) * ".dat", "w" ) do f
		for i in 1:pop.pN[2]
			print(f,pop.aGty[i].pdX[1],"\t")
			for x in pop.aGty[i].X
				print(f,x,"\t")
			end
			print(f,pop.aGty[i].pF[1],"\n")
		end
	end
end

# simplified constructor for stored populations
function tLivingPop{T, Tevo, Tenv, TaGty}( ety::Tevo,env::Tenv,aGtyFileName::String,
		repFactor::Float64,mutFactor::Float64,Xvar::T,ΔtOffset::Float64 ) where {T,
		Tevo<:atEvotype, Tenv<:atEnvironment, TaGty<:Array{<:atGenotype,1}}
	let mGty = readdlm(aGtyFileName), N = size(mGty)[1]
		tLivingPop{T,Tevo,Tenv,TaGty}( Int32[N,N,N],ety,env,
		[ tVecGty{Array{T,1}}([Int32(mGty[i,1])],view(mGty,i,2:1+Int32(mGty[i,1])),Float64[mGty[i,end]]) for i in 1:N ],
		repFactor/(2maximum([mGty[i,1] for i in 1:N])*mutFactor+ΔtOffset),
		mutFactor/(2maximum([mGty[i,1] for i in 1:N])*mutFactor+ΔtOffset),
		Xvar )
	end
end

export tLivingPop, write_aGty

# *********************************
# | PROPER EVOLUTIONARY FUNCTIONS \
# *********************************

# replication function: ( living population, replication factor, Marsenne Twister ) → replicated population
function replication!(pop::tLivingPop,R::Vector{MersenneTwister})
	G = zeros(Int32,pop.pN[2])
	# @threads
	for i in 1:pop.pN[2]
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
	# returning growth factor
	return log(pop.pN[1]/pop.pN[2])
end

# effective mutation function: ( population, variation, fitness!Function acting on genotypes, MT )
function effMutation!(pop::tLivingPop,prm::tDTMCprm,R::Vector{MersenneTwister})
	Nmutations::Int32 = 0
	# Nmutations = Atomic{Int32}(0) # @threads
	for i in 1:pop.pN[1]
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
			fitness!(pop.ety,pop.env,pop.aGty[i],prm)
			Nmutations += Int32(1)
		end
	end
	# returning the total number of mutations normalized by population
	return Nmutations/pop.pN[1]
end

# effective selection function: ( population, MT ) → selected population
function effSelection!(pop::tLivingPop, ubermode::Bool)
	popGtyRef::Array{atGenotype,1} = copy(pop.aGty)

	if ubermode
		# survival of the fittest
		survivedGty = sortperm([pop.aGty[i].pF[1] for i in 1:pop.pN[1]],rev=true)
		# sort!(pop.aGty, by= x -> x.pF[1], rev=true)
	else
		# selection with replacement of individuals
		survivedGty = rand(1:pop.pN[1],pop.pN[2])
	end
	for i in 1:pop.pN[2]
		pop.aGty[i] = popGtyRef[survivedGty[i]]
	end

	# population size renormalization
	pop.pN[1] = pop.pN[2]
end

# evolution function: ( population, evolutionary dynamics data )
function evolution!(pop::tLivingPop,evo::tEvoData,prm::tDTMCprm; ubermode::Bool=false)#::Int8
	R = let m = MersenneTwister(1)
	        [m; accumulate(Future.randjump, fill(big(10)^20, nthreads()-1), init=m)]
	    end;

	@showprogress 1 "Evolutionary Dynamics Status: " for gen in 1:evo.Ngen
		evo.growthFactor[gen] = replication!(pop,R)
		evo.mutationNumber[gen] = effMutation!(pop,prm,R)
		effSelection!(pop,ubermode)
		evo.aveFitness[gen] = sum([pop.aGty[i].pF[1] for i in 1:pop.pN[2]])/pop.pN[2]
	end
end

export replication!, effMutation!, effSelection!, evolution!


# ***********
# *  ISING  *
# ***********

# function: flipping --- for metropolis
function flipping!(isingST::tIsingSigTransEty, βJij::Array{<:Real,1}, βhi::Real, n::Array{<:Integer,2})
	# Monte Carlo step evaluated using Glauber transition rates:
	# next possibly transitioning spin (coordinates)
	i, j = rand(collect(1:isingST.L)), rand(collect(1:isingST.L))
	# transitioning?!
	if i <= isingST.li && j <= isingST.li
		if rand() < ( 1. - n[i,j]*tanh(
				βJij[isingST.Jpi[i,j]]*n[isingST.jp[i],j] + βJij[isingST.Jmi[i,j]]*n[isingST.jm[i],j] +
				βJij[isingST.Jpj[i,j]]*n[i,isingST.jp[j]] + βJij[isingST.Jmj[i,j]]*n[i,isingST.jm[j]] +
				isingST.ℍe + βhi ))/2
			n[i,j] = - n[i,j]
		end
	else
		if rand() < ( 1. - n[i,j]*tanh(
				βJij[isingST.Jpi[i,j]]*n[isingST.jp[i],j] + βJij[isingST.Jmi[i,j]]*n[isingST.jm[i],j] +
				βJij[isingST.Jpj[i,j]]*n[i,isingST.jp[j]] + βJij[isingST.Jmj[i,j]]*n[i,isingST.jm[j]] +
				isingST.ℍe ))/2
			n[i,j] = - n[i,j]
		end
	end
end

# ===================
# MonteCarlo simulation function for ising signal transduction: readout magnetization evaluation
# 	( tIsingST isingST, interaction matrix Jij, input field h ) → readout magnetization m
function metropolis(isingST::tIsingSigTransEty,Jij::Array{T,1},hi::Real,prm::tDTMCprm)::Real where {T<:Real}
	βJij::Array{Float64,1} = Jij*isingST.β;		βhi::Float64 = hi*isingST.β

	# the initial state vector
	n = Int8(isingST.he != 0.0 ? sign(isingST.he) : 1).*ones(Int8, isingST.L, isingST.L)

	# definition: time-averaged readout magnetization
	mro::Float64 = 0.0

	for is in 1:prm.Nsmpl
		for imcs in 1:prm.Nmcsps
			flipping!(isingST,βJij,βhi,n)
		end
		# evaluation: time-averaged readout magnetization
		mro += sum(n[isingST.halfL+1:isingST.halfL+isingST.li,isingST.halfL+1:isingST.halfL+isingST.li])
	end
	# evaluation: time-averaged readout magnetization
	mro /= prm.Nsmpl * isingST.li2

	return mro
end

# ===================
# MonteCarlo simulation function for ising signal transduction: readout magnetization evaluation, n evolution
# 	( tIsingST isingST, state n, interaction matrix Jij, input field h ) → readout magnetization m
function metropolis!(isingST::tIsingSigTransEty,Jij::Array{<:Real,1},hi::Real,n::Array{<:Integer,2},prm::tDTMCprm)
	βJij::Array{Float64,1} = Jij*isingST.β;		βhi::Float64 = hi*isingST.β

	# definition: time-averaged readout magnetization and readout region range
	mro::Float64 = 0.0
	roRegionRange = isingST.halfL+1:isingST.halfL+isingST.li

	for is in 1:prm.Nsmpl
		for imcs in 1:prm.Nmcsps
			flipping!(isingST,βJij,βhi,n)
		end
		# evaluation: time-averaged readout magnetization
		mro += sum(view(n,roRegionRange,roRegionRange))
	end
	# evaluation: time-averaged readout magnetization
	mro /= prm.Nsmpl * isingST.li2

	return mro
end

# ===================
# MonteCarlo simulation function for ising signal transduction: time-averaged spin config evaluation
# 	( tIsingST isingST, state n, interaction matrix Jij, input field h ) → readout magnetization m
function metropolis!(isingST::tIsingSigTransEty,Jij::Array{<:Real,1},hi::Real,aves::Array{Float64,2},prm::tDTMCprm)
	βJij::Array{Float64,1} = Jij.*isingST.β;		βhi::Float64 = hi.*isingST.β

	# initialization: initial state vector
	n = Int8(isingST.he != 0.0 ? sign(isingST.he) : 1).*ones(Int8, isingST.L, isingST.L)

	for is in 1:prm.Nsmpl
		for imcs in 1:prm.Nmcsps
			flipping!(isingST,βJij,βhi,n)
		end
		# evaluation: time-averaged spin config
		aves .+= n
	end
	# evaluation: time-averaged spin config
	aves ./= prm.Nsmpl
end

# fitness function for ising signal transduction
# 	( evotype isingST, genotype gty, environment isingSTenv )
function fitness!(isingST::tIsingSigTransEty,isingSTenv::tCompEnv{<:Array{Float64}},gty::tVecGty,prm::tDTMCprm)
	d2::Float64 = 0.0
	for iio in isingSTenv.idealInputOutput
		d2 += (metropolis( isingST,broadcast(i->exp(i),gty.X),iio[1],prm ) - iio[2])^2
	end
	gty.pF[1] = exp(-sqrt(d2)/isingSTenv.selFactor)
end

function fitness(isingST::tIsingSigTransEty,isingSTenv::tCompEnv{<:Array{Float64}},gty::tVecGty,
		prm::tDTMCprm)::Float64
	d2::Float64 = 0.0
	for iio in isingSTenv.idealInputOutput
		d2 += (metropolis( isingST,broadcast(i->exp(i),gty.X),iio[1],prm ) - iio[2])^2
	end
	return exp(-sqrt(d2)/isingSTenv.selFactor)
end

# function: showing the spin config of the ising signal transduction system
# 	( evotype isingST, environment isingSTenv, genotype gty )
function showPhenotype!(isingST::tIsingSigTransEty,isingSTenv::tCompEnv{<:Array{Float64}},gty::tVecGty,
		aAves::Array{Array{Float64,2},1},prm::tDTMCprm)
	i::Int32 = 1
	for iio in isingSTenv.idealInputOutput
		metropolis!(isingST,broadcast(i->exp(i),gty.X),iio[1],aAves[i],prm)
		i+=1
	end
end

function showGenotype!(isingST::tIsingSigTransEty,gty::tVecGty,JijMat::Array{Float64,2})
	ii::Int32 = 0
	for j in 1:2isingST.L, i in 1:2isingST.L
		JijMat[i,j] = j%2==1 ? ( i%2==0 ? exp(gty.X[ii+=1]) : -1 ) : ( i%2==1 ? exp(gty.X[ii+=1]) : 10^6+1 )
	end
end

export metropolis, metropolis!, fitness, fitness!, showPhenotype!, showGenotype!


# *******************
# TRIVIAL
# *******************

struct tTrivialEty <: atEvotype end
struct tTrivialEnv <: atEnvironment end

function fitness!(trivialEty::tTrivialEty,trivialEnv::tTrivialEnv,gty::tVecGty{Array{T,1}}) where {T<:Real}
	gty.pF[1]=1/(euclidean(gty.X,ones(Float64,gty.pdX[1]))+FITNESSOFFSET)
end

export tTrivialEty, tTrivialEnv, fitness!

end
