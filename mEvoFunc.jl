
module mEvoFunc
using mEvoTypes, Random, Base.Threads, DelimitedFiles, Dates, Distances, ProgressMeter
import Future

const FITNESSOFFSET, FITNESSTHRESHOLD = .001, .95

# initializer constructor for discrete time evolution. Same MetaGenotype for all.
function initLivingPop( N::Int32,ety::Tety,env::Tenv,aMGty::Vector{Tmgty},aGty::Vector{Tgty} ) where {
		Tety<:atEvotype,Tenv<:atEnvironment,Tmgty<:atMetaGenotype,Tgty<:atGenotype }
	for i in 1:N fitness!(env,aGty[i]) end
	return tLivingPop{Tety,Tenv,Vector{Tmgty},Vector{Tgty}}( Int32[N,N,length(aGty)],ety,env,aMGty,aGty )
end

tEty{Tx}(repFactor::Float64,mutFactor::Float64,ΔtOffset::Float64,dX::Int32,Xvar::Tx) where {Tx<:Number} =
	tEty{Tx}( [repFactor/(2dX*mutFactor+ΔtOffset)], [mutFactor/(2dX*mutFactor+ΔtOffset)], Xvar )

# function changing rep and mut -factors
function set_tEtyFactors(ety::tEty,repFactor::Float64,mutFactor::Float64,ΔtOffset::Float64,dX::Int32)
	ety.pRepFactor[1] = repFactor/(2dX*mutFactor+ΔtOffset)
	ety.pMutFactor[1] = mutFactor/(2dX*mutFactor+ΔtOffset)
end

# function: saving the genotypes for VecGty's
function write_aGty(pop::tLivingPop{<:atEvotype,<:atEnvironment,<:Vector{<:atMetaGenotype},<:Vector{<:atVecGty}})
	open( "population_" * string(now()) * ".dat", "w" ) do f
		for i in 1:pop.pN[2]
			print(f,pop.aGty[i].pMetaGty[1].dX,"\t")
			for x in pop.aGty[i].X
				print(f,x,"\t")
			end
			print(f,pop.aGty[i].pF[1],"\n")
		end
	end
end

# read population and spits metagenotype and genotype
function read_aIsingSigTransGty( aGtyFileName::String,β::Float64,hi::Float64 )
	matGty = readdlm(aGtyFileName)
	Npop = size(matGty)[1]
	sysSizes = Int32[]
	L = zeros(Int32,Npop)

	for i in 1:Npop
		L[i] = Int32(sqrt(matGty[i,1]/2))
		if !( L[i] in sysSizes )
			push!(sysSizes,L[i])
		end
	end

	aMGty = [ tIsingSigTransMGty(L,β,hi) for L in sysSizes ]
	return aMGty, [ tVecGty(
		[ aMGty[collect(1:length(sysSizes))[sysSizes .== L[i]][1]] ], matGty[i,2:1+Int32(matGty[i,1])] , Float64[matGty[i,end]] ) for i in 1:Npop ]
end

export initLivingPop, write_aGty, read_aIsingSigTransGty

# *********************************
# | PROPER EVOLUTIONARY FUNCTIONS \
# *********************************

# replication function: ( living population, replication factor, Marsenne Twister ) → replicated population
function replication!(pop::tLivingPop,R::Vector{MersenneTwister})
	G = zeros(Int32,pop.pN[2])
	# @threads
	for i in 1:pop.pN[2]
		Kr::Float64 = pop.ety.pRepFactor[1]*pop.aGty[i].pF[1]
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
function effMutation!(pop::tLivingPop,R::Vector{MersenneTwister})
	Nmutations::Int32 = 0
	# Nmutations = Atomic{Int32}(0) # @threads
	for i in 1:pop.pN[1]
		r1::Float64 = rand(R[threadid()])
		cumProb::Float64 = 0.
		xvar::Int32 = -1

		effMutProb::Float64 = pop.ety.pMutFactor[1]/(1+pop.ety.pRepFactor[1]*pop.aGty[i].pF[1])

		while cumProb < r1 && xvar < 2pop.aGty[i].pMetaGty[1].dX
			xvar += 1
			cumProb += effMutProb
		end
		if xvar < 2pop.aGty[i].pMetaGty[1].dX
			pop.aGty[i].X[xvar%pop.aGty[i].pMetaGty[1].dX+1] += xvar < pop.aGty[i].pMetaGty[1].dX ? pop.ety.Xvar : -pop.ety.Xvar
			fitness!(pop.env,pop.aGty[i])
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

function upgradeGtyX!(gty::tVecGty{Vector{TMGty},Vector{Tx}}) where {TMGty<:atMetaGenotype,Tx}
	JijRef = copy(gty.X)
	append!( gty.X, ones(Float64, 8gty.pMetaGty[1].L+8) )
	# append!( gty.X, rand(-4:.5:5, 4gty.pMetaGty[1].L+8) )
	gty.X .= 1.

	for i in 1:gty.pMetaGty[1].halfL, j in 1:gty.pMetaGty[1].L,
		gty.X[i+(j-1)*(gty.pMetaGty[1].L+2)] = JijRef[i+(j-1)*gty.pMetaGty[1].L]
	end
	for i in gty.pMetaGty[1].halfL+2:gty.pMetaGty[1].L+1, j in 1:gty.pMetaGty[1].L,
		gty.X[i+(j-1)*(gty.pMetaGty[1].L+2)] = JijRef[i-1+(j-1)*gty.pMetaGty[1].L]
	end
	for i in 1:gty.pMetaGty[1].halfL, j in gty.pMetaGty[1].L+3:2gty.pMetaGty[1].L+2,
		gty.X[i+(j-1)*(gty.pMetaGty[1].L+2)] = JijRef[i+(j-3)*gty.pMetaGty[1].L]
	end
	for i in gty.pMetaGty[1].halfL+2:gty.pMetaGty[1].L+1, j in gty.pMetaGty[1].L+3:2gty.pMetaGty[1].L+2,
		gty.X[i+(j-1)*(gty.pMetaGty[1].L+2)] = JijRef[i-1+(j-3)*gty.pMetaGty[1].L]
	end
end

# fitness!(istEnv::tCompEnv{<:Array{Float64}},gty::tVecGty)
function evoUpgrade!(pop::tLivingPop{<:atEvotype,<:atEnvironment,Vector{tIsingSigTransMGty},<:Vector{<:atVecGty}})
	for i in 1:pop.pN[2]
		if pop.aGty[i].pF[1] > FITNESSTHRESHOLD - (pop.aGty[i].pMetaGty[1].L - 6.)/2.
			upgradeGtyX!(pop.aGty[i])
			fitness!(pop.env,pop.aGty[i])

			foundMetaGty = false
			for metaGty in pop.aMetaGty
				if pop.aGty[i].pMetaGty[1].L+2 == metaGty.L
					foundMetaGty = true
					pop.aGty[i].pMetaGty[1] = metaGty
				end
			end
			if !foundMetaGty
				push!( pop.aMetaGty, tIsingSigTransMGty(pop.aGty[i].pMetaGty[1].L+Int32(2),pop.aGty[i].pMetaGty[1].β,pop.aGty[i].pMetaGty[1].he) )
				pop.aGty[i].pMetaGty[1] = pop.aMetaGty[end]
			end
		end
	end
end

# evolution function: ( population, evolutionary dynamics data )
function evolution!(pop::tLivingPop,evo::tEvoData; ubermode::Bool=false)#::Int8
	R = let m = MersenneTwister(1)
	        [m; accumulate(Future.randjump, fill(big(10)^20, nthreads()-1), init=m)]
	    end;

	@showprogress 1 "Evolutionary Dynamics Status: " for gen in 1:evo.Ngen
		evo.growthFactor[gen] = replication!(pop,R)
		evo.mutationFactor[gen] = effMutation!(pop,R)
		effSelection!(pop,ubermode)
		evo.aveFitness[gen] = sum([pop.aGty[i].pF[1] for i in 1:pop.pN[2]])/pop.pN[2]
	end
end

export replication!, effMutation!, effSelection!, evoUpgrade!, evolution!


# ***********
# *  ISING  *
# ***********

# function: flipping --- for metropolis
function flipping!(istMGty::tIsingSigTransMGty, βJij::Array{<:Real,1}, βhi::Real, n::Array{<:Integer,2})
	# Monte Carlo step evaluated using Glauber transition rates:
	# next possibly transitioning spin (coordinates)
	i, j = rand(collect(1:istMGty.L)), rand(collect(1:istMGty.L))
	# transitioning?!
	if i <= istMGty.li && j <= istMGty.li
		if rand() < ( 1. - n[i,j]*tanh(
				βJij[istMGty.Jpi[i,j]]*n[istMGty.jp[i],j] + βJij[istMGty.Jmi[i,j]]*n[istMGty.jm[i],j] +
				βJij[istMGty.Jpj[i,j]]*n[i,istMGty.jp[j]] + βJij[istMGty.Jmj[i,j]]*n[i,istMGty.jm[j]] +
				istMGty.βhe + βhi ))/2
			n[i,j] = - n[i,j]
		end
	else
		if rand() < ( 1. - n[i,j]*tanh(
				βJij[istMGty.Jpi[i,j]]*n[istMGty.jp[i],j] + βJij[istMGty.Jmi[i,j]]*n[istMGty.jm[i],j] +
				βJij[istMGty.Jpj[i,j]]*n[i,istMGty.jp[j]] + βJij[istMGty.Jmj[i,j]]*n[i,istMGty.jm[j]] +
				istMGty.βhe ))/2
			n[i,j] = - n[i,j]
		end
	end
end

# ===================
# MonteCarlo simulation function for ising signal transduction: readout magnetization evaluation
# 	( tIsingST istMGty, interaction matrix Jij, input field h ) → readout magnetization m
function metropolis(istMGty::tIsingSigTransMGty,Jij::Array{T,1},hi::T)::Real where {T<:Real}
	βJij::Array{Float64,1} = Jij*istMGty.β;		βhi::Float64 = hi*istMGty.β

	# the initial state vector
	n = Int8(istMGty.he != 0.0 ? sign(istMGty.he) : 1).*ones(Int8, istMGty.L, istMGty.L)

	# definition: time-averaged readout magnetization
	mro::Float64 = 0.0

	for is in 1:istMGty.prms.Nsmpl
		for imcs in 1:istMGty.prms.Nmcsps
			flipping!(istMGty,βJij,βhi,n)
		end
		# evaluation: time-averaged readout magnetization
		mro += sum(n[istMGty.halfL+1:istMGty.halfL+istMGty.li,istMGty.halfL+1:istMGty.halfL+istMGty.li])
	end
	# evaluation: time-averaged readout magnetization
	mro /= istMGty.prms.Nsmpl*istMGty.li2

	return mro
end

# ===================
# MonteCarlo simulation function for ising signal transduction: readout magnetization evaluation, n evolution
# 	( tIsingST istMGty, state n, interaction matrix Jij, input field h ) → readout magnetization m
function metropolis!(istMGty::tIsingSigTransMGty,Jij::Array{<:Real,1},hi::Real,n::Array{<:Integer,2})
	βJij::Array{Float64,1} = Jij*istMGty.β;		βhi::Float64 = hi*istMGty.β

	# definition: time-averaged readout magnetization and readout region range
	mro::Float64 = 0.0
	roRegionRange = istMGty.halfL+1:istMGty.halfL+istMGty.li

	for is in 1:istMGty.prms.Nsmpl
		for imcs in 1:istMGty.prms.Nmcsps
			flipping!(istMGty,βJij,βhi,n)
		end
		# evaluation: time-averaged readout magnetization
		mro += sum(view(n,roRegionRange,roRegionRange))
	end
	# evaluation: time-averaged readout magnetization
	mro /= istMGty.prms.Nsmpl * istMGty.li2

	return mro
end

# ===================
# MonteCarlo simulation function for ising signal transduction: time-averaged spin config evaluation
# 	( tIsingST istMGty, state n, interaction matrix Jij, input field h ) → readout magnetization m
function metropolis!(istMGty::tIsingSigTransMGty,Jij::Array{<:Real,1},hi::Real,aves::Array{Float64,2})
	βJij::Array{Float64,1} = Jij.*istMGty.β;		βhi::Float64 = hi.*istMGty.β

	# initialization: initial state vector
	n = Int8(istMGty.he != 0.0 ? sign(istMGty.he) : 1).*ones(Int8, istMGty.L, istMGty.L)

	for is in 1:istMGty.prms.Nsmpl
		for imcs in 1:istMGty.prms.Nmcsps
			flipping!(istMGty,βJij,βhi,n)
		end
		# evaluation: time-averaged spin config
		aves .+= n
	end
	# evaluation: time-averaged spin config
	aves ./= istMGty.prms.Nsmpl
end

# fitness function for ising signal transduction
# 	( evotype istMGty, genotype gty, environment istEnv )
function fitness!(istEnv::tCompEnv{<:Array{Float64}},gty::tVecGty)
	d2::Float64 = 0.0
	for iio in istEnv.idealInputOutput
		d2 += ( metropolis(gty.pMetaGty[1],broadcast(i->exp(i),gty.X),iio[1]) - iio[2] )^2
	end
	gty.pF[1] = exp(-sqrt(d2)/istEnv.selFactor) + (gty.pMetaGty[1].L - 6.)/2.
end

function fitness(istEnv::tCompEnv{<:Array{Float64}},gty::tVecGty)::Float64
	d2::Float64 = 0.0
	for iio in istEnv.idealInputOutput
		d2 += ( metropolis(gty.pMetaGty[1],broadcast(i->exp(i),gty.X),iio[1]) - iio[2] )^2
	end
	return exp(-sqrt(d2)/istEnv.selFactor) + (gty.pMetaGty[1].L - 6.)/2.
end

# function: showing the spin config of the ising signal transduction system
function showPhenotype!(env::tCompEnv,gty::tVecGty,aAves::Array{Array{Float64,2},1})
	i::Int32 = 0
	for iio in env.idealInputOutput
		metropolis!( gty.pMetaGty[1], broadcast(i->exp(i),gty.X), iio[1], aAves[i+=1] )
	end
end

function showGenotype!(gty::tVecGty,JijMat::Array{Float64,2})
	ii::Int32 = 0
	for j in 1:2gty.pMetaGty[1].L, i in 1:2gty.pMetaGty[1].L
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
	gty.pF[1]=1/(euclidean(gty.X,ones(Float64,gty.pMetaGty[1].dX))+FITNESSOFFSET)
end

export tTrivialEty, tTrivialEnv, fitness!

end
