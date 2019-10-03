
module mEvoFunc
using mEvoTypes, Random, Statistics, Base.Threads, DelimitedFiles, Dates, Distances, ProgressMeter
import Future

# threads random generators initialization
const THREADRNG = let m = MersenneTwister(1)
            [m; accumulate(Future.randjump, fill(big(10)^20, nthreads()-1), init=m)]
        end;

const FITNESSOFFSET, FITNESSTHRESHOLD, BASALFITNESS = .001, .95, 3.0

# initializer constructor for discrete time evolution. Same MetaGenotype for all.
function initLivingPop( N::Int32,ety::Tety,env::Tenv,aMGty::Vector{Tmgty},aGty::Vector{Tgty} ) where {
		Tety<:atEvotype,Tenv<:atEnvironment,Tmgty<:atMetaGenotype,Tgty<:atGenotype }
	for i in 1:N fitness!(env,aGty[i]) end
	return tLivingPop{Tety,Tenv,Vector{Tmgty},Vector{Tgty}}( Int32[N,N,length(aGty)],ety,env,aMGty,aGty )
end

# function. changing rep and mut -factors in tEty
function set_tEtyFactors(ety::tEty,gty::atGenotype)
	ety.pRepFactor[1] = ety.repRate/(2gty.pdG[1]*ety.mutRate+ety.ΔtOffset)
	ety.pMutFactor[1] = ety.mutRate/(2gty.pdG[1]*ety.mutRate+ety.ΔtOffset)
end

# function. changing rep and mut -factors in tEty
function set_tEtyFactors(ety::tEty,gty::tAlphaGty)
	ety.pRepFactor[1] = ety.repRate/(gty.pdg[1]*gty.pdG[1]*ety.mutRate+ety.ΔtOffset)
	ety.pMutFactor[1] = ety.mutRate/(gty.pdg[1]*gty.pdG[1]*ety.mutRate+ety.ΔtOffset)
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

function blindMutation!(gty::tAddGty,mutProb::Float64)
	Pmut = 2gty.pdG[1]*mutProb
	Pmut <= 1 || throw("probability of mutation exceeds 1")

	cumProb::Float64 = 0.0; ig::Int32 = -1; r::Float64 = rand(THREADRNG[threadid()])
	if r <= Pmut
		while cumProb < r
			ig += 1
			cumProb += mutProb
		end
		gty.G[ ig % gty.pdG[1] + 1 ] += ig % 2 == 0 ? gty.Δg : -gty.Δg
		return Int32(1)
	else
		return Int32(0)
	end
end

function blindMutation!(gty::tMltGty,mutProb::Float64)
	Pmut = 2gty.pdG[1]*mutProb
	Pmut <= 1 || throw("probability of mutation exceeds 1")

	cumProb::Float64 = 0.0; ig::Int32 = -1; r::Float64 = rand(THREADRNG[threadid()])
	if r <= Pmut
		while cumProb < r
			ig += 1
			cumProb += mutProb
		end
		gty.G[ ig % gty.pdG[1] + 1 ] *= ig % 2 == 0 ? gty.δg : 1.0/gty.δg
		return Int32(1)
	else
		return Int32(0)
	end
end

function blindMutation!(gty::tAlphaGty,mutProb::Float64)
	Pmut = gty.pdG[1]*gty.pdg[1]*mutProb
	Pmut <= 1 || throw("probability of mutation exceeds 1")

	# cumProb::Float64 = 0.0; ig::Int32 = -1; r::Float64 = rand(R[threadid()])
	cumProb::Float64 = 0.0; ig::Int32 = -1; r::Float64 = rand(THREADRNG[threadid()])
	if r <= Pmut
		while cumProb < r
			ig += 1
			cumProb += mutProb
		end
		gty.G[ ig % gty.pdG[1] + 1 ] = gty.g[ ig % gty.pdg[1] + 1 ]
		return Int32(1)
	else
		return Int32(0)
	end
end

function blindMutation!(gty::tAlphaGty,mutProb::Float64)
	Pmut = gty.pdG[1]*gty.pdg[1]*mutProb
	Pmut <= 1 || throw("probability of mutation exceeds 1")

	# cumProb::Float64 = 0.0; ig::Int32 = -1; r::Float64 = rand(R[threadid()])
	cumProb::Float64 = 0.0; ig::Int32 = -1; r::Float64 = rand(THREADRNG[threadid()])
	if r <= Pmut
		while cumProb < r
			ig += 1
			cumProb += mutProb
		end
		gty.G[ ig % gty.pdG[1] + 1 ] = gty.g[ ig % gty.pdg[1] + 1 ]
		return Int32(1)
	else
		return Int32(0)
	end
end

# function. effective mutation: blindly mutate the population's genotype according to the effective dynamical mutation rate
function effMutation!(pop::tLivingPop)
	Nmutations = zeros(Int32,nthreads())
	@threads for i in 1:pop.pN[1]
		Nmutations[threadid()] += blindMutation!(pop.aGty[i], pop.ety.pMutFactor[1]/(1+pop.ety.pRepFactor[1]*pop.aGty[i].pF[1]))
		fitness!(pop.env,pop.aGty[i])
	end
	return sum(Nmutations)/pop.pN[1] 	# normalized number of mutations
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
	return gty.pF[1] - (gty.pMetaGty[1].halfL - BASALFITNESS) > FITNESSTHRESHOLD
end

# function. new genotypic variables choice. additive genotype case
function newg!(gty::tAddGty{<:atIsingMetaGty})
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
function newg!(gty::tAlphaGty{<:atIsingMetaGty})
	# duplication of previous line genotype + some randomness
	for k in 0:1, j in 1:gty.pMetaGty[1].L, i in 1:2
		gty.G[i*(gty.pMetaGty[1].halfL+1)+(j+k*(gty.pMetaGty[1].L+2)-1)*(gty.pMetaGty[1].L+2)] = rand(THREADRNG[threadid()],gty.g)
	end
	for u in 0:1, k in 0:1, j in 1:2, i in 1:gty.pMetaGty[1].halfL
		gty.G[i+k*(gty.pMetaGty[1].halfL+1)+(j+gty.pMetaGty[1].L+u*(gty.pMetaGty[1].L+2)-1)*(gty.pMetaGty[1].L+2)] = rand(THREADRNG[threadid()],gty.g)
	end
	# there are a couple of connection missing in this routine
end

# function. upgrade genotypic variables
function upgradeGtyG!(gty::atSystemGty{<:atIsingMetaGty})
	Gref = copy(gty.G)
	append!( gty.G, ones(Float64, 8gty.pMetaGty[1].L+8) )
	gty.pdG[1] += 8gty.pMetaGty[1].L+8

	for u in 0:1, j in 1:gty.pMetaGty[1].L, k in 0:1, i in 1:gty.pMetaGty[1].halfL
		gty.G[ i + k*(gty.pMetaGty[1].halfL + 1) + (j + u*(gty.pMetaGty[1].L + 2) - 1)*(gty.pMetaGty[1].L+2) ] =
			Gref[ i + k*gty.pMetaGty[1].halfL + (j + u*gty.pMetaGty[1].L - 1)*gty.pMetaGty[1].L ]
	end
	newg!(gty)
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
		set_tEtyFactors(ety,gty)
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

# function: genetic evolution
function evolution!(pop::tLivingPop,evo::tEvoData; ubermode::Bool=false)

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

export replication!, effMutation!, effSelection!, evoUpgrade!, evolution!, upgradeGtyG!


# ***********
# *  ISING  *
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

# ===================
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

# ===================
# MonteCarlo simulation function for ising signal transduction: readout magnetization evaluation, n evolution
# 	( tIsingST istMGty, state n, interaction matrix Jij, input field h ) → readout magnetization m
function metropolis!(istMGty::tIsingSigTransMGty,Jij::Array{<:Real,1},hi::Real,n::Array{<:Integer,2})
	βJij::Array{Float64,1} = Jij*istMGty.β;		βhi::Float64 = hi*istMGty.β

	n[1:istMGty.li,1:istMGty.li] .= sign(hi)

	# definition: time-averaged readout magnetization and readout region range
	mro::Float64 = 0.0
	roRegionRange = istMGty.halfL+1:istMGty.halfL+istMGty.li

	for is in 1:istMGty.prms.Nsmpl
		for imcs in 1:istMGty.prms.Nmcsps, ilp in 1:istMGty.L2
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
function metropolis!(istMGty::tIsingSigTransMGty, Jij::Array{<:Real,1}, hi::Real, an::Vector{Array{Int8,2}}, prms::tDTMCprm)
	βJij::Array{Float64,1} = Jij.*istMGty.β;		βhi::Float64 = hi.*istMGty.β

	# initialization: initial state vector
	n = Int8(istMGty.he != 0.0 ? sign(istMGty.he) : 1).*ones(Int8, istMGty.L, istMGty.L)
	n[1:istMGty.li,1:istMGty.li] .= sign(hi)

	@showprogress 1 "Monte Carlo Dynamics Status: " for is in 1:prms.Nsmpl
		for imcs in 1:prms.Nmcsps, ilp in 1:istMGty.L2
			flipping!(istMGty,βJij,βhi,n)
		end
		an[is] = deepcopy(n)
	end
end

# fitness function for ising signal transduction
# 	( evotype istMGty, genotype gty, environment istEnv )
function fitness(istEnv::tCompEnv{<:Array{Float64}},gty::atSystemGty{<:atIsingMetaGty})::Float64
	fValues = zeros(Float64,gty.pMetaGty[1].prms.Ntrials)
	for t in 1:gty.pMetaGty[1].prms.Ntrials
		d2::Float64 = 0.0
		for iio in istEnv.idealInputOutput
			d2 += ( metropolis(gty.pMetaGty[1],broadcast(x->10^(x),gty.G),iio[1]) - iio[2] )^2
		end
		fValues[t] = exp(-sqrt(d2)/istEnv.selFactor)
	end
	return minimum(fValues) + (gty.pMetaGty[1].halfL - BASALFITNESS)
end

function fitness!(istEnv::tCompEnv{<:Array{Float64}},gty::atSystemGty{<:atIsingMetaGty})
	gty.pF[1] = fitness(istEnv,gty)
end


# **************************************
# * STATISTICS FUNCTIONS
# **************************************

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
	aan = Vector{Vector{Array{Int8,2}}}(undef,length(env.idealInputOutput))
	aSpinCorFisherz = Vector{Array{Float64,2}}(undef,length(env.idealInputOutput))

	for i in 1:length(env.idealInputOutput)
		aan[i] = [ Array{Int8}(undef, gty.pMetaGty[1].L, gty.pMetaGty[1].L) for ismpl in 1:prms.Nsmpl ]
		aSpinCorFisherz[i] = Array{Float64}(undef, gty.pMetaGty[1].L2, gty.pMetaGty[1].L2)

		metropolis!( gty.pMetaGty[1], broadcast(x->10^(x),gty.G), env.idealInputOutput[i][1], aan[i], prms )

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
		# sCor[si,sj] = tanh( mean( [ aSpinCorFisherz[i][si,sj] for i in 1:length(env.idealInputOutput)] ) )
		sCor[si,sj] = mean( [ aSpinCorFisherz[i][si,sj] for i in 1:length(env.idealInputOutput)] )
	end
end

function showJij!(gty::atSystemGty{<:atIsingMetaGty},JijMat::Array{Float64,2})
	ii::Int32 = 0
	for j in 1:2gty.pMetaGty[1].L, i in 1:2gty.pMetaGty[1].L
		JijMat[i,j] = j%2==1 ? ( i%2==0 ? 10^(gty.G[ii+=1]) : -1 ) : ( i%2==1 ? 10^(gty.G[ii+=1]) : 10^6+1 )
	end
end

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

export metropolis, metropolis!, fitness, fitness!, getSpinStat!, showJij!, getGStat!, getJijStat!


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
			print(f, aGty[i].pdG[1], "\t")
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
			print(f, aGty[i].pdG[1], "\t")
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
	open( "evoData_" * fileTag * ".dat", "w" ) do f
		for dataBatch in aData, t in 1:dataBatch.Ngen print(f,dataBatch.aveFitness[t],"\t") end, print(f,"\n")
		for dataBatch in aData, t in 1:dataBatch.Ngen print(f,dataBatch.growthFactor[t],"\t") end, print(f,"\n")
		for dataBatch in aData, t in 1:dataBatch.Ngen print(f,dataBatch.mutationFactor[t],"\t") end, print(f,"\n")
	end

	for dataBatch in aData
		for (gen, pop)  in zip(dataBatch.aGen, dataBatch.aLivingPop)
			write_tLivingPop( pop, fileTag * "_g$gen")
		end
	end
end

# to do: saving evolution parameters

export write_DTMCprm, write_MetaGty, write_aGty, write_tLivingPop, write_tEvoData

# read metagenotype
function read_MetaGty(L::Integer,prms::atMonteCarloPrm,fileTag::String)::atIsingMetaGty
	metaGtyMat = readdlm( fileTag * "_MetaGty" * ".dat" )
	return tIsingSigTransMGty(L,metaGtyMat[1],metaGtyMat[2],prms)
end

# read population and spits metagenotype and genotype
function read_aIsingSigTransAlphaGty(prms::atMonteCarloPrm,fileTag::String)
	gtyMat = readdlm( fileTag * "_atAlphaGty" * ".dat" )
	gMat = readdlm( fileTag * "_atAlphaGty_g" * ".dat" )
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
