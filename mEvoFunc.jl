
module mEvoFunc
using mEvoTypes, Random, Statistics, Base.Threads, DelimitedFiles, Dates, Distances, ProgressMeter
import Future

const FITNESSOFFSET, FITNESSTHRESHOLD, BASALFITNESS = .001, .95, 3.0

# initializer constructor for discrete time evolution. Same MetaGenotype for all.
function initLivingPop( N::Int32,ety::Tety,env::Tenv,aMGty::Vector{Tmgty},aGty::Vector{Tgty} ) where {
		Tety<:atEvotype,Tenv<:atEnvironment,Tmgty<:atMetaGenotype,Tgty<:atGenotype }
	for i in 1:N fitness!(env,aGty[i]) end
	return tLivingPop{Tety,Tenv,Vector{Tmgty},Vector{Tgty}}( Int32[N,N,length(aGty)],ety,env,aMGty,aGty )
end

tEty{Tx}(repRate::Float64,mutRate::Float64,ΔtOffset::Float64,dX::Int32,Xvar::Tx) where {Tx<:Number} =
	tEty{Tx}( repRate,mutRate,ΔtOffset,[repRate/(2dX*mutRate+ΔtOffset)],[mutRate/(2dX*mutRate+ΔtOffset)],Xvar )

# function changing rep and mut -factors
function set_tEtyFactors(ety::tEty,dX::Int32)
	ety.pRepFactor[1] = ety.repRate/(2dX*ety.mutRate+ety.ΔtOffset)
	ety.pMutFactor[1] = ety.mutRate/(2dX*ety.mutRate+ety.ΔtOffset)
end

export tEty, initLivingPop


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
				pop.aGty[pop.pN[1]] = deepcopy(pop.aGty[i]) 		# <- check whether is necessary to deepcopy
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

function upgradeGtyX!(gty::tVecGty{Vector{TMGty},Vector{Tx}},Δx::Tx) where {TMGty<:atMetaGenotype,Tx}
	Xref = copy(gty.X)
	append!( gty.X, ones(Float64, 8gty.pMetaGty[1].L+8) )
	# append!( gty.X, rand(-4:.5:5, 4gty.pMetaGty[1].L+8) )
	# gty.X .= 1.																	# <- get rid of this

	for u in 0:1, j in 1:gty.pMetaGty[1].L, k in 0:1, i in 1:gty.pMetaGty[1].halfL
		gty.X[ i + k*(gty.pMetaGty[1].halfL + 1) + (j + u*(gty.pMetaGty[1].L + 2) - 1)*(gty.pMetaGty[1].L+2) ] =
			Xref[ i + k*gty.pMetaGty[1].halfL + (j + u*gty.pMetaGty[1].L - 1)*gty.pMetaGty[1].L ]
	end

	# duplication of previous line genotype + some randomness

	for k in 0:1, j in 1:gty.pMetaGty[1].L, i in 1:2
		gty.X[i*(gty.pMetaGty[1].halfL+1)+(j+k*(gty.pMetaGty[1].L+2)-1)*(gty.pMetaGty[1].L+2)] =
			Xref[i*gty.pMetaGty[1].halfL+(j+k*(gty.pMetaGty[1].L)-1)*gty.pMetaGty[1].L] + rand(-1:1)*gty.pMetaGty[1].Δx
	end

	for u in 0:1, k in 0:1, j in 1:2, i in 1:gty.pMetaGty[1].halfL
		gty.X[i+k*(gty.pMetaGty[1].halfL+1)+(j+gty.pMetaGty[1].L+u*(gty.pMetaGty[1].L+2)-1)*(gty.pMetaGty[1].L+2)] =
			Xref[i+k*gty.pMetaGty[1].halfL+(j+gty.pMetaGty[1].L*(u+1)-3)*gty.pMetaGty[1].L] + rand(-1:1)*Δx
	end
end

# function evoUpgrade!(pop::tLivingPop{<:atEvotype,<:atEnvironment,Vector{tIsingSigTransMGty{<:atMonteCarloPrm}},Vector{<:atVecGty}})
function evoUpgrade!(pop::tLivingPop)
	for i in 1:pop.pN[2]
		if pop.aGty[i].pF[1] - (pop.aGty[i].pMetaGty[1].halfL - BASALFITNESS) > FITNESSTHRESHOLD

			# update the genotypic variables X
			upgradeGtyX!(pop.aGty[i],pop.ety.Xvar)

			# update the metagenotype
			foundMetaGty = false
			for metaGty in pop.aMetaGty
				if pop.aGty[i].pMetaGty[1].L+2 == metaGty.L
					foundMetaGty = true
					pop.aGty[i].pMetaGty[1] = metaGty
				end
			end
			if !foundMetaGty
				push!( pop.aMetaGty,tIsingSigTransMGty(
					pop.aGty[i].pMetaGty[1].L+Int32(2), pop.aGty[i].pMetaGty[1].β, pop.aGty[i].pMetaGty[1].he, pop.aGty[i].pMetaGty[1].prms
					))
				pop.aGty[i].pMetaGty[1] = pop.aMetaGty[end]
				set_tEtyFactors(pop.ety,pop.aGty[i].pMetaGty[1].dX)
			end

			# update the fitness
			fitness!(pop.env,pop.aGty[i])
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
		evoUpgrade!(pop)
		evo.aveFitness[gen] = mean( [pop.aGty[i].pF[1] for i in 1:pop.pN[2]] )

		if evo.aveFitness[gen] >= evo.pAveFt[1]
			push!(evo.aLivingPop,deepcopy(pop))
			evo.pAveFt[1] += evo.aveFtinc
			println("\nEvolutionary achievement at generation $gen: ⟨f⟩ = ", evo.aveFitness[gen] )
		end
	end
end

export replication!, effMutation!, effSelection!, evoUpgrade!, evolution!, upgradeGtyX!


# ***********
# *  ISING  *
# ***********

# function: flipping --- for metropolis
function flipping!(istMGty::tIsingSigTransMGty, βJij::Array{<:Real,1}, βhi::Real, n::Array{<:Integer,2})
	# Monte Carlo step evaluated using Glauber transition rates:
	# next possibly transitioning spin (coordinates)
	i, j = rand(collect(1:istMGty.L)), rand(collect(1:istMGty.L))
	# transitioning?!
	if i > istMGty.li || j > istMGty.li
		if rand() < ( 1. - n[i,j]*tanh(
				βJij[istMGty.Jpi[i,j]]*n[istMGty.jp[i],j] + βJij[istMGty.Jmi[i,j]]*n[istMGty.jm[i],j] +
				βJij[istMGty.Jpj[i,j]]*n[i,istMGty.jp[j]] + βJij[istMGty.Jmj[i,j]]*n[i,istMGty.jm[j]] +
				istMGty.βhe ))/2
			n[i,j] = - n[i,j]
		end
	# else
	# 	if rand() < ( 1. - n[i,j]*tanh(
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
	n = Int8(istMGty.he != 0.0 ? sign(istMGty.he) : 1).*ones(Int8, istMGty.L, istMGty.L)
	n[1:istMGty.li,1:istMGty.li] .= sign(hi)

	# definition: time-averaged readout magnetization
	mro::Float64 = 0.0

	for is in 1:istMGty.prms.Nsmpl
		for imcs in 1:istMGty.prms.Nmcsps, ilp in 1:istMGty.L2
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
function fitness!(istEnv::tCompEnv{<:Array{Float64}},gty::tVecGty)
	fValues = zeros(Float64,gty.pMetaGty[1].prms.Ntrials)
	for t in 1:gty.pMetaGty[1].prms.Ntrials
		d2::Float64 = 0.0
		for iio in istEnv.idealInputOutput
			d2 += ( metropolis(gty.pMetaGty[1],broadcast(x->10^(x),gty.X),iio[1]) - iio[2] )^2
		end
		fValues[t] = exp(-sqrt(d2)/istEnv.selFactor)
	end
	gty.pF[1] = minimum(fValues) + (gty.pMetaGty[1].halfL - BASALFITNESS)
end

function fitness(istEnv::tCompEnv{<:Array{Float64}},gty::tVecGty)::Float64
	fValues = zeros(Float64,gty.pMetaGty[1].prms.Ntrials)
	for t in 1:gty.pMetaGty[1].prms.Ntrials
		d2::Float64 = 0.0
		for iio in istEnv.idealInputOutput
			d2 += ( metropolis(gty.pMetaGty[1],broadcast(x->10^(x),gty.X),iio[1]) - iio[2] )^2
		end
		fValues[t] = exp(-sqrt(d2)/istEnv.selFactor)
	end
	return minimum(fValues) + (gty.pMetaGty[1].halfL - BASALFITNESS)
end

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
function getSpinStat!(env::tCompEnv,gty::tVecGty,aSpinAve::Array{Array{Float64,2},1},
		aSpinCov::Array{Array{Float64,2},1},sCor::Array{Float64,2},prms::tDTMCprm)
	aan = Vector{Vector{Array{Int8,2}}}(undef,length(env.idealInputOutput))
	aSpinCorFisherz = Vector{Array{Float64,2}}(undef,length(env.idealInputOutput))

	for i in 1:length(env.idealInputOutput)
		aan[i] = [ Array{Int8}(undef, gty.pMetaGty[1].L, gty.pMetaGty[1].L) for ismpl in 1:prms.Nsmpl ]
		aSpinCorFisherz[i] = Array{Float64}(undef, gty.pMetaGty[1].L2, gty.pMetaGty[1].L2)

		metropolis!( gty.pMetaGty[1], broadcast(x->10^(x),gty.X), env.idealInputOutput[i][1], aan[i], prms )

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

function showJij!(gty::tVecGty,JijMat::Array{Float64,2})
	ii::Int32 = 0
	for j in 1:2gty.pMetaGty[1].L, i in 1:2gty.pMetaGty[1].L
		JijMat[i,j] = j%2==1 ? ( i%2==0 ? 10^(gty.X[ii+=1]) : -1 ) : ( i%2==1 ? 10^(gty.X[ii+=1]) : 10^6+1 )
	end
end

function getXStat!(pop::tLivingPop,XAve::Vector{Float64},XCov::Array{Float64,2},XCor::Array{Float64,2})
	mapreduce( x -> x ==  pop.aGty[1].pMetaGty[1].dX, &, [ pop.aGty[i].pMetaGty[1].dX for i in 2:pop.pN[2] ] ) ||
		throw(DimensionMismatch("inconsistent dimensions"))

	for x in 1:pop.aGty[1].pMetaGty[1].dX
		XAve[x] = mean( [pop.aGty[i].X[x] for i in 1:pop.pN[2]] )
	end

	for y in 1:pop.aGty[1].pMetaGty[1].dX, x in 1:pop.aGty[1].pMetaGty[1].dX
		XCov[x,y] = myCov( [ pop.aGty[i].X[x] for i in 1:pop.pN[2] ], XAve[x], [ pop.aGty[i].X[y] for i in 1:pop.pN[2] ], XAve[y] )
	end

	for y in 1:pop.aGty[1].pMetaGty[1].dX, x in 1:pop.aGty[1].pMetaGty[1].dX
		XCor[x,y] = XCov[x,y]/sqrt(XCov[x,x]*XCov[y,y])
	end
end

function getJijStat!(pop::tLivingPop,JijAve::Vector{Float64},JijCov::Array{Float64,2},JijCor::Array{Float64,2})
	mapreduce( x -> x ==  pop.aGty[1].pMetaGty[1].dX, &, [ pop.aGty[i].pMetaGty[1].dX for i in 2:pop.pN[2] ] ) ||
		throw(DimensionMismatch("inconsistent dimensions"))

	for x in 1:pop.aGty[1].pMetaGty[1].dX
		JijAve[x] = mean( [ 10^pop.aGty[i].X[x] for i in 1:pop.pN[2]] )
	end

	for y in 1:pop.aGty[1].pMetaGty[1].dX, x in 1:pop.aGty[1].pMetaGty[1].dX
		JijCov[x,y] = myCov( [ 10^pop.aGty[i].X[x] for i in 1:pop.pN[2] ], JijAve[x], [ 10^pop.aGty[i].X[y] for i in 1:pop.pN[2] ], JijAve[y] )
	end

	for y in 1:pop.aGty[1].pMetaGty[1].dX, x in 1:pop.aGty[1].pMetaGty[1].dX
		JijCor[x,y] = JijCov[x,y]/sqrt(JijCov[x,x]*JijCov[y,y])
	end
end

function getJij!(gty::tVecGty,JijMat::Array{Float64,2})
	JijMat .= 0
	for x in 1:gty.pMetaGty[1].L, y in 1:gty.pMetaGty[1].L
		JijMat[ x+(y-1)*gty.pMetaGty[1].L, gty.pMetaGty[1].jp[x]+(y-1)*gty.pMetaGty[1].L ] = 10.0^gty.X[ gty.pMetaGty[1].Jpi[x,y] ]
		JijMat[ x+(y-1)*gty.pMetaGty[1].L, gty.pMetaGty[1].jm[x]+(y-1)*gty.pMetaGty[1].L ] = 10.0^gty.X[ gty.pMetaGty[1].Jmi[x,y] ]
		JijMat[ x+(y-1)*gty.pMetaGty[1].L, x+(gty.pMetaGty[1].jp[y]-1)*gty.pMetaGty[1].L ] = 10.0^gty.X[ gty.pMetaGty[1].Jpj[x,y] ]
		JijMat[ x+(y-1)*gty.pMetaGty[1].L, x+(gty.pMetaGty[1].jm[y]-1)*gty.pMetaGty[1].L ] = 10.0^gty.X[ gty.pMetaGty[1].Jmj[x,y] ]
	end
end

export metropolis, metropolis!, fitness, fitness!, getSpinStat!, showJij!, getXStat!, getJijStat!


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
	open( fileTag * "_aGty" * ".dat", "w" ) do f
		for i in 1:Npop
			print(f, aGty[i].pMetaGty[1].dX, "\t")
			for x in aGty[i].X
				print(f, x, "\t")
			end
			print(f, aGty[i].pF[1], "\n")
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
end

# to do: saving evolution parameters

export write_DTMCprm, write_MetaGty, write_aGty, write_tLivingPop, write_tEvoData

# read metagenotype
function read_MetaGty(L::Integer,prms::atMonteCarloPrm,fileTag::String)::atIsingMetaGty
	metaGtyMat = readdlm( fileTag * "_MetaGty" * ".dat" )
	return tIsingSigTransMGty(L,metaGtyMat[1],metaGtyMat[2],prms)
end

# read population and spits metagenotype and genotype
function read_aIsingSigTransGty(prms::atMonteCarloPrm,fileTag::String)
	gtyMat = readdlm( fileTag * "_aGty" * ".dat" )
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
	aGty = [ tVecGty([aMetaGty[collect(1:length(sysSizes))[sysSizes .== L[i]][1]]],
		gtyMat[i,2:1+Int32(gtyMat[i,1])] , Float64[gtyMat[i,end]] ) for i in 1:Npop ]
	return aMetaGty, aGty, Npop
end

export read_aIsingSigTransGty

# *******************
# TRIVIAL
# *******************

function fitness!(trivialEty::tTrivialEty,trivialEnv::tTrivialEnv,gty::tVecGty{Array{T,1}}) where {T<:Real}
	gty.pF[1]=1/(euclidean(gty.X,ones(Float64,gty.pMetaGty[1].dX))+FITNESSOFFSET)
end

export fitness!

end
