
module mIsingEvolution
using mAbstractTypes, mEvolution, Random, PyPlot

export tIsingSigTransEty, tCompEnv, tVecGty, metropolis, fitness

# ===================
# type: ising signal transduction
struct tIsingSigTransEty <: atEvotype
	L::Int32				# system size
	lJij::Int32				# interaction matrix number of entries

	β::Float64 				# inverse temperature [ critical inverse temperature ≃ 1/2.3 ≃ 0.43 ]
	he::Float64 			# global external field H

	L2::Int32;	halfL::Int32
	li::Int32;	li2::Int32			# input/redout region size set to one tenth of the system size
	ℍe::Float64

	# definition: useful nearest neighbour coordinate vectors
	jp::Array{Int32,1}; jm::Array{Int32,1}
	Jpi::Array{Int32,2}; Jmi::Array{Int32,2}; Jpj::Array{Int32,2}; Jmj::Array{Int32,2}
end

# constructor: ising signal transduction
tIsingSigTransEty(L::Int32, β::Real, he::Real) = tIsingSigTransEty(
	L, 2L^2, β, he, L^2, L÷2, L÷10, (L÷10)^2, he*β,
	Int32[i%L+1 for i in 1:L],
	Int32[(L-2+i)%L+1 for i in 1:L],
	Int32[i+2(j-1)*L for i in 1:L, j in 1:L],
	Int32[(L-2+i)%L+1+2(j-1)*L for i in 1:L, j in 1:L],
	Int32[i+(2j-1)*L for i in 1:L, j in 1:L],
	Int32[i+(2(L-2+j)%L+1)*L for i in 1:L, j in 1:L]
)

# ===================
# MonteCarlo simulation function for ising signal transduction: readout magnetization evaluation
# 	( tIsingST isingST, interaction matrix Jij, input field h ) → readout magnetization m
function metropolis(isingST::tIsingSigTransEty, Jij::Array{T,1}, hi::Real)::Real where {T<:Real}
	Nmcs::Int32 =  50*isingST.L2		# number of Monte Carlo steps
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
	Nmcs::Int32 =  50*isingST.L2		# number of Monte Carlo steps
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
	Nmcs::Int32 =  50*isingST.L2		# number of Monte Carlo steps
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
	gty.F[1] = 1/(sqrt(invf)+0.1)	# fitness offset = 0.1
end

# constructor of tLivingPop based on ising fitness
function tLivingPop{T}(N::Int32,ety::tIsingSigTransEty,env::atEnv,aGty::Array{<:atGenotype{T},1}) where {T}
	for i in 1:N fitness!(ety,env,aGty[i]) end
	tLivingPop{T}( Int32[N,N,length(aGty)],ety,env,aGty )
end

end
