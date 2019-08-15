
module isingFunc
export tIsingST, metropolis

using Random, PyPlot

# ===================
# type: ising signal transduction
struct tIsingST
	L::Int32				# system size
	β::Float32 				# inverse temperature [ critical inverse temperature ≃ 1/2.3 ≃ 0.43 ]
	he::Float32 			# global external field H

	L2::Int32;	halfL::Int32
	li::Int32;	li2::Int32			# input/redout region size set to one tenth of the system size
	ℍe::Float32

	# definition: useful nearest neighbour coordinate vectors
	jp::Array{Int32,1}; jm::Array{Int32,1}
	Jpi::Array{Int32,2}; Jmi::Array{Int32,2}; Jpj::Array{Int32,2}; Jmj::Array{Int32,2}
end

# ===================
# constructor: ising signal transduction
tIsingST(L, β, he) = tIsingST(
	L, β, he, L^2, L÷2, L÷10, (L÷10)^2, he*β,
	Int32[i%L+1 for i in 1:L],
	Int32[(L-2+i)%L+1 for i in 1:L],
	Int32[i+2(j-1)*L for i in 1:L, j in 1:L],
	Int32[(L-2+i)%L+1+2(j-1)*L for i in 1:L, j in 1:L],
	Int32[i+(2j-1)*L for i in 1:L, j in 1:L],
	Int32[i+(2(L-2+j)%L+1)*L for i in 1:L, j in 1:L]
)

# ===================
# MonteCarlo simulation function for ising signal transduction: readout magnetization evaluation
# 	( tIsingST isingST, state n, interaction matrix Jij, input field h ) → readout magnetization m
function metropolis(isingST::tIsingST, n::Array{Int8,2}, Jij::Array{T,1}, hi::AbstractFloat)::T where {T<:AbstractFloat}
	Nmcs::Int32 =  50*isingST.L2		# number of Monte Carlo steps
	Nsamplings::Int16 = 20				# number of samplings points
	𝕁::Array{Float32,1} = Jij*isingST.β;	ℍi::Float32 = hi*isingST.β

	# definition: time-averaged readout magnetization
	mro::Float32 = 0.0

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

# MonteCarlo simulation function for ising signal transduction: time-averaged spin config evaluation
# 	( tIsingST isingST, state n, interaction matrix Jij, input field h ) → readout magnetization m
function metropolis(isingST::tIsingST, n::Array{Int8,2}, Jij::Array{T,1}, hi::AbstractFloat, aves::Array{T,2}) where {T<:AbstractFloat}
	Nmcs::Int32 =  50*isingST.L2		# number of Monte Carlo steps
	Nsamplings::Int16 = 20				# number of samplings points
	𝕁::Array{Float32,1} = Jij*isingST.β;	ℍi::Float32 = hi*isingST.β

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

end

SYSTEMSIZE = 30

# the (initial) state vector
# n = rand([-1,1], 	L, L)
n = ones(Int8, SYSTEMSIZE, SYSTEMSIZE)

# the (initial) interaction matrix
# interactionmatrix = ones(Float32,2*SYSTEMSIZE^2)
interactionmatrix = rand(Float32,2*SYSTEMSIZE^2)

isingST = tIsingST(SYSTEMSIZE,0.7,0.5)

aves = Array{Float32}(undef, SYSTEMSIZE, SYSTEMSIZE)

@time println(metropolis(isingST,n,interactionmatrix,-10.0))
@time metropolis(isingST,n,interactionmatrix,-10.0,aves)
# Juno.@run metropolis(SYSTEMSIZE,interactionmatrix)
