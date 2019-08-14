
using Random, PyPlot#, Plots

SYSTEMSIZE = 30

# the (initial) state vector
# n = rand([-1,1], 	L, L)

# the (initial) interaction matrix
# interactionmatrix = ones(Float32,2*SYSTEMSIZE^2)
interactionmatrix = rand(Float32,2*SYSTEMSIZE^2)

# ising 2d function simulator ( state n, system size L, interaction matrix Jij, input field h) → readout magnetization m
function ising2d(L::Integer, Jij::Array{T,1}, hi::AbstractFloat)::AbstractFloat where {T<:AbstractFloat}
	n = ones(Int16, L, L)

	L2::Int16 = L^2;	halfL::Int16 = L÷2

	β::Float32 = 0.7 			# critical inverse temperature ≃ 1/2.3 ≃ 0.43
	he::Float32 = 0.5 			# global external field H
	Nmcs::Int32 =  50*L2		# number of Monte Carlo steps
	Nsamplings::Int16 = 20		# number of samplings points

	li::Int16 = L÷10;	li2::Int16 = li^2			# input/redout region size set to one tenth of the system size

	𝕁::Array{Float32,1} = Jij*β;	ℍe::Float32 = he*β;		ℍi::Float32 = hi*β

	# definition: time-averaged readout magnetization, time-averaged spin config, and overall magnetization
	mro::Float32 = 0.0
	# aves = zeros(Int16, L, L)
	# m = zeros(Float32, Nsamplings)

	# definition: useful nearest neighbour coordinate vectors
	jp = collect(Int16, 2:L);	sizehint!(jp,L); 	jp = push!(jp,1)
	jm = collect(Int16, 1:L-1);	sizehint!(jm,L); 	jm = pushfirst!(jm,L)

	Jpi = Int16[i+2(j-1)*L for i in 1:L, j in 1:L]
	Jmi = Int16[jm[i]+2(j-1)*L for i in 1:L, j in 1:L]
	Jpj = Int16[i+(2j-1)*L for i in 1:L, j in 1:L]
	Jmj = Int16[i+(2*jm[j]-1)*L for i in 1:L, j in 1:L]

	for is in 1:Nsamplings
		for imcs in 1:Nmcs
			# Monte Carlo step evaluated using Glauber transition rates:
			# next possibly transitioning spin (coordinates)
			i, j = rand(collect(1:L)), rand(collect(1:L))
			# transitioning?!
			if i <= li && j <= li
				if rand() < ( 1. - n[i,j]*tanh( 𝕁[Jpi[i,j]]*n[jp[i],j] + 𝕁[Jmi[i,j]]*n[jm[i],j] +
						𝕁[Jpj[i,j]]*n[i,jp[j]] + 𝕁[Jmj[i,j]]*n[i,jm[j]] + ℍe + ℍi ))/2
					n[i,j] = - n[i,j]
				end
			else
				if rand() < ( 1. - n[i,j]*tanh( 𝕁[Jpi[i,j]]*n[jp[i],j] + 𝕁[Jmi[i,j]]*n[jm[i],j] +
						𝕁[Jpj[i,j]]*n[i,jp[j]] + 𝕁[Jmj[i,j]]*n[i,jm[j]] + ℍe ))/2
					n[i,j] = - n[i,j]
				end
			end
		end

		# evaluation: time-averaged readout magnetization, time-averaged spin config, and overall magnetization
		mro += sum(n[halfL+1:halfL+li,halfL+1:halfL+li])
		# aves += n
		# m[is] = sum(n)/L2
	end

	# aves /= Nsamplings
	mro /= Nsamplings * li2

	# println(n)

	# pyplot() # Switch to using the PyPlot.jl backend
	# return plot(collect(1:Nsamplings), m)
	# matshow(aves,cmap="Greys_r"); gcf()
	# matshow(n,cmap="Greys_r"); gcf()

	return mro
end

@time println(ising2d(SYSTEMSIZE,interactionmatrix,-10.0))
# @time ising2d(SYSTEMSIZE,interactionmatrix,-10.0)

# Juno.@run ising2d(SYSTEMSIZE,interactionmatrix)
