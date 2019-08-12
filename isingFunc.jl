
using Random, PyPlot#, Plots

const SYSTEMSIZE = 40

# the (initial) state vector

# n = rand([-1,1], 	L, L)

# the (initial) interaction matrix
# Jij = randn(2*L^2)
interactionmatrix = ones(Float64,2*SYSTEMSIZE^2)

# ising 2d function simulator ( state n, system size L, interaction matrix Jij)
function ising2d(L::Int64, Jij::Array{Float64,1}) #n::Array{Int64,2},
	n = ones(Int64, L, L)

	β::Float64 = 0.7 			# critical inverse temperature ≃ 1/2.3 ≃ 0.43
	L2::Int64 = L^2
	L3::Int64 = L^3
	Nmcs::Int64 = 100*L2		# number of Monte Carlo steps
	Nsamplings::Int64 = 100		# number of samplings points

	𝕁::Array{Float64,1} = Jij*β

	# magnetization vector
	m = zeros(Float64, Nsamplings)

	# average spin configuration
	aves = zeros(Int64, L, L)

	# useful nearest neighbour coordinate vectors
	jp = collect(Int16, 2:L);	sizehint!(jp,L); 	jp = push!(jp,1)
	jm = collect(Int16, 1:L-1);	sizehint!(jm,L); 	jm = pushfirst!(jm,L)

	for is in 1:Nsamplings
		for imcs in 1:Nmcs
			# Monte Carlo step evaluated using Glauber transition rates:
			# next possibly transitioning spin (coordinates)
			i, j = rand(collect(1:L)), rand(collect(1:L))
			# transitioning?!
			if rand() < ( 1. - n[i,j]*tanh( 𝕁[i+2(j-1)*L]*n[jp[i],j] + 𝕁[jm[i]+2(j-1)*L]*n[jm[i],j] + 𝕁[i+(2j-1)*L]*n[i,jp[j]] + 𝕁[i+(2*jm[j]-1)*L]*n[i,jm[j]] ) )/2
				n[i,j] = - n[i,j]
			end
		end
		aves += n

		m[is] = sum(n)/L2
	end

	aves /= Nsamplings

	# pyplot() # Switch to using the PyPlot.jl backend
	# return plot(collect(1:Nsamplings), m)
	matshow(aves,cmap="Greys_r"); gcf()
end

@time ising2d(SYSTEMSIZE,interactionmatrix)

# Juno.@run ising2d(SYSTEMSIZE,interactionmatrix)
