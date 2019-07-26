
using Random, Plots

const SYSTEMSIZE = 40

# the (initial) state vector

# n = rand([-1,1], L, L)

# the (initial) interaction matrix
# Jij = randn(2*L^2)
interactionmatrix = ones(Float64,2*SYSTEMSIZE^2)

# ising 2d function simulator ( state n, system size L, interaction matrix Jij)
function ising2d(L::Int64, Jij::Array{Float64,1}) #n::Array{Int64,2},
	n = ones(Int64, L, L)

	β::Float64 = 1 # critical inverse temperature = 1/2.3
	Nsamplings::Int64 = 100
	L2::Int64 = L^2
	L3::Int64 = L^3
	Namcs::Int64 = 100*L2
	Γ = zeros(Int64,2)

	𝕁::Array{Float64,1} = Jij*β

	# magnetization vector
	m = zeros(Float64, Nsamplings)

	# useful nearest neighbour coordinate vectors
	Γp = collect(Int16, 2:L);	sizehint!(Γp,L); 	Γp = push!(Γp,1)
	Γm = collect(Int16, 1:L-1);	sizehint!(Γm,L); 	Γm = pushfirst!(Γm,L)

	for is in 1:Nsamplings
		for imcs in 1:Namcs
			i, j = rand(collect(1:L)), rand(collect(1:L))
			if rand() < ( 1. - n[i,j]*tanh( 𝕁[i+2(j-1)*L]*n[Γp[i],j] + 𝕁[Γm[i]+2(j-1)*L]*n[Γm[i],j] + 𝕁[i+(2j-1)*L]*n[i,Γp[j]] + 𝕁[i+(2*Γm[j]-1)*L]*n[i,Γm[j]] ) )/2
				n[i,j] = - n[i,j]
			end
		end
		m[is] = sum(n)/L2
	end

	pyplot() # Switch to using the PyPlot.jl backend
	return plot(collect(1:Nsamplings), m)
end

@time ising2d(SYSTEMSIZE,interactionmatrix)

# Juno.@run ising2d(SYSTEMSIZE,interactionmatrix)
