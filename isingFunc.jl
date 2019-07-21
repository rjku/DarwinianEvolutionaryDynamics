
using Random, Plots

const L = 20

# the (initial) state vector
n = rand([-1,1], L, L)
Jij = randn(2*L^2)

show(Jij)

function ising2d(state::Array{Int64,2}, systemSize::Int64, interactionMatrix::Array{Float64,1})
	β::Float64 = 1/4 # critical inverse temperature = 1/2.3
	Nsteps::Int64 = 100000
	L2::Int64 = systemSize^2
	L3::Int64 = systemSize^3
	Γ = zeros(Int64,2)

	𝕁::Array{Float64,1} = interactionMatrix*β
	m = zeros(Float64, Nsteps)

	# useful nearest neighbour coordinate vectors
	Γp = collect(Int8, 2:L);	sizehint!(Γp,L); 	Γp = push!(Γp,1)
	Γm = collect(Int8, 1:L-1);	sizehint!(Γm,L); 	Γm = pushfirst!(Γm,L)

	for n1 in 1:Nsteps
		for n2 in 1:L3
			i, j = rand(collect(1:L)), rand(collect(1:L))
			if rand() < .5 ( 1. - state[i,j]*tanh( 𝕁[i+2(j-1)*systemSize]*n[Γp[i],j] + 𝕁[Γm[i]+2(j-1)*systemSize]*n[Γm[i],j] + 𝕁[i+(2j-1)*systemSize]*n[i,Γp[j]] + 𝕁[i+(2Γm[j]-1)*systemSize]*n[i,Γm[j]] ) )
				n[i,j] = - n[i,j]
			end
		end
		m[n1] = sum(n)/L2
	end

	pyplot() # Switch to using the PyPlot.jl backend
	return plot(collect(1:Nsteps), m)
end

@time ising2d(n,L,Jij)

println(n)
