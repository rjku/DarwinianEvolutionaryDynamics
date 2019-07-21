
using Random, PyPlot

t1 = time_ns()

const L = 5

# the (initial) state vector
n = rand([-1,1], L, L)

# println(n)

function ising2d(initialState::Array{Int8,2}, systemSize::Int64, 𝕁::Array{Float64,1})
	β::Float64 = 1/4 # critical inverse temperature = 1/2.3
	Nsteps::Int64 = 100
	L2::Int64 = systemSize^2

	# useful nearest neighbour coordinate vectors
	Γp = collect(Int8, 2:L);	sizehint!(Γp,L); 	Γp = push!(Γp,1)
	Γm = collect(Int8, 1:L-1);	sizehint!(Γm,L); 	Γm = pushfirst!(Γm,10)

	W = Float32[ exp(-2J*β*r) for r in -4:4 ]

end


# magnetization vector
m = zeros(Float64, N)

DH = Int8
for n1 in 1:N
	for n2 in 1:N, j in 1:L, i in 1:L
		DH = n[i,j]*(n[Γp[i],j] + n[Γm[i],j] + n[i,Γp[j]] + n[i,Γm[j]]) + 5
		if rand() < W[DH]
			n[i,j] = - n[i,j]
		end
	end
	m[n1] = sum(n)/L2
end

t2 = time_ns()

println("elapsed time is: ", (t2-t1)/1.0e9, " seconds")

# time vector
t = collect(1:N)

plot(t, m)

# display(matshow(n))
