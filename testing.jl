
using Distances

a = rand(1:3, 5)
b = rand(1:3, 5)

display( Distances.hamming( a,b ) )

display( findall( e -> e == 5, a ) )

struct Ciccio
	height::Float64
	weight::Float64
end

function uno(T::Type,f1,f2)
	return T(f1,f2)
end

ac = uno(Ciccio,1.,2.)

push!(ac, Ciccio(2,1))

println(ac)

testVec = zeros(7)
testMat = zeros(7,4)

changer!(view(testMat,:,1))

testVec .+= 1

display(testVec)
display(testMat)

testArrVec = [ Float64[] for i in 1:5 ]
i = rand(1:5)
append!(testArrVec[i],Vector{Float64}(undef, 3))
display(testArrVec)

function fill!(vec)
	for i in eachindex(vec)
		vec[i] = 0
	end
end

fill!(@view testArrVec[i][end-2:end])
display(testArrVec)
