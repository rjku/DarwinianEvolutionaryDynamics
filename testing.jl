
struct ciccio{T<:AbstractArray}
	weight::T
end

struct Pasticcio{T<:ciccio}
	size::T
end

testVec = [1.0]

ac = [ ciccio(testVec), ciccio([testVec[1]]), ciccio(deepcopy(testVec)) ]
push!( ac, deepcopy(ac[1]) )

testVec[1] += 1.0

display(testVec)
display( [ ac[i].weight for i in eachindex(ac) ])

function feed_ciccio!(c)
	c = ciccio([10])
end

feed_ciccio!(c)
# display(ac[1].weight)

p = Pasticcio( ac[1] )

pcopy = deepcopy(p)

display(( p.size.weight[1], pcopy.size.weight[1]))
