
struct ciccio
	weight::Vector{Float64}
end

testVec = [1.0]
testVec[1] += 1.0

ac = [ ciccio(testVec), ciccio(deepcopy(testVec)) ]

display(testVec)

display(ac[1].weight)
display(ac[2].weight)

# function feed_ciccio!(ac)
# 	for i in 1:1
# 		cp = ciccio([12])
# 		ac[2] = cp
# 	end
# end
#
# feed_ciccio!(ac)
# display(ac[1].weight)
