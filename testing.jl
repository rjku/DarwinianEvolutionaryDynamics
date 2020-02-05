
function changer!(vec)
	for i in 1:length(vec)
		vec[i] = rand()
	end
end

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
