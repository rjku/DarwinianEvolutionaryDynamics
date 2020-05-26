
using mUtils

testMat = [ [1,2] [5,2] ]

N = 3

testMat = rand([-1,0,1],N,N)
display(testMat)

testVec = rand(N)
display( testVec )

x = testMat \ testVec

display( sum(testMat .* x',dims=2) ≈ testVec )

display( testMat .* testVec' .* testMat )

display( transpose(transpose(testMat) .* testVec) )

display( transpose(transpose(testMat) .* testVec) )

newTestVec = testMat[3,:] .* testVec

display( newTestVec )

display(ReLog(testMat))


map(e->log(e),testVec)
map(log,testVec)

function main(args)
    @show parse(Int64,args[1])
end

main(ARGS)

a = Array{Matrix{Float64}}(undef, 2)

a[1] = rand(2,2)
a[2] = rand(3,4)

display(a)

folderVar = Base.run(`ls`)

typeof(folderVar)

folderVar
