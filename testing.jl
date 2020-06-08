
using mUtils

struct Ciccio{Tad<:Vector{<:Dict{<:AbstractArray,Float64}}}
    aTb::Tad
end

function Ciccio(testArray,fn)
    table = Dict( e => fn(e) for e in testArray )
    Ciccio([table])
end

foo(x) = Float64(sum(x.^2))
testArray = [ [i,j] for i in 0:3, j in 0:3 ]

cic = Ciccio(testArray,foo)
display(cic.aTb[2][[2,2]])

using SparseArrays
A = sparse( [1,4,3], [2,7,2], [7.0,8.8,2.0], 10, 10)

display(nzrange(A,2))

testArray = [ rand(1:5,5), rand(1:5,5) ]

display( testArray )
display( testArray[2][ findall( e -> e == 1, testArray[1] ) ] )

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

using Revise
import mGraphs

testVec = rand(12)
display( testVec )

square = mGraphs.FiniteSquareLattice(3)
weighted = mGraphs.EdgeWeightedSquareLattice(3, testVec)

display( mGraphs.neighbors(square,1) )
display( mGraphs.neighbors(weighted,5) )

display( Array(mGraphs.transitionMatrix(square)) )

maximum([1,2,3])

testVec = rand(3)
testMat = rand([0,1],3,3)

display( testMat )
display( testMat * testVec )
