
N = 3
p = rand( [1,2], N )
b = [ rand(4) for i in 1:N ]

display( p )
display( b )
display( b ./ p )

using EvolutionaryDynamics

aGty = [ TabularGenotype( rand(1:10,2) ) for i in 1:3]
aGtyCopy = copy( aGty )
aGtyCopy2 = deepcopy( aGty )

aGty[1] = aGtyCopy[2]

display(aGty)
display(aGtyCopy)
display(aGtyCopy2)

testTensor = rand(5,3,5,3)

display( testTensor )

display( sum(testTensor, dims=2)[1,1,1,1] )

sum( sum(testTensor, dims=1), dims=2 ) ≈ sum( testTensor, dims=[1,2] )

display(sum( sum(testTensor, dims=1), dims=2 ))

testVec1 = rand(10)
testVec2 = rand(10)

for (i,v) in enumerate(testVec1), (j,w) in enumerate(testVec2)
    println(i, " ", j)
end

for d in 1:length( size( testTensor ) )
    for i in 1:size( testTensor )[d]
        println( i )
    end
end

for i in 1:size(testTensor)[1], j

testMat = Matrix{Int64}(undef, 3, 3)

for i in 1:3, j in 1:3
    testMat[i,j] = i + (j-1)*3
end

display( testMat )

testMat = rand(4,4)

testVec = [1,2,3,4]

display( testMat )
display( testVec .^ 2 )
display( testMat .* testVec' )
display( sum( testVec .* testMat, dims=1 ) )
display( testVec' * testVec )
display( testVec .* testMat .* testVec' )
display( sum(testMat, dims=1) )
display( testVec .* testMat ./ sum( testVec .* testMat, dims=1 ) )

display( [ (testVec[i] - testVec[j])/testVec[i] for i in 1:2, j in 1:4 ] )

display( sum( testVec .* testMat, dims=1 ) )

display( typeof(testVec) )
