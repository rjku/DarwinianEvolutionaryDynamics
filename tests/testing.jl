
N = 3
p = rand( [1,2], N )
b = [ rand(4) for i in 1:N ]

display( p )
display( b )
display( b ./ p )

using mEvoTypes

aGty = [ TabularGenotype( rand(1:10,2) ) for i in 1:3]
aGtyCopy = copy( aGty )
aGtyCopy2 = deepcopy( aGty )

aGty[1] = aGtyCopy[2]

display(aGty)
display(aGtyCopy)
display(aGtyCopy2)

testTensor = rand(5,3,5,3)

sum( sum(testTensor, dims=1), dims=2 ) ≈ sum( testTensor, dims=[1,2] )

display(sum( sum(testTensor, dims=1), dims=2 ))

struct testType
    x::Float64
end

import Base: +
+(X::testType,Y::testType) = testType( X.x + Y.x )

ciccio = +(testType(1),testType(2))

ciccio.x
