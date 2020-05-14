
testVec = rand(5)
display(testVec)

Base.log(a::AbstractArray{<:Real}) = map(log,a)
Base.log(b::Real,a::AbstractArray{<:Real}) = map(e->log(b,e),a)
log10(a::AbstractArray{<:Real}) = map(e->log(10.0,e),a)

display(log10(testVec))

map(e->log(e),testVec)
map(log,testVec)

function main(args)
    @show parse(Int64,args[1])
end

main(ARGS)
