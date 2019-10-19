
module mUtils

function myCov(X::AbstractArray,XAve::Number,Y::AbstractArray,YAve::Number)
	N = length(X)
	N == length(Y) || throw(DimensionMismatch("inconsistent dimensions"))
	cov = 0.0
	for (x, y) in zip(X,Y)
		cov += (x - XAve)*(y - YAve)
	end
	return cov/(N-1)
end

function chopArray!(a::AbstractArray, δ::Float64=10^-12)
	for (i, e) in enumerate(a)
		if abs(e) < δ
			a[i] = 0.0
		end
	end
end

end
