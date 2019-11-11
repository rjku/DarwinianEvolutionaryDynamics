
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

function chopArray!(a::Array{<:Complex}, δ::Float64=10^-12)
	for (i, e) in enumerate(a)
		if abs(real(e)) < δ && abs(imag(e)) < δ
			a[i] = 0.0 + im*0.0
		elseif abs(real(e)) < δ && abs(imag(e)) >= δ
			a[i] = 0.0 + im*imag(e)
		elseif abs(real(e)) >= δ && abs(imag(e)) < δ
			a[i] = real(e) + im*0.0
		end
	end
end

function MPpdf(λ::Number,r::Number)
	λm = (1 - sqrt(r))^2;	λp = (1 + sqrt(r))^2

	return [ (λ >= λm) && (λ <= λp) ? sqrt( (λp - λ)*(λ - λm) )/( 2pi*r*λ ) : 0.0, λm, λp ]
end

export myCov, chopArray!, MPpdf

end
