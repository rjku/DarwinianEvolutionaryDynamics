
module mUtils

using JLD, HDF5

# ========== #
# STATISTICS #
# ========== #

function myCov(X::AbstractArray,XAve::Number,Y::AbstractArray,YAve::Number)
	N = length(X)
	N == length(Y) || throw(DimensionMismatch("inconsistent dimensions"))
	cov = 0.0
	for (x, y) in zip(X,Y)
		cov += (x - XAve)*(y - YAve)
	end
	return cov/(N-1)
end

function cumulativeCount(θ::Real,a::AbstractArray)
	N = 0
	for e in a
		if e >= θ
			N += 1
		end
	end
	return N/length(a)
end

function MPpdf(λ::Number,r::Number)
	λm = (1 - sqrt(r))^2;	λp = (1 + sqrt(r))^2
	return [ (λ >= λm) && (λ <= λp) ? sqrt( (λp - λ)*(λ - λm) )/( 2pi*r*λ ) : 0.0, λm, λp ]
end

export myCov, cumulativeCount, MPpdf

# ==== #
# MATH #
# ==== #

Base.log(a::AbstractArray{<:Real}) = map(log,a)
Base.log(b::Real,a::AbstractArray{<:Real}) = map(e->log(b,e),a)
log10(a::AbstractArray{<:Real}) = map(e->log(10.0,e),a)

Base.exp(a::AbstractArray{<:Real}) = map(exp,a)

# === #
# I/O #
# === #

function readJLD!(fileName,typeName,pointer)
	file = jldopen(fileName, "r")
	push!(pointer,read(file[typeName]))
	close(file)
end


# =========== #
# MISCELLANEA #
# =========== #

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

function get_binIndex(val::T,aBins::Vector{<:T}) where T
	for i in 1:length(aBins)-2
		if val >= aBins[i] && val < aBins[i+1]
			return i
		end
	end
	if val >= aBins[end-1] && val <= aBins[end]
		return length(aBins)-1
	else
		return 0
	end
end

export chopArray!, get_binIndex
export readJLD!

end
