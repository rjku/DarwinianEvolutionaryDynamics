
function popStatistics(pop::Population)
	aG = [ pop.aGty[i].G[1] for i in 1:pop.pN[2] ]
	p = zeros( pop.ety.G.Nv )
	for g in aG
		p[g] += 1
	end
	return p/pop.pN[2]
end

function popStatistics(aGty::Vector{<:AbstractGenotype},G::AbstractGraph)
	aG = [ gty.G[1] for gty in aGty ]
	p = zeros( G.Nv )
	for g in aG
		p[g] += 1
	end
	return p/length(aGty)
end

function popStatisticsVaryingTabularI(aGty::Vector{<:AbstractGenotype},G::AbstractGraph)
	aG = [ gty.G for gty in aGty ]

	p = zeros( G.Nv )
	b = [ zeros(4) for i in 1:G.Nv ]

	for g in aG
		p[g[1]] += 1
		if g[2] > 0
			b[g[1]][g[2]] += 1
		end
	end
	return p/length(aGty), b ./ p
end

export popStatistics
