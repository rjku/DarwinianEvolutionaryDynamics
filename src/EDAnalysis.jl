
function popStatistics(pop::Population)
	aG = [ pop.aGty[i].G[1] for i in 1:pop.pN[2] ]
	p = zeros( pop.ety.graph.Nv )
	for g in aG
		p[g] += 1
	end
	return p/pop.pN[2]
end

function popStatistics(aGty::Vector{<:AbstractGenotype},G::AbstractGraph)
	aG = [ gty.genome[1] for gty in aGty ]
	p = zeros( graph.Nv )
	for g in aG
		p[g] += 1
	end
	return p/length(aGty)
end

function popStatisticsVaryingTabularI(aGty::Vector{<:AbstractGenotype},G::AbstractGraph)
	aG = [ gty.genome for gty in aGty ]

	p = zeros( graph.Nv )
	b = [ zeros(4) for i in 1:graph.Nv ]

	for g in aG
		p[g[1]] += 1
		if g[2] > 0
			b[g[1]][g[2]] += 1
		end
	end
	return p/length(aGty), b ./ p
end

export popStatistics
