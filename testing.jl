
using Base.Threads, Random, BenchmarkTools, LinearAlgebra, SparseArrays, ForwardDiff, Distributions, MATLAB, DelimitedFiles #, PyPlot

abstract type atEvotype end
abstract type atMetaGenotype end
abstract type atGenotype end
abstract type atPhenotype end
abstract type atGtyMutEty <: atEvotype end
abstract type atPntMutEty <: atEvotype end
abstract type atIsingMetaGty <: atMetaGenotype end
abstract type atChannelMetaGty <: atMetaGenotype end
abstract type atSystemGty{Tmgty} <: atGenotype end		# genotypes endowed with a pointer to a atMetaGenotype: pMetaGty

# type: disordered channel metagenotype <: atChannelMetaGty
struct tDisChnMGty <: atChannelMetaGty
	L::Int32					# system size
	L2::Int32;	L2mL::Int32
	dG::Int32;	Nvbl::Int32		# genotype length and viable transitions
	p::Float64 					# viability probability
	V::Vector{Vector{Int32}}	# [ t(e), o(e) ] vector
	kout::Float64				# output rate constant
end

# constructor: disordered channel
function tDisChnMGty(L::Int32, p::Float64, kout::Float64)
	L2mL = L*(L-1)

	viability = rand(Bernoulli(p),2L2mL);	Nvbl = sum(viability)
	V = [ Vector{Int32}(undef, 2Nvbl), Vector{Int32}(undef, 2Nvbl) ]

	g = 0
	for (i,v) in enumerate(viability)
		if i <= L2mL && Bool(v)
			g+=1
			V[1][g] = i + (i-1)÷(L-1) + 1;	V[2][g] = i + (i-1)÷(L-1);
			V[1][g+Nvbl] = V[2][g];			V[2][g+Nvbl] = V[1][g];
		elseif i > L2mL && Bool(v)
			g+=1
			V[1][g] = i - L2mL + L;		V[2][g] = i - L2mL;
			V[1][g+Nvbl] = V[2][g];		V[2][g+Nvbl] = V[1][g];
		end
	end

	tDisChnMGty(L, L^2, L2mL, 2Nvbl, Nvbl, p, V, kout)
end

# type: genotype with continuous bounded genotypic variables
struct tCntGty{Tmgty<:atMetaGenotype,Tpmgty<:Vector{Tmgty},Tag<:AbstractArray} <: atSystemGty{Tmgty}
	pMetaGty::Tpmgty
	G::Tag
	gbounds::Tag
	pF::Vector{Float64}
end

# constructor: undefined fitness
tCntGty(pMGty::Vector{Tmgty},G::Tag,gbounds::Tag) where {Tmgty<:atMetaGenotype,Tag<:AbstractArray} =
	tCntGty{Tmgty,Vector{Tmgty},Tag}(pMGty,G,gbounds,[0.0])

struct tCompEnv{T<:Vector{<:AbstractVector}}
	IOidl::T
	selFactor::Float64
end

# function: constructor of stochastic matrix for channel metagenotype
# The last column is kept null as it is reserved for input transitions
function getW(gty::atSystemGty{<:atChannelMetaGty})
	W = sparse(gty.pMetaGty[1].V..., gty.G, gty.pMetaGty[1].L2+1, gty.pMetaGty[1].L2+1)
	W[end,gty.pMetaGty[1].L2mL+1:gty.pMetaGty[1].L2] .= gty.pMetaGty[1].kout
	for i in 1:gty.pMetaGty[1].L2
		W[i,i] = -sum(W[:,i])
	end
	return W
end

function responseFD(gty::atSystemGty{<:atChannelMetaGty},W::AbstractMatrix,Input::Vector{<:Real} )
	length(Input) == gty.pMetaGty[1].L || throw( "inconsistent dimensions between input currents and system" )
	W[1:gty.pMetaGty[1].L,end] = Input
	W[end,end] = -sum(W[1:end-1,end])

	Λ(q::Vector) = det( Array(W .* ( 1.0 .+ sparse( fill(gty.pMetaGty[1].L2+1,gty.pMetaGty[1].L), gty.pMetaGty[1].L2mL+1:gty.pMetaGty[1].L2, map(e -> exp(e) - 1.0, q[1:end-1]),gty.pMetaGty[1].L2+1,gty.pMetaGty[1].L2+1 ) ) ) - q[end] .* I )
	dΛ(q::Vector) = ForwardDiff.gradient(Λ,q)

	q0 = zeros(Float64,gty.pMetaGty[1].L+1);	y = dΛ(q0);

	if y[end] != 0.0
		return [ -y[i]/y[end] for i in 1:gty.pMetaGty[1].L ]
	else
		q0[end] = 10^-14;	y = dΛ(q0);
		if y[end] != 0.0
			return [ -y[i]/y[end] for i in 1:gty.pMetaGty[1].L ]
		else
			throw( "a₁(0) = 0. Irreducible stochastic matrix ... or numerical error" )
		end
	end
end

# function: constructor of normalized stochastic matrix for disordered channel metagenotype
# The last column is kept null as it is reserved for input transitions. The last row encodes the normalization.
function getWnrmd(gty::atSystemGty{<:atChannelMetaGty})
	W = sparse(gty.pMetaGty[1].V..., gty.G, gty.pMetaGty[1].L2+1, gty.pMetaGty[1].L2+1)
	W[end,gty.pMetaGty[1].L2mL+1:gty.pMetaGty[1].L2] .= gty.pMetaGty[1].kout
	for i in 1:gty.pMetaGty[1].L2
		W[i,i] = -sum(W[:,i])
	end
	W[end,:] = ones(Float64,gty.pMetaGty[1].L2+1)
	return W
end

function response(gty::atSystemGty{<:atChannelMetaGty},W::AbstractMatrix,Input::Vector{<:Real})
	length(Input) == gty.pMetaGty[1].L || throw( "inconsistent dimensions between input currents and system" )
	W[1:gty.pMetaGty[1].L,end] = Input
	b = zeros(Float64,gty.pMetaGty[1].L2+1);	b[end] = 1.0

	return gty.pMetaGty[1].kout .* ( ( W \ b )[gty.pMetaGty[1].L2mL+1:gty.pMetaGty[1].L2] )
end

function fitness(istEnv::tCompEnv{<:Vector{<:Vector{<:Vector{<:Real}}}},gty::atSystemGty{<:atChannelMetaGty})::Float64
	W = getWnrmd(gty)

	d2::Float64 = 0.0
	for io in istEnv.IOidl
		d2 += sum( (response(gty, W, io[1]) - io[2]).^2 )
	end
	return exp(-d2/istEnv.selFactor)
end

function chopArray!(a::AbstractArray, δ::Float64=10^-12)
	for (i, e) in enumerate(a)
		if abs(e) < δ
			a[i] = 0.0
		end
	end
end

struct tFluxPattern
	Nq::Int32
	V::Vector{Vector{Int32}}	# [ t(e), o(e) ] vector
	jave::Vector{Vector{Float64}}
	jcov::Vector{Array{Float64,2}}
	jcor::Array{Float64,2}
end

# function tFluxPattern(env::atCompEnv,gty::atSystemGty{<:atChannelMetaGty})
function tFluxPattern(env,gty::atSystemGty{<:atChannelMetaGty})
	Nq = gty.pMetaGty[1].dG + 2gty.pMetaGty[1].L

	tFluxPattern( Nq,
		[ Vector{Int32}(undef, Nq), Vector{Int32}(undef, Nq) ], [ Vector{Float64}(undef, Nq) for i in eachindex(env.IOidl) ],
		[ Array{Float64}(undef, Nq, Nq) for i in eachindex(env.IOidl) ], Array{Float64}(undef, Nq, Nq)
	)
end

function getCrntStat!(env::tCompEnv,gty::atSystemGty{<:atChannelMetaGty},fluxPtrn::tFluxPattern)
	length(fluxPtrn.jave) == length(env.IOidl) || throw( "inconsistent dimensions between _fluxPtrn.j_ currents and _env.IOidl_" )
	aCrntCorFisherz = Vector{Array{Float64,2}}(undef,length(env.IOidl))

	W = getW(gty)
	Nq = gty.pMetaGty[1].dG + 2gty.pMetaGty[1].L
	q0 = zeros(Float64, Nq + 1)
	for (i, io) in enumerate(env.IOidl)
		W[1:gty.pMetaGty[1].L,end] = io[1]
		W[end,end] = -sum(W[1:end-1,end])

		Λ(q::Vector) = det( Array(W .* ( 1.0 .+
			sparse( gty.pMetaGty[1].V[1][1:gty.pMetaGty[1].Nvbl], gty.pMetaGty[1].V[2][1:gty.pMetaGty[1].Nvbl], map(e -> exp(e) - 1.0, q[1:gty.pMetaGty[1].Nvbl]), gty.pMetaGty[1].L2+1, gty.pMetaGty[1].L2+1 ) .+
			sparse( gty.pMetaGty[1].V[2][1:gty.pMetaGty[1].Nvbl], gty.pMetaGty[1].V[1][1:gty.pMetaGty[1].Nvbl], map(e -> exp(e) - 1.0, q[gty.pMetaGty[1].Nvbl+1:gty.pMetaGty[1].dG]), gty.pMetaGty[1].L2+1, gty.pMetaGty[1].L2+1) .+
			sparse( fill(gty.pMetaGty[1].L2+1,gty.pMetaGty[1].L), collect(gty.pMetaGty[1].L2mL+1:gty.pMetaGty[1].L2), map(e -> exp(e) - 1.0, q[gty.pMetaGty[1].dG+1:gty.pMetaGty[1].dG+gty.pMetaGty[1].L]), gty.pMetaGty[1].L2+1, gty.pMetaGty[1].L2+1) .+
			sparse( collect(1:gty.pMetaGty[1].L), fill(gty.pMetaGty[1].L2+1,gty.pMetaGty[1].L), map(e -> exp(e) - 1.0, q[gty.pMetaGty[1].dG+gty.pMetaGty[1].L+1:Nq]), gty.pMetaGty[1].L2+1, gty.pMetaGty[1].L2+1) ) ) -
			q[end] .* I )
		dΛ(q::Vector) = ForwardDiff.gradient(Λ,q)
		ddΛ(q::Vector) = ForwardDiff.hessian(Λ,q)

		jave = dΛ(q0)
		if jave[end] == 0.0
			q0[end] = 10^-14
			jave = dΛ(q0)
		end
		if jave[end] == 0.0
			throw( "a₁(0) = 0. Irreducible stochastic matrix ... or numerical error" )
		end

		jcov = ddΛ(q0)
		chopArray!(jcov)

		fluxPtrn.jave[i] .= [ -jave[j]/jave[end] for j in 1:Nq ]
		fluxPtrn.jcov[i] .= [ -( jcov[ji,jj] + jcov[ji,end]*fluxPtrn.jave[i][jj] + jcov[jj,end]*fluxPtrn.jave[i][ji] + jcov[end,end]*fluxPtrn.jave[i][ji]*fluxPtrn.jave[i][jj] )/jave[end] for ji in 1:Nq, jj in 1:Nq ]
		aCrntCorFisherz[i] = broadcast( r -> log( (1 + r)/(1 - r) )/2, [ fluxPtrn.jcov[i][ji,jj] / sqrt( fluxPtrn.jcov[i][ji,ji] * fluxPtrn.jcov[i][jj,jj] ) for ji in 1:Nq, jj in 1:Nq ] )

		q0[end] = 0.0
	end

	fluxPtrn.jcor .= [ tanh( mean( [ aCrntCorFisherz[i][ji,jj] for i in eachindex(env.IOidl) ] ) ) for ji in 1:Nq, jj in 1:Nq ]

	fluxPtrn.V[1] .= vcat( gty.pMetaGty[1].V[1], fill(gty.pMetaGty[1].L2+1,gty.pMetaGty[1].L), collect(1:gty.pMetaGty[1].L) )
	fluxPtrn.V[2] .= vcat( gty.pMetaGty[1].V[2], collect(gty.pMetaGty[1].L2mL+1:gty.pMetaGty[1].L2), fill(gty.pMetaGty[1].L2+2,gty.pMetaGty[1].L) )

	return 0
end

L = Int32(2); p = 1.0; Nio = 2
testChannel = tDisChnMGty(L, p, 1.0)
testGty = tCntGty([testChannel],[ rand() for i in 1:testChannel.dG ],[0.0,1.0] )
testEnv = tCompEnv([ [ Float64[rand(Bool) + 0.001 for i in 1:L], Float64[rand(Bool) for i in 1:L] ] for i in 1:Nio ],1.0)
testFluxPtrn = tFluxPattern(testEnv, testGty)

getCrntStat!(testEnv,testGty,testFluxPtrn)
display( testFluxPtrn.jcov[2] )
display( testFluxPtrn.jcor )

responseFD(testGty,getW(testGty),testEnv.IOidl[1][1]) ≈ response(testGty,getWnrmd(testGty),testEnv.IOidl[1][1])
display(responseFD(testGty,getW(testGty),testEnv.IOidl[1][1]))
display(response(testGty,getWnrmd(testGty),testEnv.IOidl[1][1]))
println(fitness(testEnv,testGty))

display(filter( e -> e < 0.0, testFluxPtrn.jcov[1] ))
display(filter( e -> e < 0.0, testFluxPtrn.jcov[1] ))
display(filter( e -> e < 0.0, diag(testFluxPtrn.jcov[1]) ))
display(filter( e -> e < 0.0, diag(testFluxPtrn.jcov[2]) ))

# using PyPlot
# imshow(testFluxPtrn.jcor,cmap="viridis"); colorbar();

# LWidths = 5.0 .* map(e->(e),fluxPtrn.j[1]) ./ ( maximum(fluxPtrn.j[1]) ) .+ 10^(-8)
# LWidths = 10*rand(length(crntAve.J[1]))

# mat"""
# G = digraph($(fluxPtrn.V[2]),$(fluxPtrn.V[1]),$(fluxPtrn.j[1]))
# G.Edges.EdgeColor = $(LWidths)
# plot(G,'Layout','force','LineWidth',$(LWidths))
# """

# testing graphs

# using LightGraphs, GraphPlot, MATLAB, MetaGraphs
#
# G₁ = Graph(3) # graph with 3 vertices
#
# # make a triangle
# add_edge!(G₁, 1, 2)
# add_edge!(G₁, 1, 3)
# add_edge!(G₁, 2, 3)
#
# gplot(G₁, nodelabel=1:3)
#
# s = [1, 1, 1, 1, 2, 2, 3, 4, 4, 5, 6];
# t = [2 3 4 5 3 6 6 5 7 7 7];
# display(s)
#
# s = [1 1 1 1 2 2 3 4 4 5 6];
#
# mat"""
# weights = [50 10 20 80 90 90 30 20 100 40 60];
# G = graph($s,$t,weights)
#
# LWidths = 5*G.Edges.Weight/max(G.Edges.Weight);
# plot(G,'LineWidth',LWidths)
# """
#
# x = range(-10.0, stop=10.0, length=500)
# mat"plot($x, sin($x))"  # evaluate a MATLAB function
#
# y = range(2.0, stop=3.0, length=500)
# mat"""
#     $u = $x + $y
# 	$v = $x - $y
# """
# @show u v               # u and v are accessible from Julia
