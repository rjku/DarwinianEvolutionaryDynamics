
module mGraphs

using SparseArrays

abstract type Graph end
	# Graphs are characterized by:
	# 	nv, #vertices; ne, #edges; V = [ Vi, Vj ], adjacent vertices.

abstract type Digraph <: Graph end
	# Digraphs are characterized by:
	# 	nv, #vertices; ne, #edges; V = [ T, O ], where T, vector for the target vertices, and O, vector for the origin vertices.

export Graph, Digraph


# ----------------
# HasWeights Trait

abstract type HasWeights end
struct HasNoWeights <: HasWeights end
struct HasVertexWeights <: HasWeights end
struct HasEdgeWeights <: HasWeights end

# default behavior
HasWeights(::Type) = HasNoWeights()

# neighbors(G::Graph, v::Integer) = vcat( G.V[1][findall( e -> e == v, G.V[2] )], G.V[2][findall( e -> e == v, G.V[1] )] )
#
# function neighbors(G::EdgeWeightedGraph, v::Integer)
# 	ej = findall( e -> e == v, G.V[2] )
# 	ei = findall( e -> e == v, G.V[1] )
#
# 	return vcat( G.V[1][ej], G.V[2][ei] ), G.W[ vcat(ej,ei) ]
# end

# =============
# SIMPLE GRAPHS
# =============

abstract type Lattice <: Graph end

abstract type SquareLattice <: Lattice end

struct FiniteSquareLattice <: SquareLattice
	Nv::Int32
	Ne::Int32
	V::Vector{Vector{Int32}}

	L::Int32
end

function FiniteSquareLattice(L::Integer)
	# vertices and edges as arranged as
	#
	#	3 --- 6 --- 9		.  9  .  12 .
	# 	|	  |		|		2	  4		6
	#	2 --- 5 --- 8		.  8  .  11 .
	# 	|	  |		|		1	  3		5
	# 	1 --- 4 --- 7		.  7  .  10 .

	Nv = Int32(L^2)
	L2mL = Int32(Nv-L)
	Ne = 2L2mL

	Vi = Vector{Int32}(undef, Ne)
	Vj = Vector{Int32}(undef, Ne)

	for ie in 1:L2mL
		# "vertical" edges
		Vi[ie] = Int32( ie + (ie-1)÷(L-1) );
		Vj[ie] = Int32( Vi[ie] + 1 );

		# "horizontal" edges
		Vi[ie+L2mL] = Int32( ie )
		Vj[ie+L2mL] = Int32( ie + L )
	end

	FiniteSquareLattice( Nv, Ne, [Vi, Vj], Int32(L) )
end

export Lattice, SquareLattice, FiniteSquareLattice


# ===============
# WEIGHTED GRAPHS
# ===============

struct VertexWeightedSquareLattice{Tw,Taw<:AbstractArray{Tw}} <: SquareLattice
	Nv::Int32
	Ne::Int32
	V::Vector{Vector{Int32}}
	Wv::Taw

	L::Int32
end

HasWeights(::Type{<:VertexWeightedSquareLattice}) = HasVertexWeights()

function VertexWeightedSquareLattice(L::Integer, T::DataType)
	SquareLattice = FiniteSquareLattice(L)
	VertexWeightedSquareLattice( SquareLattice.Nv, SquareLattice.Ne, SquareLattice.V, Vector{T}(undef, SquareLattice.Ne), Int32(L) )
end

function VertexWeightedSquareLattice(L::Integer, Wv::AbstractArray)
	SquareLattice = FiniteSquareLattice(L)
	VertexWeightedSquareLattice( SquareLattice.Nv, SquareLattice.Ne, SquareLattice.V, Wv, Int32(L) )
end

struct EdgeWeightedSquareLattice{Tw,Taw<:AbstractArray{Tw}} <: SquareLattice
	Nv::Int32
	Ne::Int32
	V::Vector{Vector{Int32}}
	We::Taw

	L::Int32
end

HasWeights(::Type{<:EdgeWeightedSquareLattice}) = HasEdgeWeights()

function EdgeWeightedSquareLattice(L::Integer, T::DataType)
	SquareLattice = FiniteSquareLattice(L)
	EdgeWeightedSquareLattice( SquareLattice.Nv, SquareLattice.Ne, SquareLattice.V, Vector{T}(undef, SquareLattice.Ne), Int32(L) )
end

function EdgeWeightedSquareLattice(L::Integer, We::AbstractArray)
	SquareLattice = FiniteSquareLattice(L)
	EdgeWeightedSquareLattice( SquareLattice.Nv, SquareLattice.Ne, SquareLattice.V, We, Int32(L) )
end

export EdgeWeightedSquareLattice, EdgeWeightedSquareLattice

# =======
# METHODS
# =======

neighbors(G::T, v::Integer) where T<:Graph = _neighbors(HasWeights(T),G,v)

_neighbors(::HasNoWeights, G::Graph, v::Integer) = vcat( G.V[1][findall( e -> e == v, G.V[2] )], G.V[2][findall( e -> e == v, G.V[1] )] )

function _neighbors(::HasEdgeWeights, G::Graph, v::Integer)
	ej = findall( e -> e == v, G.V[2] )
	ei = findall( e -> e == v, G.V[1] )

	return vcat( G.V[1][ej], G.V[2][ei] ), G.We[ vcat(ej,ei) ]
end


transitionMatrix(G::T) where T<:Graph = _transitionMatrix(HasWeights(T),G)
_transitionMatrix(::HasNoWeights,G::Graph) = throw(error("Graph $(G) has no weights"))

function _transitionMatrix(::HasEdgeWeights,G::Graph)
	Vi = vcat( G.V[1], G.V[2] )
	Vj = vcat( G.V[2], G.V[1] )
	We = vcat( G.We, G.We )

	T = sparse(Vi, Vj, We, G.Nv, G.Nv)
	for i in 1:G.Nv
		T[i,i] = -sum(T[:,i])
	end

	return T
end

matrixForm(G::T) where T<:SquareLattice = _matrixForm(HasWeights(T),G)
_matrixForm(::HasNoWeights,G::SquareLattice) = throw(error("Graph $(G) has no weights"))

_matrixForm(::HasVertexWeights,G::SquareLattice) = reshape(G.Wv,(G.L,G.L))

export neighbors, transitionMatrix

end # module
