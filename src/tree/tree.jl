# TODO: Parameter/Criteria for stopping
# 1. Low residual (good stop)
# 2. To many extra peaks (bad stop/killed)
#    Condition: norm(max.(fitted-data)./data) ?
#    Or is there anything better than 2 norm
using Combinatorics
using LinearAlgebra
using Base.Threads
using CrystalShift: evaluate!

abstract type AbstractTree end

struct Tree{T, CP<:AbstractVector{T}, DP<:Integer} <: AbstractTree
    nodes::CP
    depth::DP # Store for convenience
end

function Tree(phases::AbstractVector{<:CrystalPhase}, depth::Int)
    # Construct tree with certain depth
    nodes = Node[]
	root = Node()
	push!(nodes, root)

	id = 2

    for d in 1:depth
		phase_combs = combinations(phases, d)
		nodes_at_level = get_nodes_at_level(nodes, d-1)
		for phases in phase_combs
			new_node = Node(phases, id)
			id += 1
			for old_node in nodes_at_level
			    if is_immidiate_child(old_node, new_node)
				    add_child!(old_node, new_node)
				    push!(nodes, new_node)
					break
                end
		    end
		end
    end
	Tree(nodes, depth)
end

Base.view(tree::AbstractTree, i) = Base.view(tree.nodes, i)
Base.size(t::AbstractTree) = size(t.nodes)
Base.size(t::AbstractTree, dim::Int) = size(t.nodes, dim)
Base.getindex(t::AbstractTree, i::Int) = Base.getindex(t.nodes, i)
Base.getindex(t::AbstractTree, I::Vector{Int}) = [t[i] for i in I]

get_nodes_at_level(tree::AbstractTree, level::Int) = get_nodes_at_level(tree.nodes, level)
get_node_with_id(tree::AbstractTree, id::Int) = get_node_with_id(tree.nodes, id)
get_node_with_id(tree::AbstractTree, ids::AbstractVector{<:Int}) = get_node_with_id(tree.nodes, ids)

function bft(t::AbstractTree)
    # Breadth-first traversal, return an array of
	# Node with the b-f order
	traversal = Int[]
	for i in 1:t.depth
	    for j in eachindex(t.nodes)
            if get_level(t.nodes[j]) == i
				push!(traversal, j)
			end
	    end
	end
	@view t[traversal]
end

bft(t::AbstractTree, level::Int) = get_nodes_at_level(t.nodes, level)

function dft(t::Tree)
    # Depth-first traversal, return an array of Node with
	# the D-F order
end

# subtree
function remove_subtree!(nodes::AbstractVector{<:Node}, root_of_subtree::Node)
    # Given a vector of node and a node, remove
	# all the node that are child of the node
	# TODO: Should remove there relationship as well??
	# TODO: Will GC take care?
	to_be_removed = Int[]
    for (idx, node) in enumerate(nodes)
		if is_child(root_of_subtree, node) && root_of_subtree != node
			push!(to_be_removed, idx)
		end
	end
	deleteat!(nodes, to_be_removed)
end

function get_optimized_nodes_at_level(tree::Tree, level::Int)
    get_optimized_nodes_at_level(tree.nodes, level)
end

function get_optimized_nodes_at_level(nodes::AbstractVector{<:Node}, level::Int)
    idx = [i for i in eachindex(nodes) if get_level(nodes[i])==level && nodes[i].is_optimized]
    return @view nodes[idx] 
end