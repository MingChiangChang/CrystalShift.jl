# Do breadth-first-search
# Recursive?

mutable struct Node{PM<:PhaseModel, CN<:AbstractVector,
	        R<:AbstractVector, K<:AbstractVector, I<:Real}
	phase_model::PM
	child_node::CN

	id::Int
	recon::R
	residual::K
	inner::I
	is_optimized::Bool
end

# TODO: allow PM() to create empty phasemodel object
Node() = Node(PhaseModel(), Node[], 1, Float64[], Float64[], 0., false) # empty root this will not do amorphous identification
Node(bg::BackgroundModel) = Node(PhaseModel(bg), Node[], 1, Float64[], Float64[], 0., false) # Root for doing amorphous identification
# Node(CP::CrystalPhase, id::Int) = Node(PhaseModel(CP), Node[], id, Float64[], Float64[], 0., false)
# Node(CPs::AbstractVector{<:CrystalPhase}, id::Int) = Node(PhaseModel(CPs), Node[], id, Float64[], Float64[], 0., false)
# Node(PM::PhaseModel, id::Int) = Node(PM, Node[], id, Float64[], Float64[], 0., false)

get_id(CP::CrystalPhase) = CP.id


function Node(PM::PhaseModel, id::Int)
	Node(PhaseModel(CrystalPhase.(PM.CPs, (PM.CPs[1].σ, ), (PM.CPs[1].profile, )),
								  PM.wildcard, PM.background),
								  Node[], id, Float64[], Float64[], 0., false)
end


function Node(CP::CrystalPhase, id::Int)
	Node(PhaseModel(CP), Node[], id, Float64[], Float64[], 0., false)
end

function Node(CPs::AbstractVector{<:CrystalPhase}, id::Int)
	Node(PhaseModel(CPs), Node[], id, Float64[], Float64[], 0., false)
end

function Node(node::Node, phases::AbstractVector{<:CrystalPhase},
			x::AbstractVector, y::AbstractVector, isOptimized::Bool = true)
	check_same_phase(node, phases) || error("Phases must be the same as in the node")
	recon = phases.(x)
	Node(PhaseModel(phases), node.child_node, node.id, 
	     recon, y.-recon, cos_angle(recon, y), isOptimized)
end

function Node(node::Node, PM::PhaseModel,
              x::AbstractVector, y::AbstractVector, isOptimized::Bool = true)
	check_same_phase(node, PM) || error("Phases must be the same as in the node")
	recon = PM(x)
	Node(PM, node.child_node, node.id,
	     recon, y.-recon, cos_angle(recon, y), isOptimized)
end

function Node(PM::PhaseModel, child_nodes::AbstractVector, id::Int64,
		x::AbstractVector, y::AbstractVector)
	recon = PM(x)
	Node(PM, child_nodes, id, recon, y.-recon, cos_angle(y, recon), false)
end
# TODO: Include background into PhaseModel
# function Node(CPs::AbstractVector{<:CrystalPhase},
# 	          child_nodes::AbstractVector,
# 			  x::AbstractVector, y::AbstractVector, id::Int) 
# 	recon = CPs.(x)
#     Node(CPs, child_nodes, id, recon, y.-recon, cos_angle(y, recon), false)
# end


Base.size(n::Nothing) = 0
Base.size(node::Node) = size(node.phase_model.CPs)
Base.getindex(n::Node, i::Int) = Base.getindex(n.phase_model, i)
Base.getindex(n::Node, I::Vector{Int}) = [n[i] for i in I]
get_phase_ids(n::Node) = get_phase_ids(n.phase_model)
get_phase_ids(n::Nothing) = []



function check_same_phase(node::Node, phases::AbstractVector{<:CrystalPhase})
	check_same_phase(node.phase_model, phases)
end

function check_same_phase(node::Node, PM::PhaseModel)
	check_same_phase(node.phase_model, PM)
end

function check_same_phase(PM1::PhaseModel, PM2::PhaseModel)
	check_same_phase(PM1.CPs, PM2.CPs)
end
function check_same_phase(PM::PhaseModel, 
	                      phase_comb::AbstractVector{<:AbstractPhase})
	check_same_phase(PM.CPs, phase_comb)
end

function check_same_phase(phase_comb::AbstractVector{<:AbstractPhase},
	                      PM::PhaseModel)
    check_same_phase(PM.CPs, phase_comb)
end

check_same_phase(::Nothing, ::Nothing) = true

function check_same_phase(phase_comb1::AbstractVector{<:AbstractPhase},
	                      phase_comb2::AbstractVector{<:AbstractPhase})
    # ids_1 = Set([phase_comb1[i].id for i in eachindex(phase_comb1)])
	# ids_2 = Set([phase_comb2[i].id for i in eachindex(phase_comb2)])

	names_1 = Set([phase_comb1[i].name for i in eachindex(phase_comb1)])
	names_2 = Set([phase_comb2[i].name for i in eachindex(phase_comb2)])

	names_1 == names_2 # ids_1 == ids_2 &&
end

function Base.show(io::IO, node::Node)
	println("Node ID: $(node.id)")
	if !isnothing(node.phase_model.CPs)
		println("Phases:")
		for phase in node.phase_model
			println("    $(phase.name)")
		end
	end
	println("Number of child nodes: $(size(node.child_node))")
	println("Inner product: $(node.inner)")
	if node.is_optimized
	    println("Optimized: Yes")
	else
		println("Optimized: No")
	end
	println("")
end

function Base.:(==)(a::Node, b::Node)
    return a.phase_model == b.phase_model
end

function is_child(parent::Node, child::Node)
	return issubset([p.id for p in parent.phase_model],
	               [p.id for p in child.phase_model])
end

function get_child_node_indicies(parent::Node, l_nodes)
    indicies = Int64[]
	for (idx, node) in enumerate(l_nodes)
        if is_child(parent, node)
            push!(indicies, idx)
		end
	end
	indicies
end

function is_immidiate_child(parent::Node, child::Node)
    return (issubset([p.id for p in parent.phase_model],
	         [p.id for p in child.phase_model]) &&
			  (get_level(parent)-get_level(child) == -1))
end

function add_child!(parent::Node, child::Node)
    push!(parent.child_node, child)
end

function remove_child!(parent::Node, child::Node)
	for (idx, cn) in enumerate(parent.child_node)
        if cn == child
	        deleteat!(parent.child_node, idx)
		end
	end
end

get_level(node::Node) = isnothing(node.phase_model.CPs) ? 0 : size(node.phase_model, 1)
get_inner(nodes::AbstractVector{<:Node}) = [node.inner for node in nodes]
get_child_ids(node::Node) = [node.child_node[i].id for i in eachindex(node.child_node)]
get_ids(nodes::AbstractVector{<:Node}) = [node.id for node in nodes]


# O(n) for now, can improve to O(1)
function get_nodes_at_level(nodes::AbstractVector{<:Node}, level::Int)
	idx = [i for i in eachindex(nodes) if get_level(nodes[i])==level]
    return @view nodes[idx]
end

function get_node_with_id(nodes::AbstractVector, id::Int)
    for i in eachindex(nodes)
		if !isempty(get_phase_ids(nodes[i])) && get_phase_ids(nodes[i])[1] == id
			return @view nodes[i:i]
		end
	end
end


function get_node_with_id(nodes::AbstractVector, ids::AbstractVector{<:Int})
	indices = Vector{Int}()
    for i in eachindex(nodes)
		if !isempty(get_phase_ids(nodes[i])) && get_phase_ids(nodes[i])[1] in ids
			push!(indices, i)
		end
	end

	return @view nodes[indices]
end


function get_node_with_exact_ids(nodes::AbstractVector, ids::AbstractVector)
    for i in eachindex(nodes)
		if get_phase_ids(nodes[i]) == ids
            return i, @view nodes[i]
		end
	end
end

(node::Node)(x::AbstractVector) = node.phase_model(x)
cos_angle(node::Node, x::AbstractVector) = cos_angle(node(x), x)


function fit!(node::Node, x::AbstractVector, y::AbstractVector,
	std_noise::Real, mean::AbstractVector, std::AbstractVector,
	maxiter=32, regularization::Bool=true)
    optimized_phases, residuals = fit_phases(node.phase_model, x, y,
									std_noise, mean, std,
									maxiter=maxiter,
									regularization=regularization)
    node.phase_model = optimized_phases
	return residuals
end

cos_angle(x1::AbstractVector, x2::AbstractVector) = x1'x2/(norm(x1)*norm(x2))

function cos_angle(node1::Node, node2::Node, x::AbstractArray)
    x1, x2 = node1(x), node2(x)
    cos_angle(x1, x2)
end

function sum_recon(nodes::AbstractVector{<:Node})
    s = zeros(size(nodes[1].recon, 1))
    for node in nodes
        s += node.recon
	end
	s./maximum(s) # Remove normalizzation
end

function residual(node::Node, θ::AbstractVector,
	              x::AbstractVector, y::AbstractVector)
	residual = copy(y)

	for phase in node.phase_model
	    full_θ = get_eight_params(phase, θ)
	    residual .-= CrystalPhase(phase, θ).(x)
		plot!(x, residual)
	end

	return norm(residual)
end

# function residual(node::Node, x::AbstractVector, y::AbstractVector)
# 	θ = Float64[]
#     for phase in node.phase_model
# 		θ = [θ; get_free_params(phase)]
# 	end
#     residual(node, θ, x, y)
# end
