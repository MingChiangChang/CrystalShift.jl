struct MPTree{NS<:AbstractVector{<:Node},
                AS<:AbstractSet,
                CP<:AbstractVector{<:AbstractPhase},
                AR<:AbstractVector,
                G,
                T,
                H } <: AbstractTree

    nodes::NS
    phase_combinations::AS # Keeping track of phase combinations
    phases::CP
    x::AR
    id_lookup_tabel::G
    simulated_patterns::T
    mp_res::H
end

function MPTree(CPs::AbstractVector{<:AbstractPhase}, x::AbstractVector,
    n_patterns::Int64=10, max_strain::Float64=0.05, w_min::Float64=0.1, w_max::Float64=0.3)

    id_lookup_table = Dict{Int64, Int64}()
    simulated_patterns = zeros(length(CPs), n_patterns, length(x))

    for i in eachindex(CPs)
        id_lookup_table[CPs[i].id] = i
        simulated_patterns[i, 1, :] .= evaluate!(simulated_patterns[i, 1, :], CPs[i:i], x)
        for j in 2:n_patterns
            simulated_patterns[i, j, :] = _simulate_pattern(CPs[i], x, max_strain, w_min, w_max)
        end
    end

    for i in axes(simulated_patterns, 1)
        for j in axes(simulated_patterns, 2)
            simulated_patterns[i, j, :] ./= norm(simulated_patterns[i, j, :])
        end
    end

    mp_res = zeros((length(CPs), n_patterns))

    MPTree(Node[], Set(), CPs, x, id_lookup_table, simulated_patterns, mp_res)
end


function search!(mptree::MPTree, x::AbstractVector, y::AbstractVector, y_uncer::AbstractVector,
                setting::MPTreeSearchSettings)

    depth, k = setting.depth, setting.k
    result = Vector{Vector{<:Node}}(undef, depth+1)

    if setting.amorphous
        bg = BackgroundModel(x, EQ(), 8., 10., rank_tol=1e-3)
        push!(mptree.nodes, Node(bg))
    else
        push!(mptree.nodes, Node())
    end

    # expand!(LT, LT.nodes[1])

    for level in 1:depth+1
        nodes = get_nodes_at_level(mptree, level-1)
        level_result = Vector{Node}(undef, length(nodes))

        @threads for i in eachindex(nodes)
            if !isnothing(nodes[i].phase_model.background) || !isempty(nodes[i].phase_model.CPs)
                pm = optimize!(nodes[i].phase_model, x, y, y_uncer, setting.opt_stn)
                if pm isa Tuple # The uncertainty flag returns two result, awful solution for now
                    pm = pm[1]
                end
                nodes[i] = Node(nodes[i], pm, x, y, true)
            end
            level_result[i] = nodes[i]
        end

        result[level] = level_result

        if level == 1
            top_k_nodes = result[1]
        else
            top_k_nodes = get_top_k_nodes(result, k, level)
        end

        if level != depth+1
            for node in top_k_nodes
                mp_expand!(mptree, node, x, setting.mp_top_k, setting.background, setting.background_length)
            end
        end
    end
    return result
end

function mp_expand!(mptree::MPTree, node::Node, x::AbstractVector, k::Integer, background::Bool, background_length::Float64)
    starting_id = get_max_id(mptree) + 1
    mp_expand!(mptree, node, x, k, background, background_length, starting_id)
end

function mp_expand!(mpt::MPTree, node::Node, x::AbstractVector, k::Integer, background::Bool, background_length::Float64, starting_id::Int)
    child_nodes = create_child_nodes_mp(mpt, node, x, k, background, background_length, starting_id)
    attach_child_nodes!(node, child_nodes)
    push!(mpt.nodes, child_nodes...)
    # println(Set.(get_phase_ids.(child_nodes))...) # TODO: This fails when k is larger than the number of phase
    if !isempty(child_nodes)
        push!(mpt.phase_combinations, Set.(get_phase_ids.(child_nodes))...)
    end
    return child_nodes
end

function create_child_nodes_mp(mpt::MPTree, node::Node, x::AbstractVector, k::Integer,
                                background::Bool, background_length::Float64, starting_id::Int)

    new_nodes = Node[]
    if isnothing(node.phase_model.CPs) || isempty(node.phase_model.CPs)# Root node case
        for i in eachindex(mpt.phases)
            push!(new_nodes, Node(add_phase(node.phase_model, mpt.phases[i], x, background, background_length),
                                 starting_id+i-1))
        end
        return new_nodes
    end

    top_indices = get_top_max_inner!(mpt.mp_res, mpt.simulated_patterns, node.residual)[1:k]
    for i in eachindex(top_indices)#eachindex(LT.phases)
        if is_allowed_new_phase(mpt, node, mpt.phases[top_indices[i]])
            push!(new_nodes, Node(add_phase(node.phase_model, mpt.phases[top_indices[i]], x, background, background_length),
                                 starting_id+i-1))
        end
    end
    return new_nodes


end

function get_top_max_inner!(mp_res, simulated_patterns, residual)
    @einsum mp_res[i, j] = simulated_patterns[i, j, k] * residual[k]
    max_inners = maximum(mp_res, dims=2)
    max_inners = dropdims(max_inners, dims=2)
    sortperm(max_inners, rev=true)
end


# function expand!(LT::Lazytree, node::Node, x::AbstractVector, background::Bool, l::Real)
#     starting_id = get_max_id(LT) + 1
#     expand!(LT, node, x, background, l, starting_id)
# end
# 
# 
# function expand!(LT::Lazytree, node::Node, x::AbstractVector, background::Bool, l::Real, starting_id::Int)
#     child_nodes = create_child_nodes(LT, node, x, background, l, starting_id)
#     attach_child_nodes!(node, child_nodes)
#     push!(LT.nodes, child_nodes...)
#     # println(Set.(get_phase_ids.(child_nodes))...) # TODO: This fails when k is larger than the number of phase
#     if !isempty(child_nodes)
#         push!(LT.phase_combinations, Set.(get_phase_ids.(child_nodes))...)
#     end
#     return child_nodes
# end
