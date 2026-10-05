struct Lazytree{NS<:AbstractVector{<:Node},
                AS<:AbstractSet,
                CP<:AbstractVector{<:AbstractPhase},
                AR<:AbstractVector,
                } <: AbstractTree
    nodes::NS
    phase_combinations::AS # Keeping track of phase combinations
    phases::CP
    x::AR
end

function Lazytree(CPs::AbstractVector{<:AbstractPhase}, x::AbstractVector)
    Lazytree(Node[], Set(), CPs, x)
end


function expand!(LT::Lazytree, node::Node, x::AbstractVector, background::Bool, l::Real)
    starting_id = get_max_id(LT) + 1
    expand!(LT, node, x, background, l, starting_id)
end


function expand!(LT::Lazytree, node::Node, x::AbstractVector, background::Bool, l::Real, starting_id::Int)
    child_nodes = create_child_nodes(LT, node, x, background, l, starting_id)
    attach_child_nodes!(node, child_nodes)
    push!(LT.nodes, child_nodes...)
    # println(Set.(get_phase_ids.(child_nodes))...) # TODO: This fails when k is larger than the number of phase
    if !isempty(child_nodes)
        push!(LT.phase_combinations, Set.(get_phase_ids.(child_nodes))...)
    end
    return child_nodes
end


function add_phase(PM::PhaseModel, phase::AbstractPhase, x::AbstractVector, background::Bool, l::Real)
    if background
        bg = BackgroundModel(x, EQ(), l)
    else
        bg = nothing
    end
    if isnothing(PM.CPs)  || isempty(PM.CPs)
        return PhaseModel([phase], PM.wildcard, bg)
    else
        PhaseModel(vcat(PM.CPs, phase), PM.wildcard, bg)
    end
end


function create_child_nodes(LT::Lazytree, node::Node, x::AbstractVector, background::Bool, l::Real, starting_id::Int)
    new_nodes = Node[]
    for i in eachindex(LT.phases)
        if isnothing(node.phase_model.CPs) || is_allowed_new_phase(LT, node, LT.phases[i])
            push!(new_nodes, Node(add_phase(node.phase_model, LT.phases[i], x, background, l),
                                 starting_id+i-1))
        end
    end
    return new_nodes
end


function is_allowed_new_phase(LT::AbstractTree, node::Node, phase::AbstractPhase)
    current_phases = get_phase_ids(node)
    new_phase = phase.id
    return phase.id ∉ get_phase_ids(node) && Set(vcat(current_phases, new_phase)) ∉ LT.phase_combinations
end


# If having immutable nodes are largely advantageous, try creating all possible nodes
# and then so subarray for child nodes
# OR should we keep having child node recorded in nodes?????
function attach_child_nodes!(node::Node, child_nodes::AbstractVector{<:Node})
    node.child_node = @view child_nodes[1:end] # made this possible by making Node mutable
end


# O(kn) method
function search!(LT::Lazytree, x::AbstractVector, y::AbstractVector, y_uncer::AbstractVector,
                 ts_stn::TreeSearchSettings)
    depth, k = ts_stn.depth, ts_stn.k
    result = Vector{Vector{<:Node}}(undef, depth+1)

    if ts_stn.amorphous
        bg = BackgroundModel(x, EQ(), 8., 10., rank_tol=1e-3)
        push!(LT.nodes, Node(bg))
    elseif !isnothing(ts_stn.default_phase)
        push!(LT.nodes, Node([ts_stn.default_phase], 0))
    else
        push!(LT.nodes, Node())
    end

    # expand!(LT, LT.nodes[1])

    for level in 1:depth+1
        nodes = get_nodes_at_level(LT, level-1)
        level_result = Vector{Node}(undef, length(nodes))

       @threads for i in eachindex(nodes)
            if !isnothing(nodes[i].phase_model.background) || !isempty(nodes[i].phase_model.CPs)
                pm = optimize_node(nodes[i].phase_model, x, y, y_uncer, ts_stn)
                if pm isa Tuple # The uncertainty flag returns two result, awful solution for now
                    pm = pm[1]
                end
                nodes[i] = Node(nodes[i], pm, x, y, true)
            end
            level_result[i] = nodes[i]
        end

        result[level] = level_result

        if level == 1
            top_k = result[1]
        else
            top_k = get_top_k_nodes(result, k, level)
        end

        if level != depth+1
            for i in top_k
                expand!(LT, i, x, ts_stn.background, ts_stn.background_length)
            end
        end
    end
    return result
end

function search!(LT::Lazytree, x::AbstractVector, y::AbstractVector, ts_stn::TreeSearchSettings)
    y_uncer = zero(y)
    search!(LT, x, y, y_uncer, ts_stn)
end

get_top_k_nodes(result::AbstractVector, k::Int, level::Int) = get_top_nodes(result[level], k)
get_top_k_nodes(result::AbstractVector, k::AbstractVector, level::Int) = get_top_nodes(result[level], k[level-1])

# TODO: Remove this and only keep those using Setting objects
function search!(LT::Lazytree, x::AbstractVector, y::AbstractVector, y_uncer::AbstractVector,
                 depth::Integer, k::ScalarOrVecInt, amorphous::Bool, background::Bool, background_length::Real,
                 std_noise::Real, mean::AbstractVector, std::AbstractVector;
                 method::OptimizationMethod = LM, objective::String = "LS",
                 optimize_mode::OptimizationMode = Simple, em_loop_num::Integer =8,
                 maxiter::Integer = 32, regularization::Bool = true, λ::Float64=1., verbose::Bool = false, tol::Real = SEARCH_DEFAULT_TOL)
    opt_stn = OptimizationSettings{eltype(mean)}(std_noise, mean, std, maxiter, regularization, method, objective, optimize_mode, em_loop_num, λ, verbose, tol)
    ts_stn = TreeSearchSettings(depth, k, amorphous, background, background_length, nothing, opt_stn)
    search!(LT, x, y, y_uncer, ts_stn)
end

function search!(LT::Lazytree, x::AbstractVector, y::AbstractVector, 
    depth::Integer, k::ScalarOrVecInt, amorphous::Bool, background::Bool, background_length::Real,
    std_noise::Real, mean::AbstractVector, std::AbstractVector;
    method::OptimizationMethod = LM, objective::String = "LS",
    optimize_mode::OptimizationMode = Simple, em_loop_num::Integer =8,
    maxiter::Integer = 32, regularization::Bool = true, λ::Float64=1., verbose::Bool = false, tol::Real = SEARCH_DEFAULT_TOL)

    y_uncer = zero(y)
    opt_stn = OptimizationSettings{eltype(mean)}(std_noise, mean, std, maxiter, regularization, method, objective, optimize_mode, em_loop_num, λ, verbose, tol)
    ts_stn = TreeSearchSettings(depth, k, amorphous, background, background_length, nothing, opt_stn)
    search!(LT, x, y, y_uncer, ts_stn)
end

# O(k^2 n) method
# expand the node and recussively call search_k2n on the top-k nodes
function search_k2n!(LT::Lazytree, x::AbstractVector, y::AbstractVector, y_uncer::AbstractVector, ts_stn::TreeSearchSettings)

    if ts_stn.amorphous
        bg = BackgroundModel(x, EQ(), 8., 10., rank_tol=1e-3)
        push!(LT.nodes, Node(bg))
    else
        push!(LT.nodes, Node())
    end

    result = Vector{Node}()
    search_k2n!(result, LT, LT.nodes[1], x, y, y_uncer, ts_stn)
    @threads for i in eachindex(result)
        if isassigned(result, i) && !result[i].is_optimized
            pm = optimize_node(result[i].phase_model, x, y, y_uncer, ts_stn)
            result[i] = Node(result[i], pm, x, y, true)
        end
    end
    result[[isassigned(result, i) for i in eachindex(result)]]
end

function search_k2n!(LT::Lazytree, x::AbstractVector, y::AbstractVector, ts_stn::TreeSearchSettings)
    y_uncer = zero(y)
    search_k2n!(LT, x, y, y_uncer, ts_stn)
end

# Doing a mixed version of depth-first search and breadth-first search
# doing best-first strategy at each level than dig futher down
function search_k2n!(result::AbstractVector, LT::Lazytree, node::Node, x::AbstractVector, y::AbstractVector, y_uncer::AbstractVector, ts_stn::TreeSearchSettings)

    if size(node)[1] == ts_stn.depth
        push!(result, node)
        return
    end

    child_nodes = expand!(LT, node, x, ts_stn.background, ts_stn.background_length)
    @threads for i in eachindex(child_nodes)
        pm = optimize_node(child_nodes[i].phase_model, x, y, y_uncer, ts_stn)
        child_nodes[i] = Node(child_nodes[i], pm, x, y, true)
    end

    append!(result, child_nodes)

    ts_stn.k isa AbstractVector && error("Vector k is not supported for search_k2n!")
    top_k = get_top_nodes(child_nodes, ts_stn.k)
    for j in eachindex(top_k)
        search_k2n!(result, LT, top_k[j], x, y, y_uncer, ts_stn)
    end
end


################## Helper functions #################
# Optimize a node's phase model with `optimize!`, or with `full_optimize!` (peak heights
# refit too) when `ts_stn.full_opt_stn` is set. Nodes without phases (e.g. the amorphous
# root) have no peak heights, so they always use `optimize!`.
# Note: `full_optimize!` ignores `y_uncer` and `opt_stn.em_loop_num`.
function optimize_node(pm::PhaseModel, x::AbstractVector, y::AbstractVector, y_uncer::AbstractVector,
                       ts_stn::TreeSearchSettings)
    f_stn = ts_stn.full_opt_stn
    if isnothing(f_stn) || isnothing(pm.CPs) || isempty(pm.CPs)
        return optimize!(pm, x, y, y_uncer, ts_stn.opt_stn)
    end
    opt_stn = ts_stn.opt_stn
    full_optimize!(pm, x, y, opt_stn.priors.std_noise, opt_stn.priors.mean_θ, opt_stn.priors.std_θ;
                   method=opt_stn.method, optimize_mode=opt_stn.optimize_mode,
                   objective=opt_stn.objective, regularization=opt_stn.regularization,
                   loop_num=f_stn.loop_num, peak_shift_iter=opt_stn.maxiter,
                   mod_peak_num=f_stn.mod_peak_num, peak_mod_mean=f_stn.peak_mod_mean,
                   peak_mod_std=f_stn.peak_mod_std, peak_mod_iter=f_stn.peak_mod_iter,
                   analytic_peak_mod=f_stn.analytic_peak_mod,
                   λ=opt_stn.λ, verbose=opt_stn.verbose, tol=opt_stn.tol)
end

function get_top_nodes(nodes::AbstractVector{Node}, k::Int)
    residuals = [norm(nodes[i].residual) for i in eachindex(nodes)]
    if k > length(nodes)
        return @view nodes[sortperm(residuals)]
    end
    return @view nodes[sortperm(residuals)[1:k]]
end

function get_max_id(LT::AbstractTree)
    ids = [LT.nodes[i].id for i in eachindex(LT.nodes)]
    maximum(ids)
end