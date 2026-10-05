const ScalarOrVecInt = Union{Integer, AbstractVector{<:Integer}}

"""
    FullOptimizeSettings(; loop_num=8, mod_peak_num=32, peak_mod_mean=[1.],
                           peak_mod_std=[.5], peak_mod_iter=32, analytic_peak_mod=true)

Peak-height refinement options for tree search. Passing one as `full_opt_stn` to
`TreeSearchSettings` makes the search optimize each node with `full_optimize!`
(which also refits peak heights) instead of `optimize!`. The fields are the
`full_optimize!` keywords of the same names; the remaining options (priors,
method, objective, `maxiter` as `peak_shift_iter`, ...) come from `opt_stn`.
"""
struct FullOptimizeSettings
    loop_num::Int
    mod_peak_num::Int
    peak_mod_mean::Vector{Float64}
    peak_mod_std::Vector{Float64}
    peak_mod_iter::Int
    analytic_peak_mod::Bool
end

function FullOptimizeSettings(; loop_num::Int = 8, mod_peak_num::Int = 32,
                              peak_mod_mean::AbstractVector = [1.], peak_mod_std::AbstractVector = [.5],
                              peak_mod_iter::Int = 32, analytic_peak_mod::Bool = true)
    FullOptimizeSettings(loop_num, mod_peak_num, peak_mod_mean, peak_mod_std, peak_mod_iter, analytic_peak_mod)
end

struct TreeSearchSettings{V} <: AbstractTreeSearchSettings
    depth::Integer
    k::ScalarOrVecInt
    amorphous::Bool # Amorphous
    background::Bool
    background_length::Real
    default_phase::Union{Nothing, CrystalPhase}
    opt_stn::OptimizationSettings{V}
    full_opt_stn::Union{Nothing, FullOptimizeSettings} # nothing: plain optimize!

    function TreeSearchSettings(depth::Integer, k::ScalarOrVecInt, amorphous::Bool, background::Bool,
                                background_length::Real, default_phase::Union{Nothing, CrystalPhase},
                                opt_stn::OptimizationSettings{V},
                                full_opt_stn::Union{Nothing, FullOptimizeSettings} = nothing) where V
        if !isnothing(full_opt_stn) && opt_stn.optimize_mode isa _WithUncer
            error("full_opt_stn does not support the WithUncer optimize_mode")
        end
        new{V}(depth, k, amorphous, background, background_length, default_phase, opt_stn, full_opt_stn)
    end
end

function TreeSearchSettings{Float64}(; full_opt_stn::Union{Nothing, FullOptimizeSettings} = nothing)
    opt_stn = OptimizationSettings{Float64}()
    TreeSearchSettings(2, 3, false, false, 5., nothing, opt_stn, full_opt_stn)
end

function TreeSearchSettings{T}(depth::Integer,
     k::ScalarOrVecInt,
     amorphous::Bool,
     background::Bool,
     background_length::Real,
     opt_stn::OptimizationSettings{T};
     full_opt_stn::Union{Nothing, FullOptimizeSettings} = nothing) where T
    TreeSearchSettings(depth, k, amorphous, background, background_length, nothing, opt_stn, full_opt_stn)
end

function TreeSearchSettings{T}(depth::Integer,
     k::ScalarOrVecInt,
     amorphous::Bool,
     background::Bool,
     background_length::Real,
     default_phase::Union{Nothing, CrystalPhase},
     opt_stn::OptimizationSettings{T};
     full_opt_stn::Union{Nothing, FullOptimizeSettings} = nothing) where T
    TreeSearchSettings(depth, k, amorphous, background, background_length, default_phase, opt_stn, full_opt_stn)
end

struct MPTreeSearchSettings{V} <: AbstractTreeSearchSettings
    depth::Integer
    k::ScalarOrVecInt
    mp_top_k::ScalarOrVecInt
    amorphous::Bool # Amorphous
    background::Bool
    background_length::Real
    default_phase::Union{Nothing, CrystalPhase}
    opt_stn::OptimizationSettings{V}
end

function MPTreeSearchSettings{Float64}()
    opt_stn = OptimizationSettings{Float64}()
    MPTreeSearchSettings(2, 3, 2, false, false, 5., nothing, opt_stn)
end

function MPTreeSearchSettings{T}(
    depth::Integer,
    k::ScalarOrVecInt,
    mp_top_k::ScalarOrVecInt,
    amorphous::Bool,
    background::Bool,
    background_length::Real,
    opt_stn::OptimizationSettings{T}) where T
    MPTreeSearchSettings(depth, k, mp_top_k, amorphous, background, background_length, nothing, opt_stn)
end