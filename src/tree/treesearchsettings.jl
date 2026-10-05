const ScalarOrVecInt = Union{Integer, AbstractVector{<:Integer}}

struct TreeSearchSettings{V} <: AbstractTreeSearchSettings
    depth::Integer
    k::ScalarOrVecInt
    amorphous::Bool # Amorphous
    background::Bool
    background_length::Real
    default_phase::Union{Nothing, CrystalPhase}
    opt_stn::OptimizationSettings{V}
end

function TreeSearchSettings{Float64}()
    opt_stn = OptimizationSettings{Float64}()
    TreeSearchSettings(2, 3, false, false, 5., nothing, opt_stn)
end

function TreeSearchSettings{T}(depth::Integer,
     k::ScalarOrVecInt,
     amorphous::Bool,
     background::Bool,
     background_length::Real,
     opt_stn::OptimizationSettings{T}) where T
    TreeSearchSettings(depth, k, amorphous, background, background_length, nothing, opt_stn)
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