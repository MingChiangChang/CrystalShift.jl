const Background = Union{BackgroundModel, Nothing}
const OptionalPhases = Union{AbstractVector{<:AbstractPhase}, Nothing}

struct PhaseModel{A<:OptionalPhases, W<:OptionalPhases, B<:Background} 
    CPs::A
    wildcard::W
    background::B
end

"""
Update Phasemodel object with new parameters
Update the Phases and background seperately
The parameters for background always comes last
"""

function PhaseModel()
    PhaseModel(CrystalPhase[], nothing, nothing)
end

function PhaseModel(CPs::AbstractVector{<:AbstractPhase})
    PhaseModel(CPs, nothing, nothing)
end

function PhaseModel(W::AbstractVector{<:Wildcard})
    PhaseModel(nothing, W, nothing)
end

function PhaseModel(W::Wildcard)
    PhaseModel(nothing, [W], nothing)
end

function PhaseModel(BG::Background)
    PhaseModel(nothing, nothing, BG)
end

# function PeakModCP(pm::PhaseModel)
#     PeakModCP.(pm.CPs)
# end

function CrystalPhase(pm::PhaseModel)
    pm.CPs
end

PhaseModel(CP::AbstractPhase) = PhaseModel([CP])
PhaseModel(CP::AbstractPhase, BG::Background) = PhaseModel([CP], nothing, BG)
PhaseModel(CP::AbstractPhase, W::Wildcard, BG::Background) = PhaseModel([CP], [W], BG)
PhaseModel(W::Wildcard, BG::Background) = PhaseModel(nothing, [W], BG)

function PhaseModel(PM::PhaseModel, θ::AbstractVector)
    # bg_param_num = get_param_nums(PM.background)
    # θ_cp = θ[1:end - bg_param_num ]
    # θ_bg = θ[end - bg_param_num + 1 : end]
    θ, CPs = reconstruct_CPs!(θ, PM.CPs)
    θ, W = reconstruct_WC!(θ, PM.wildcard)
    θ, BG  = reconstruct_BG!(θ, PM.background)
    PhaseModel(CPs, W, BG)
end

Base.size(PM::PhaseModel) = size(PM.CPs)
Base.size(PM::PhaseModel, dim::Integer) = size(PM.CPs, dim)
Base.length(PM::PhaseModel) = length(PM.CPs)
Base.iterate(PM::PhaseModel) = iterate(PM.CPs)
Base.iterate(PM::PhaseModel, state) = iterate(PM.CPs, state)
Base.:(==)(PM1::PhaseModel, PM2::PhaseModel) = [p.id for p in PM1.CPs] == [p.id for p in PM2.CPs]

function reconstruct_CPs!(θ::AbstractVector, CPs::AbstractVector{<:AbstractPhase})
    start = 1
    new_CPs = Vector{AbstractPhase}(undef, length(CPs))
    for i in eachindex(CPs)
		new_CPs[i] = CrystalPhase(CPs[i], θ[start:start + get_param_nums(CPs[i])-1])
		start += get_param_nums(CPs[i])
	end
    θ[start:end], new_CPs
end

reconstruct_CPs!(θ::AbstractVector, CPs::Nothing) = θ, nothing
reconstruct_WC!(θ::AbstractVector, W::Nothing) = θ, nothing

function reconstruct_WC!(θ::AbstractVector, W::AbstractVector{<:Wildcard})
    start = 1
    new_W = Vector{Wildcard}(undef, length(W))
    for i in eachindex(W)
		new_W[i] = Wildcard(W[i], θ[start:start + get_param_nums(W[i])-1])
		start += get_param_nums(W[i])
	end
    θ[start:end], new_W
end

function reconstruct!(pm::PhaseModel, θ::AbstractVector)
    θ, CPs = reconstruct_CPs!(θ, pm.CPs)
    θ, W = reconstruct_WC!(θ, pm.wildcard)
    θ, background = reconstruct_BG!(θ, pm.background)
    isempty(θ) || error("θ should be empty after reconstructing the PhaseModel object.")
    return PhaseModel(CPs, W, background)
end

Nothing(n::Nothing, a::Any) = nothing
get_param_nums(B::Nothing) = 0
get_free_params(B::Nothing) = Float64[]

evaluate!(y::AbstractVector, B::Nothing, θ::AbstractVector, x::AbstractVector) = y
evaluate!(y::AbstractVector, B::Nothing, x::AbstractVector) = y
evaluate_residual!(B::Nothing, θ::AbstractVector, x::AbstractVector, r::AbstractVector) = r
evaluate_residual!(B::Nothing, x::AbstractVector, r::AbstractVector) = r

get_param_nums(PM::PhaseModel) = get_param_nums(PM.CPs) + get_param_nums(PM.wildcard) + get_param_nums(PM.background)
get_free_params(PM::PhaseModel) = vcat(get_free_params(PM.CPs), get_free_params(PM.wildcard), get_free_params(PM.background))

get_phase_ids(PM::PhaseModel) = get_phase_ids(PM.CPs)
get_phase_ids(CPs::AbstractVector) = [p.id for p in CPs]

# TODO: Combine these function with the ones in CrystalPhase
function (PM::PhaseModel)(x::AbstractVector, y::AbstractVector)
    evaluate!(y, PM, x)
end

function (PM::PhaseModel)(x::AbstractVector)
    y = zero(x)
    PM(x, y)
end

function evaluate(pm::PhaseModel, θ::AbstractVector, x::AbstractVector)
    y = zeros(Real, size(x))
    evaluate!(y, pm, θ, x)
end

function evaluate!(y::AbstractVector, PM::PhaseModel, x::AbstractVector)
    evaluate!(y, PM.CPs, x)
    evaluate!(y, PM.wildcard, x)
    evaluate!(y, PM.background, x)
    y
end

function evaluate!(y::AbstractVector, PM::PhaseModel, θ::AbstractVector,
                   x::AbstractVector)
    PhaseModel(PM, θ)(x, y)
end

function evaluate_residual!(PM::PhaseModel, θ::AbstractVector,
                            x::AbstractVector, r::AbstractVector)
    evaluate_residual!(PhaseModel(PM, θ), x, r)
end

function evaluate_residual!(PM::PhaseModel, x::AbstractVector, r::AbstractVector)
    evaluate_residual!(PM.CPs, x, r)
    evaluate_residual!(PM.wildcard, x, r)
    evaluate_residual!(PM.background, x, r)
end

function get_PeakModCP(PM::PhaseModel, x::AbstractVector, allowed_num::Int64)
    IMs = PeakModCP.(PM.CPs, (x, ), (allowed_num, ))
    evaluate!(IMs[1].const_basis, PM.background, x)
    IMs
end