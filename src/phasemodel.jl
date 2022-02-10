const Background = Union{BackgroundModel, Nothing}

struct PhaseModel{A<:AbstractVector{CrystalPhase}, B<:Background} 
    CPs::A
    background::B
end

"""
Update Phasemodel object with new parameters
Update the Phases and background seperately
The parameters for background always comes last
"""
function PhaseModel(CPs::AbstractVector{<:CrystalPhase})
    PhaseModel(CPs, nothing)
end

PhaseModel(CP::CrystalPhase) = PhaseModel([CP])
PhaseModel(CP::CrystalPhase, BG::Background) = PhaseModel([CP], BG)

function PhaseModel(PM::PhaseModel, θ::AbstractVector)
    # bg_param_num = get_param_nums(PM.background)
    # θ_cp = θ[1:end - bg_param_num ]
    # θ_bg = θ[end - bg_param_num + 1 : end]
    θ, CPs = reconstruct_CPs!(θ, PM.CPs)
    θ, BG  = reconstruct_BG!(θ, PM.background)
    PhaseModel(CPs, BG)
end

function reconstruct!(pm::PhaseModel, θ::AbstractVector)
    θ, CPs = reconstruct_CPs!(θ, pm.CPs)
    θ, background = reconstruct_BG!(θ, pm.background)
    isempty(θ) || error("θ should be empty after reconstructing the PhaseModel object.")
    return PhaseModel(CPs, background)
end

Nothing(n::Nothing, a::Any) = nothing
get_param_nums(B::Nothing) = 0
get_free_params(B::Nothing) = Float64[]
evaluate!(y::AbstractVector, BG::Nothing, x::AbstractVector) = y
evaluate_residual!(y::AbstractVector, BG::Nothing, x::AbstractVector) = y

get_param_nums(PM::PhaseModel) = get_param_nums(PM.CPs) + get_param_nums(PM.background)
get_free_params(PM::PhaseModel) = vcat(get_free_params(PM.CPs), get_free_params(PM.background))

# TODO: Combine these function with the ones in CrystalPhase
function (PM::PhaseModel)(x::AbstractVector, y::AbstractVector)
    evaluate!(y, PM, x)
end

function evaluate!(y::AbstractVector, PM::PhaseModel, x::AbstractVector)
    evaluate!(y, PM.CPs, x)
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
    r
end

function evaluate_residual!(PM::PhaseModel, x::AbstractVector, r::AbstractVector)
    evaluate_residual!(PM.CPs, x, r)
    evaluate_residual!(PM.background, x, r)
    r
end

function evaluate_residual!(CP::CrystalPhase, θ::AbstractVector,
              x::AbstractVector, r::AbstractVector)
    evaluate_residual!(CrystalPhase(CP, θ), x, r)
    r
end

function evaluate_residual!(CPs::AbstractVector{<:CrystalPhase},
              x::AbstractVector, r::AbstractVector)
    @simd for i in eachindex(CPs)
        evaluate_residual!(CPs[i], x, r)
    end
    r
end

function evaluate_residual!(CP::CrystalPhase, x::AbstractVector, r::AbstractVector)
    @simd for i in eachindex(CP.peaks)
        q = (CP.cl)(CP.peaks[i]) * 10 # account for unit difference
        @. r -= CP.act * CP.peaks[i].I * CP.profile((x-q)/CP.σ) # Main bottle neck
    end
    r
end

function evaluate_residual!(CPs::AbstractVector{<:CrystalPhase},
             θ::AbstractVector, x::AbstractVector, r::AbstractVector)
    s = 1
    for i in eachindex(CPs)
        num_of_param = get_param_nums(CPs[i])
        θ_temp = @view θ[s : s+num_of_param-1]
        evaluate_residual!(CPs[i], θ_temp, x, r)
        s += num_of_param
    end
    r
end

# Without preallocation, useful at times....
function (CP::CrystalPhase)(x::Real)
    y = zero(x)
    @simd for i in eachindex(CP.peaks)
        q = (CP.cl)(CP.peaks[i]) * 10 # account for unit difference
        y += CP.act * CP.peaks[i].I * CP.profile((x-q)/CP.σ) # Main bottle neck
    end
    y
end

function (CPs::AbstractVector{<:CrystalPhase})(x::Real)
    y = zero(x)
    @simd for i in eachindex(CPs) #
        y += CPs[i](x)
    end
    y
end


function evaluate(CP::CrystalPhase, θ::AbstractVector,
                      x::AbstractVector)
    CrystalPhase(CP, θ).(x)
end

function evaluate(CPs::AbstractVector{<:CrystalPhase},
                  θ::AbstractVector, x::AbstractVector)
    y = zeros(size(x))
    s = 1
    for i in eachindex(CPs)
        num_of_param = get_param_nums(CPs[i])
        θ_temp = @view θ[s:s+num_of_param-1]
        y += evaluate(CPs[i], θ_temp, x)
        s += num_of_param
    end
    y
end