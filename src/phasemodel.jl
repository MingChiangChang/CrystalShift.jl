struct Phasemodel
    CPs::AbstractVector{<:CrystalPhase}
    background::BackgroundModel
end

function Phasemodel(PM::Phasemodel, θ::AbstractVector)
    
end

get_param_nums(PM::Phasemodel) = get_param_nums(PM.CPs) + get_param_nums(PM.background)
get_free_params(PM::Phasemodel) = vcat(get_free_params(PM.CPs) + get_free_params(PM.background))

function (PM::Phasemodel)(x::AbstractVector, y::AbstractVector)
    evaluate!(y, PM, x)
end

function evaluate!(y::AbstractVector, PM::Phasemodel, x::AbstractVector)
    evaluate!(y, PM.CPs, x)
    evaluate!(y, PM.background, x)
    y
end

function evaluate!(y::AbstractVector, PM::Phasemodel, θ::AbstractVector,
                   x::AbstractVector)
    Phasemodel(PM, θ)(x, y)
end

function evaluate_residual!(CP::CrystalPhase, θ::AbstractVector,
              x::AbstractVector, r::AbstractVector)
    evaluate_residual!(CrystalPhase(CP, θ), x, r)
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