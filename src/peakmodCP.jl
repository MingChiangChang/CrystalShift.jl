struct PeakModCP{T, V<:AbstractMatrix{T}, L<:AbstractVector{T}, K} <: AbstractPhase
    basis::V
    const_basis::L
    peak_int::K
end

# test
function PeakModCP(CP::CrystalPhase, x::AbstractVector, allowed_num::Int64)
    peak_num = min(allowed_num, size(CP.peaks, 1))
    basis = zeros(Float64, size(x, 1), peak_num)
    evaluate!(basis, CP, CP.peaks[1:peak_num], x)
    const_basis = zeros(Float64, size(x, 1))
    for i in peak_num+1:size(CP.peaks, 1)
        evaluate!(const_basis, CP, CP.peaks[i], x)
    end
    peak_int =  ones(Float64, peak_num)
    return PeakModCP(basis, const_basis, peak_int)
end

get_param_nums(IM::PeakModCP) = length(IM.peak_int)
get_free_params(IM::PeakModCP) = IM.peak_int

function change_peak_int!(CP::CrystalPhase, IM::PeakModCP)
    change_peak_int!(CP, IM.peak_int)
end

function change_peak_int!(CP::CrystalPhase, peak_int::AbstractVector)
    for i in 1:length(peak_int)
        CP.peaks[i] = Peak(CP.peaks[i].h, CP.peaks[i].k, CP.peaks[i].l,
                        CP.peaks[i].q, CP.peaks[i].I * peak_int[i])
    end
end

function reconstruct_CPs!(θ::AbstractVector, CPs::AbstractVector{<:PeakModCP})
    start = 1
    new_CPs = Vector{PeakModCP}(undef, length(CPs))
    @simd for i in eachindex(CPs)
		new_CPs[i] = PeakModCP(CPs[i], θ[start:start + get_param_nums(CPs[i])-1])
		start += get_param_nums(CPs[i])
	end
    θ[start:end], new_CPs
end


function evaluate!(y::AbstractVector, IM::PeakModCP, x::AbstractVector) 
    y .+= IM.basis * IM.peak_int
    @. y += IM.const_basis
    y
end

function evaluate!(y::AbstractVector, IM::PeakModCP, θ::AbstractVector,
                  x::AbstractVector)
    @. IM.peak_int = θ
    y .+= IM.basis * IM.peak_int
    @. y += IM.const_basis
    y
end

function evaluate_residual!(IM::PeakModCP, x::AbstractVector, r::AbstractVector)
    r .-= IM.basis * IM.peak_int
    @. r -= IM.const_basis
    r
end

function evaluate_residual!(IM::PeakModCP, θ::AbstractVector,
                            x::AbstractVector, r::AbstractVector)
    r .-= IM.basis[:,i] *θ
    @. r -= IM.const_basis
    r
end

function PeakModCP(IM::PeakModCP, θ::AbstractVector)
    PeakModCP(IM.basis, IM.const_basis, θ)
end