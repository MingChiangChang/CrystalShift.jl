struct int_mod{T, V<:AbstractMatrix{T}, L<:AbstractVector{T}, K} <: AbstractPhase
    basis::V
    const_basis::L
    peak_int::K
end

function int_mod(CP::CrystalPhase, x::AbstractVector, allowed_num::Int64)
    peak_num = min(allowed_num, size(CP.peaks, 1))
    basis = zeros(Float64, size(x, 1), peak_num)
    evaluate!(basis, CP, CP.peaks[1:peak_num], x)
    const_basis = zeros(Float64, size(x, 1))
    for i in peak_num+1:size(CP.peaks, 1)
        evaluate!(const_basis, CP, CP.peaks[i], x)
    end
    peak_int =  ones(Float64, peak_num)
    return int_mod(basis, const_basis, peak_int)
end

get_param_nums(IM::int_mod) = length(IM.peak_int)
get_free_params(IM::int_mod) = IM.peak_int

function change_peak_int!(CP::CrystalPhase, IM::int_mod)
    change_peak_int!(CP, IM.peak_int)
end

function change_peak_int!(CP::CrystalPhase, peak_int::AbstractVector)
    for i in 1:length(peak_int)
        CP.peaks[i] = Peak(CP.peaks[i].h, CP.peaks[i].k, CP.peaks[i].l,
                        CP.peaks[i].q, CP.peaks[i].I * peak_int[i])
    end
end

function reconstruct_CPs!(θ::AbstractVector, CPs::AbstractVector{<:int_mod})
    start = 1
    new_CPs = Vector{int_mod}(undef, length(CPs))
    @simd for i in eachindex(CPs)
		new_CPs[i] = int_mod(CPs[i], θ[start:start + get_param_nums(CPs[i])-1])
		start += get_param_nums(CPs[i])
	end
    θ[start:end], new_CPs
end


function evaluate!(y::AbstractVector, IM::int_mod, x::AbstractVector) 
    @simd for i in eachindex(IM.peak_int)
        y .+= IM.basis[:,i] .* IM.peak_int[i]
    end
    y += IM.const_basis
    y
end

function evaluate!(y::AbstractVector, IM::int_mod, θ::AbstractVector,
                  x::AbstractVector) 
    @. IM.peak_int = θ
    @simd for i in eachindex(IM.peak_int)
        @. y += IM.basis[:,i] * IM.peak_int[i]
    end
    @. y += IM.const_basis
    y
end

function evaluate_residual!(IM::int_mod, x::AbstractVector, r::AbstractVector)
    @simd for i in eachindex(IM.peak_int)
        @. r -= IM.basis[:,i] * IM.peak_int[i]
    end
    @. r -= IM.const_basis
    r
end

function int_mod(IM::int_mod, θ::AbstractVector)
    int_mod(IM.basis, IM.const_basis, θ)
end