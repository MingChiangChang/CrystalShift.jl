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

# CrystalPhase(CP, θ) shares the peaks vector with CP, so give a phase its own copy
# before modifying peak intensities in place with change_peak_int!
copy_peaks(CP::CrystalPhase) = with_peaks(CP, CP.peaks)

# CP (lattice, activation, width, profile) with its own copy of `peaks`
function with_peaks(CP::CrystalPhase, peaks::AbstractVector)
    CrystalPhase(CP.cl, CP.origin_cl, copy(peaks), CP.param_num, CP.id, CP.name,
                 CP.act, CP.σ, CP.profile, CP.norm_constant)
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
    r .-= IM.basis * θ
    @. r -= IM.const_basis
    r
end

function PeakModCP(IM::PeakModCP, θ::AbstractVector)
    PeakModCP(IM.basis, IM.const_basis, θ)
end

"""
    LinearPeakMod(IMs, y, mean_θ, std_θ)

Peak-height subproblem of `full_optimize!` in a form that is solved without
automatic differentiation. The model y ≈ B*w + c is linear in the height factors
w = exp(u), so the least-squares term only depends on the sufficient statistics
BᵀB, Bᵀ(y - c) and ‖y - c‖², which are precomputed here once. Each iteration of
`optimize!(::LinearPeakMod, ...)` then only works with n×n quantities
(n = total number of free peak heights) instead of length(x)-sized arrays.
B is the column-concatenation of `IM.basis` and c the sum of `IM.const_basis`
over all `IMs`. Priors are the same as in the `PeakModCP` route (`extend_priors`).
"""
struct LinearPeakMod{T<:Real, P<:PeakModCP}
    IMs::Vector{P}
    BtB::Matrix{T}       # BᵀB
    Btr::Vector{T}       # Bᵀ(y - c)
    rtr::T               # ‖y - c‖²
    mean_log_θ::Vector{T}
    std_θ::Vector{T}
end

function LinearPeakMod(IMs::AbstractVector{<:PeakModCP}, y::AbstractVector,
                       mean_θ::AbstractVector, std_θ::AbstractVector)
    B = reduce(hcat, [IM.basis for IM in IMs])
    r = y .- sum(IM.const_basis for IM in IMs)
    full_mean_θ, full_std_θ = extend_priors(mean_θ, std_θ, IMs)
    LinearPeakMod(collect(IMs), B'B, B'r, dot(r, r), log.(full_mean_θ), full_std_θ)
end

get_param_nums(P::LinearPeakMod) = length(P.Btr)

# Split the stacked height factors back into one PeakModCP per phase
function reconstruct_IMs(P::LinearPeakMod, w::AbstractVector)
    start = 1
    map(P.IMs) do IM
        n = get_param_nums(IM)
        new_IM = PeakModCP(IM, w[start:start+n-1])
        start += n
        new_IM
    end
end