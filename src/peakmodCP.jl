struct PeakModCP{T, V<:AbstractVector{T}, C, L, CL, P, K, M, N} <: AbstractPhase
    cl::C # crystal object
    origin_cl::CL # save for later comparison
    peaks::V # Vector of peak object
    peak_int::L

    id::Int64 # Just index
    name::String # For printing

    act::K # Activation
    σ::M # Width of peaks
    profile::P # Peak profile
               # Do import from PhaseMapping if
               # other peak profiles are needed
    norm_constant::N
end

Base.Bool(CP::PeakModCP) = true
Base.Bool(CPs::AbstractVector{<:PeakModCP}) = true
get_param_nums(CP::PeakModCP) = length(CP.peak_int)
# get_param_nums(CPs::AbstractVector{<:PeakModCP}) = sum(get_param_nums.(CPs))

function get_intensity(p::Peak)
    return p.I
end

function PeakModCP(_stn::String, wid_init::Real=.1, allow_peak_mod::Int64=10,
                 profile::PeakProfile=PseudoVoigt(0.5))
    f = split(_stn, '\n')
    lattice_info = split(f[1], ',')
    id = parse(Int64, lattice_info[1])
    crystal = get_crystal(cast(lattice_info[4:end], Float64))
    peaks, norm_constant = get_peaks(f[2:end])
    sort!(peaks, rev=true)
    allowed_peak = min(allow_peak_mod, length(peaks))
    peak_int = zeros(Float64, allowed_peak)
    @. peak_int = get_intensity(peaks)[1:allowed_peak]
    name = String(lattice_info[2])
    act = 1.0

    PeakModCP(crystal, crystal, peaks, peak_int,
               id, name, act, wid_init, profile, norm_constant)
end

function get_free_params(CP::PeakModCP)
    CP.peak_int
end

get_param_nums(CP::PeakModCP) = length(CP.peak_int)

function PeakModCP(CP::CrystalPhase, allow_peak_mod::Int64=10)
    allowed_peak = min(allow_peak_mod, length(CP.peaks))
    peak_int = zeros(Float64, allowed_peak)
    @. peak_int = get_intensity(CP.peaks[1:allowed_peak])
    PeakModCP(CP.cl, CP.origin_cl, CP.peaks, peak_int, CP.id, CP.name, CP.act, CP.σ, CP.profile, CP.norm_constant)
end

function PeakModCP(CP::CrystalPhase)
    peak_int = zeros(Float64, length(CP.peaks))
    @. peak_int = get_intensity(CP.peaks)
    PeakModCP(CP.cl, CP.origin_cl, CP.peaks, peak_int, CP.id, CP.name, CP.act, CP.σ, CP.profile, CP.norm_constant)
end

function CrystalPhase(CP::PeakModCP)
    for i in 1:length(CP.peak_int)
        CP.peaks[i] = Peak(CP.peaks[i].h, CP.peaks[i].k, CP.peaks[i].l,
                        CP.peaks[i].q, CP.peak_int[i])
    end
    CrystalPhase(CP.cl, CP.origin_cl, CP.peaks, CP.id, CP.name, CP.act, CP.σ, CP.profile, CP.norm_constant)
end

function PeakModCP(CP::PeakModCP, θ::AbstractVector) 
    CP = PeakModCP(CP.cl, CP.origin_cl, CP.peaks, θ, 
                CP.id, CP.name,
                CP.act, CP.σ, CP.profile, CP.norm_constant)
    return CP
end

function reconstruct_CPs!(θ::AbstractVector, CPs::AbstractVector{<:PeakModCP})
    start = 1
    new_CPs = Vector{PeakModCP}(undef, length(CPs))
    for i in eachindex(CPs)
		new_CPs[i] = PeakModCP(CPs[i], θ[start:start + get_param_nums(CPs[i])-1])
		start += get_param_nums(CPs[i])
	end
    θ[start:end], new_CPs
end

function (CP::PeakModCP)(x::Real) 
    y = zero(x)
    n = length(CP.peak_int)
    @simd for i in eachindex(CP.peaks)
        q = (CP.cl)(CP.peaks[i]) * 10 # account for unit difference
        if i <= n
            y += CP.act * CP.peak_int[i] * CP.profile((x-q)/CP.σ) # Main bottle neck
        else
            y += CP.act * CP.peaks[i].I * CP.profile((x-q)/CP.σ)
        end
    end
    y
end

function evaluate(CP::PeakModCP, θ::AbstractVector, x::AbstractVector) 
    PeakModCP(CP, θ).(x)
end

function evaluate_residual!(CP::PeakModCP, x::AbstractVector, r::AbstractVector) 
    n = length(CP.peak_int)
    for i in eachindex(CP.peak_int)
        q = (CP.cl)(CP.peaks[i]) * 10 # account for unit difference
        @. r -= CP.act * CP.peak_int[i] * CP.profile((x-q)/CP.σ) # Main bottle neck
    end
    for i in n+1:length(CP.peaks)
        q = (CP.cl)(CP.peaks[i]) * 10
        @. r -= CP.act * CP.peaks[i].I * CP.profile((x-q)/CP.σ)
    end
    r
end

function evaluate!(y::AbstractVector, CP::PeakModCP, θ::AbstractVector,  x::AbstractVector) 
    PeakModCP(CP, θ)(x, y)
end

function evaluate!(y::AbstractVector, CP::PeakModCP, x::AbstractVector) 
    n = length(CP.peak_int)
    for i in eachindex(CP.peaks)
        q = (CP.cl)(CP.peaks[i]) * 10 # account for unit difference
        if i <= n
            @. y += CP.act * CP.peak_int[i] * CP.profile((x-q)/CP.σ) # Main bottle neck
        else
            @. y += CP.act * CP.peaks[i].I * CP.profile((x-q)/CP.σ)
        end
    end
    y
end
function (CP::PeakModCP)(x::AbstractVector, y::AbstractVector)
    evaluate!(y, CP, x)
end