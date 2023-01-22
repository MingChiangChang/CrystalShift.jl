abstract type AbstractPhase end

struct CrystalPhase{T, V<:AbstractVector{T}, C, CL, P, K, M, N} <: AbstractPhase
    cl::C # crystal object
    origin_cl::CL # save for comparison
    peaks::V # Vector of peak object

    id::Int64 # For indexing phases
    name::String # Phase name

    act::K # Activation parameter
    σ::M # Peak width parameter
    profile::P # Peak profile
    norm_constant::N # Normalization constant. Stored when phase is created.
                     # Use to calculate phase fraction
end

function Base.show(io::IO, CP::CrystalPhase)
    println("Phase name: $(CP.name), ID: $(CP.id)")
    println("Optimization parameters:")
    println("Activation: $(CP.act), Peak width: $(CP.σ)")
    println("Normalization: $(CP.norm_constant)")
    println("Lattice information:")
    println(CP.cl)
    return true
end

Base.Bool(CP::AbstractPhase) = true
Base.Bool(CP::CrystalPhase) = true
Base.Bool(CPs::AbstractVector{<:CrystalPhase}) = true
get_param_nums(CP::CrystalPhase) = CP.cl.free_param + 2 + get_param_nums(CP.profile)
get_param_nums(APs::AbstractVector{<:AbstractPhase}) = sum(get_param_nums.(APs))
# get_param_nums(CPs::AbstractVector{<:CrystalPhase}) = sum(get_param_nums.(CPs))

## Constrcutors
function CrystalPhase(io::IOStream,
                      width_init::Real=.1,
                      profile::PeakProfile=FixedPseudoVoigt(0.5))
    if Sys.iswindows()
        s = split(read(io, String), "#\r\n")
    else
        s = split(read(io, String), "#\n")
    end

    CrystalPhase.(String.(s[1:end-1]), (width_init), (profile,))
end

function CrystalPhase(_stn::String, wid_init::Real=.1,
                      profile::PeakProfile=PseudoVoigt(0.5))
    f = split(_stn, '\n')
    lattice_info = split(f[1], ',')
    id = parse(Int64, lattice_info[1])
    crystal = get_crystal(cast(lattice_info[4:end], Float64))
    peaks, norm_constant = get_peaks(f[2:end])
    sort!(peaks, rev=true)
    name = String(lattice_info[2])
    act = 1.0

    CrystalPhase(crystal, crystal, peaks, id, name, act, wid_init, profile, norm_constant)
end

function CrystalPhase(CP::CrystalPhase, wid_init::Real=.1, profile::PeakProfile=FixedPseudoVoigt(0.5))
    cl_params = get_free_lattice_params(CP.origin_cl)
    profile_param_num = get_param_nums(CP.profile)
    cl = get_intrinsic_crystal_type(typeof(CP.origin_cl))
    profile_type = get_intrinsic_profile_type(typeof(CP.profile))
    t = eltype(cl_params)
    c = CrystalPhase(cl{t}(cl_params...), CP.origin_cl, CP.peaks, CP.id, CP.name,
                0.1, wid_init, profile, CP.norm_constant)
    return c
end

function CrystalPhase(CP::CrystalPhase, θ::AbstractVector)
    fp = CP.cl.free_param
    profile_param_num = get_param_nums(CP.profile)
    cl = get_intrinsic_crystal_type(typeof(CP.cl))
    profile_type = get_intrinsic_profile_type(typeof(CP.profile))
    t = eltype(θ)
    if profile_param_num > 0
        c = CrystalPhase(cl{t}(θ[1:fp]...), CP.origin_cl, CP.peaks, CP.id, CP.name,
        # θ[fp+1], θ[fp+2], CP.profile, CP.norm_constant)
                    θ[fp+1], θ[fp+2], profile_type{t}(θ[fp+3:fp+2+profile_param_num]...),
                    CP.norm_constant)
    else
        c = CrystalPhase(cl{t}(θ[1:fp]...), CP.origin_cl, CP.peaks, CP.id, CP.name,
        θ[fp+1], θ[fp+2], CP.profile, CP.norm_constant)
    end
    return c
end

function CrystalPhase(path::AbstractString,
                     q_range::Tuple, id::Integer,
                     wid_init::Real=.1,
                     profile::PeakProfile=FixedPseudoVoigt(0.5))
    cif = CifParser(path)
    lattice = CIFFile(path).SGLattice()

    cif_dict = cif.as_dict()
    key = _get_key(cif_dict)
    info_dict = cif_dict[key]

    # Getting names
    phase_name = _get_phase_name(info_dict)
    crystal_sys = _get_crystal_system(info_dict)
    name = phase_name * "_" * crystal_sys
    lattice_params = _get_lattice_parameters(info_dict)
    crystal = get_crystal(lattice_params)

    # Getting peaks
    peaks, norm_constant = _get_peaks(lattice, q_range )
    sort!(peaks, rev=true)

    CrystalPhase(crystal, crystal, peaks, id, name, 1.0,
                 wid_init, profile, norm_constant)
end

# Helper function
function _get_peaks(lattice, q_range)
    xtal = Xtal("test", lattice)
    xrd = PowderDiffraction(xtal).data

    peaks = Vector{Peak}()
    for peak in keys(xrd)
        q = xrd[peak]["qpos"] * 10
        I = xrd[peak]["r"]
        if q_range[1] < q < q_range[2] && I>0.0001
           push!(peaks, Peak(peak[1], peak[2], peak[3], q, I))
        end
    end
    normalize_peaks!(peaks)
end

function get_intrinsic_profile_type(profile_type::Type)
    if profile_type <: PseudoVoigt
        return PseudoVoigt
    elseif profile_type <: Lorentz
        return Lorentz
    elseif profile_type <: Gauss
        return Gauss
    end
end

function get_intrinsic_crystal_type(cl::Type)
    if cl <: Cubic
        return Cubic
    elseif cl <: Tetragonal
        return Tetragonal
    elseif cl <: Monoclinic
        return Monoclinic
    elseif cl <: Hexagonal
        return Hexagonal
    elseif cl <: Orthorhombic
        return Orthorhombic
    elseif cl <: Rhombohedral
        return Rhombohedral
    elseif cl <: Triclinic
        return Triclinic
    end
end

function get_moles(CP::CrystalPhase)
    CP.act/CP.norm_constant
end

function get_fraction(CPs::AbstractVector{<:CrystalPhase})
    moles = zeros(size(CPs, 1))
    @. moles = get_moles(CPs)

    moles./sum(moles)
end

function get_free_params(CP::CrystalPhase)
    cl_params = get_free_lattice_params(CP.cl)
    parms = Vector{eltype(cl_params)}(undef, get_param_nums(CP))
    parms[1:CP.cl.free_param] .= get_free_lattice_params(CP.cl)
    parms[CP.cl.free_param+1:CP.cl.free_param+2] .= [CP.act, CP.σ]
    p = get_free_params(CP.profile)
    if !isempty(p)
        parms[CP.cl.free_param+3:end] .= p
    end
    return parms
end

function get_free_params(CPs::AbstractVector{<:AbstractPhase})
    p = Vector{Float64}()
    for cp in CPs
        push!(p, get_free_params(cp)...) # Preallocation?
    end
    p
end

function get_eight_param end
"""Returns a vector of the eight parameters determining a phase structure,
i.e. three sidelengths and three angles for the unit cell, the activation
coefficient, as well as the peak width of the empirical spectrum.

Args:
    - an AbstractPhase object.
    - (optional) a parameter vector θ holding a subset of phase parameters that
        are required to specify a phase, i.e. only one side length for a Cubic.
    - (optional) a real-valued fill angle, fills the unspecified angles with π/2
        by default. For uncertainty computations, we use fill_angle = 0.

Returns:
    - Vector of length eight an real elements, containing the phase parameters.
"""
const FILL_ANGLE = pi / 2

function get_eight_params(CP::AbstractPhase)
    vcat([CP.cl.a, CP.cl.b, CP.cl.c, CP.cl.α, CP.cl.β, CP.cl.γ, CP.act, CP.σ], get_free_params(CP.profile))
end
# ignore fill angle if the method is independent of it:
# get_eight_params(CP::AbstractPhase, θ::AbstractVector, ::Real) = get_eight_params(CP, θ)

function get_eight_params(CP::AbstractPhase, θ::AbstractVector, fill_angle::Real = FILL_ANGLE)
    get_eight_params(CP.cl, θ, fill_angle)
end

function get_eight_params(crystal::Cubic, θ::AbstractVector, fill_angle::Real = FILL_ANGLE)
    [θ[1], θ[1], θ[1], fill_angle, fill_angle, fill_angle, θ[2], θ[3]]
end

function get_eight_params(crystal::Tetragonal, θ::AbstractVector, fill_angle::Real = FILL_ANGLE)
    [θ[1], θ[1], θ[2], fill_angle, fill_angle, fill_angle, θ[3], θ[4]]
end

function get_eight_params(crystal::Orthorhombic, θ::AbstractVector, fill_angle::Real = FILL_ANGLE)
    [θ[1], θ[2], θ[3], fill_angle, fill_angle, fill_angle, θ[4], θ[5]]
end

function get_eight_params(crystal::Rhombohedral, θ::AbstractVector, fill_angle::Real = FILL_ANGLE)
    [θ[1], θ[1], θ[1], θ[2], θ[2], θ[2], θ[3], θ[4]]
end

function get_eight_params(crystal::Hexagonal, θ::AbstractVector, fill_angle::Real = FILL_ANGLE)
    # about third fill angle: pi/2 * 4/3 = 2*pi/3
    [θ[1], θ[1], θ[2], fill_angle, fill_angle, 4/3*fill_angle, θ[3], θ[4]]
end

function get_eight_params(crystal::Monoclinic, θ::AbstractVector, fill_angle::Real = FILL_ANGLE)
    [θ[1], θ[2], θ[3], fill_angle, θ[4], fill_angle, θ[5], θ[6]]
end

get_eight_params(crystal::Triclinic, θ::AbstractVector, fill_angle::Real) = θ

function get_free_lattice_params(CPs::AbstractVector{<:AbstractPhase})
    p = Vector{Float64}()
    for cp in CPs
        push!(p, get_free_lattice_params(cp)...) # Preallocation?
    end
    p
end

function get_eight_params(CP::AbstractVector{<:AbstractPhase}, θ::AbstractVector, fill_angle::Real = FILL_ANGLE)
    start = 1
    params = zeros(Float64, 8*length(CP))
    for i in eachindex(CP)
        free_param_nums = get_param_nums(CP[i])
        params[(i-1)*8+1 : i*8] .= get_eight_params(CP[i], θ[start:start+free_param_nums-1], fill_angle)
        start += free_param_nums
    end
    params
end

get_free_lattice_params(CP::AbstractPhase) = get_free_lattice_params(CP.cl)
get_free_lattice_params(CP::AbstractPhase, θ::AbstractVector) = get_free_lattice_params(CP.cl, θ)
get_free_lattice_params(cl::Cubic, θ::AbstractVector) = [θ[1]]
get_free_lattice_params(cl::Tetragonal, θ::AbstractVector) = [θ[1], θ[3]]
get_free_lattice_params(cl::Hexagonal, θ::AbstractVector) = [θ[1], θ[3]]
get_free_lattice_params(cl::Orthorhombic, θ::AbstractVector) = [θ[1], θ[2], θ[3]]
get_free_lattice_params(cl::Rhombohedral, θ::AbstractVector) = [θ[1], θ[4]]
get_free_lattice_params(cl::Monoclinic, θ::AbstractVector) = [θ[1], θ[2], θ[3], θ[5]]
get_free_lattice_params(cl::Triclinic, θ::AbstractVector) = θ


collect_crystals(CPs::AbstractVector{<:CrystalPhase}) = [CP.cl for CP in CPs]

# Preallocating
# Functor comes in handy but use evaluate! when you can to be clear


# TODO: Vectorize all peak location calculation
function (CP::AbstractPhase)(x::AbstractVector, y::AbstractVector)
    evaluate!(y, CP, x)
end

function (CPs::AbstractVector{<:AbstractPhase})(x::AbstractVector,
                                               y::AbstractVector)
    @simd for i in eachindex(CPs)
        CPs[i](x, y)
    end
    y
end

function evaluate!(y::AbstractVector, CP::AbstractPhase, peak::Peak, x::AbstractVector)
    q = (CP.cl)(peak) * 10 # account for unit difference
    @. y += CP.act * peak.I * CP.profile((x-q)/CP.σ)
    y
end

function evaluate!(y::AbstractMatrix, CP::AbstractPhase, peaks::AbstractVector{<:Peak}, x::AbstractVector)
    peak_locs = (CP.cl).(CP.peaks) .* 10
    @simd for i in eachindex(peaks)
        @. y[:,i] = CP.act * peaks[i].I * CP.profile((x-peak_locs[i])/CP.σ)
    end
    y
end

function evaluate!(y::AbstractVector, CP::CrystalPhase, x::AbstractVector)
    peak_locs = (CP.cl).(CP.peaks) .* 10
    @fastmath @inbounds @simd for i in eachindex(CP.peaks)
        @. y += CP.act * CP.peaks[i].I * CP.profile((x-peak_locs[i])/CP.σ) # Main bottle neck
    end
    y
end

function evaluate!(y::AbstractVector, CPs::AbstractVector{<:AbstractPhase}, x::AbstractVector)
    for CP in CPs
        evaluate!(y, CP, x)
    end
    y
end

function evaluate!(y::AbstractVector, CP::CrystalPhase, θ::AbstractVector,
                   x::AbstractVector)
    CrystalPhase(CP, θ)(x, y)
end

function evaluate!(y::AbstractVector, CPs::AbstractVector{<:AbstractPhase},
                   θ::AbstractVector, x::AbstractVector)
    s = 1
    for i in eachindex(CPs)
        num_of_param = get_param_nums(CPs[i])
        θ_temp = @view θ[s : s+num_of_param-1]
        evaluate!(y, CPs[i], θ_temp, x)
        s += num_of_param
    end
    y
end

# function new_evaluate!(y::AbstractVector, CP::CrystalPhase, x::AbstractVector)
#     q = get_peak_pos.((CP.cl, ), CP.peaks) .* 10 # TODO: preallocate
#     temp = zero(q)
#     @. y += get_evalutation_at((CP, ), x, (q, ), (temp, )) # Main bottle neck
#     y
# end

# function get_evalutation_at(CP::CrystalPhase, x::Real, q::AbstractVector, temp::AbstractVector)
#     temp .= get_intensities(CP) .* CP.profile.(get_distance(CP, x, q))
#     return CP.act * sum(temp)
# end

function get_distance(CP::CrystalPhase, x::Real, q::AbstractVector)
    (x.-q)/CP.σ
end

function get_intensities(CP::CrystalPhase)
    [CP.peaks[i].I for i in eachindex(CP.peaks)]
end

# function new_evaluate!(y::AbstractVector, CPs::AbstractVector{<:CrystalPhase}, x::AbstractVector)
#     for CP in CPs
#         new_evaluate!(y, CP, x)
#     end
#     y
# end

function evaluate_residual!(CP::CrystalPhase, θ::AbstractVector,
              x::AbstractVector, r::AbstractVector)
    evaluate_residual!(CrystalPhase(CP, θ), x, r)
end

function evaluate_residual!(CPs::AbstractVector{<:AbstractPhase},
              x::AbstractVector, r::AbstractVector)
    @inbounds @simd for i in eachindex(CPs)
        evaluate_residual!(CPs[i], x, r)
    end
    r
end

function evaluate_residual!(CPs::AbstractVector{<:CrystalPhase},
    x::AbstractVector, r::AbstractVector)
    @inbounds @simd for i in eachindex(CPs)
        if isinf(CPs[i].cl.volume)
            return Inf
        end
        evaluate_residual!(CPs[i], x, r)
    end
    r
end

function evaluate_residual!(CP::CrystalPhase, x::AbstractVector, r::AbstractVector)
    # println(typeof(CP.cl))
    peak_locs = (CP.cl).(CP.peaks) .* 10
    # y = zero(r)
    @inbounds @simd for i in eachindex(CP.peaks)
        if isinf(peak_locs[i])
            return Inf
        end
        @. r -= CP.act * CP.peaks[i].I * CP.profile((x-peak_locs[i])/CP.σ)
        # Note: For Dual type, plus and minus are more expensive than multiply..
        # Slight improvement for 3 phase case
        # @. y += CP.peaks[i].I * CP.profile((x-peak_locs[i])/CP.σ)
    end
    # @. r -= y * CP.act
    r
end

function evaluate_residual!(CPs::AbstractVector{<:AbstractPhase},
             θ::AbstractVector, x::AbstractVector, r::AbstractVector)
    s = 1
    @inbounds @simd for i in eachindex(CPs)
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
    @inbounds @simd for i in eachindex(CP.peaks)
        q = (CP.cl)(CP.peaks[i]) * 10 # account for unit difference
        y += CP.act * CP.peaks[i].I * CP.profile((x-q)/CP.σ) # Main bottle neck
    end
    y
end

function (CPs::AbstractVector{<:AbstractPhase})(x::Real)
    y = zero(x)
    @inbounds @simd for i in eachindex(CPs) #
        y += CPs[i](x)
    end
    y
end


function evaluate(CP::CrystalPhase, θ::AbstractVector,
                      x::AbstractVector)
    CrystalPhase(CP, θ).(x)
end

function evaluate(CPs::AbstractVector{<:AbstractPhase},
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
