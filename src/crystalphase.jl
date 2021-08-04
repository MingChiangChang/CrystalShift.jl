using PhaseMapping: Lorentz

struct CrystalPhase{T, V<:AbstractVector{T}, C, P, K, M}
    cl::C # crystal object
    peaks::V # Vector of peak object

    id::Int # Just index
    name::String # For printing

    act::K # Activation
    σ::M # Width of peaks
    profile::P # Peak profile
               # Do import from PhaseMapping if
               # other peak profiles are needed
end

get_param_nums(CP::CrystalPhase) = CP.cl.free_param + 2

function CrystalPhase(_stn::String, wid_init::Real=.1,
                      profile=Lorentz())
    f = split(_stn, '\n')
    lattice_info = split(f[1], ',')
    id = parse(Int64, lattice_info[1])
    crystal = get_crystal(cast(lattice_info[4:end], Float64))
    peaks = get_peaks(f[2:end])
    name = String(lattice_info[2])
    act = 1.0

    CrystalPhase(crystal, peaks, id, name, act, wid_init, profile)
end

# Create a new CP object with the new parameters
function CrystalPhase(CP::CrystalPhase, θ::AbstractVector)
    crystal = typeof(CP.cl)
    CrystalPhase(crystal(θ[1:end-2]...), CP.peaks, CP.id, CP.name,
                θ[end-1], θ[end], CP.profile)
end

# TODO ... unpacking does not work for arrays?
function get_peaks(lines)
    peaks = Vector{Peak}()
    for line in lines
        push!(peaks, Peak(String(line)))
    end
    peaks
end

# Reconstruct spectrum
function (CP::CrystalPhase)(x::AbstractVector)
    y = zero(x)
    @simd for i in eachindex(CP.peaks)
        q = (CP.cl)(CP.peaks[i])*100 # TODO fix python script to make input to be nm
        println(q)
        y += CP.act * CP.peaks[i].I * CP.profile.((x.-q)/CP.σ)
    end
    y
end

function (CPs::AbstractVector{<:CrystalPhase})(x::AbstractVector)
    y = zero(x)
    @simd for i in eachindex(CPs)
        y += CPs[i](x)
    end
    y
end

function reconstruct!(CP::CrystalPhase, θ::AbstractVector, x::AbstractVector)
    # Pop the first (# free param of CP) and create a new phase
    # reconstruction
    num_of_param = get_param_nums(CP)
    y = CrystalPhase(CP, θ)(x)
    deleteat!(θ, collect(1:num_of_param))
    return y
end

function reconstruct!(CPs::AbstractVector{<:CrystalPhase},
                      θ::AbstractVector, x::AbstractVector)
    y = zeros(size(x))
    for i in eachindex(CPs)
        y += reconstruct(CPs[i], θ, x)
    end
end
