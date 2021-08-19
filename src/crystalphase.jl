using PhaseMapping: Lorentz
using ForwardDiff: Dual

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

function Base.show(io::IO, CP::CrystalPhase)
    println("Phase name: $(CP.name), ID: $(CP.id)")
    println("Optimization parameters:")
    println("Activation: $(CP.act), Peak width: $(CP.σ)")
    println("Lattice information:")
    println(CP.cl)
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
    θ_temp = θ[1:get_param_nums(CP)]
    params = get_six_params(crystal, θ_temp)
    #println(isCubic(params...))
    #println(params, get_crystal(params, false))
    c = CrystalPhase(get_crystal(params, false), CP.peaks, CP.id, CP.name,
                θ_temp[end-1], θ_temp[end], CP.profile)
    return c
end

function get_six_params(crystal::Type, θ::AbstractVector)
    if crystal <: Cubic
        return [θ[1], θ[1], θ[1], pi/2, pi/2, pi/2]
    elseif crystal <: Tetragonal
        return [θ[1], θ[1], θ[2], pi/2, pi/2, pi/2]
    elseif crystal <: Orthorhombic
        return [θ[1], θ[2], θ[3], pi/2, pi/2, pi/2]
    elseif crystal <: Rhombohedral
        return [θ[1], θ[1], θ[1], θ[2], θ[2], θ[2]]
    elseif crystal <: Hexagonal
        return [θ[1], θ[1], θ[2], pi/2, pi/2, 2*pi/3]
    elseif crystal <: Monoclinic
        return [θ[1], θ[2], θ[3], pi/2, θ[4], pi/2]
    elseif crystal <: Triclinic
        return θ
    end
end

function get_peaks(lines)
    peaks = Vector{Peak}()
    for line in lines
        push!(peaks, Peak(String(line)))
    end
    peaks
end

function get_parameters(CP::CrystalPhase)
    return [get_free_params(CP.cl)..., CP.act, CP.σ]
end

function get_parameters(CPs::AbstractVector{<:CrystalPhase})
    p = Vector{Float64}()
    for cp in CPs
        push!(p, get_parameters(cp)...)
    end
    p
end

# Reconstruct spectrum
function (CP::CrystalPhase)(x::AbstractVector)
    y = zero(x)
    @simd for i in eachindex(CP.peaks)
        q = (CP.cl)(CP.peaks[i]) * 10 # account for unit difference
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
    #println(θ[1:num_of_param])
    y = CrystalPhase(CP, θ[1:num_of_param])(x)
    deleteat!(θ, collect(1:num_of_param))
    return y
end

function reconstruct!(CPs::AbstractVector{<:CrystalPhase},
                      θ::AbstractVector, x::AbstractVector)
    y = zeros(size(x))
    for i in eachindex(CPs)
        y += reconstruct!(CPs[i], θ, x)
    end
    #plt = plot(x, y)
    #savefig("recon$(sum(y)).png")
    y
end
