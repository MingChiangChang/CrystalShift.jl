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
    params = get_eight_params(CP.cl, θ)
    c = CrystalPhase(get_crystal(params[1:6], false), CP.peaks, CP.id, CP.name,
                     params[7], params[8], CP.profile)
    return c
end

get_eight_params(crystal::Cubic, θ::AbstractVector) = [θ[1], θ[1], θ[1], pi/2, pi/2, pi/2, θ[2], θ[3]]
get_eight_params(crystal::Tetragonal, θ::AbstractVector) = [θ[1], θ[1], θ[2], pi/2, pi/2, pi/2, θ[3], θ[4]]
get_eight_params(crystal::Orthorhombic, θ::AbstractVector) = [θ[1], θ[2], θ[3], pi/2, pi/2, pi/2, θ[4], θ[5]]
get_eight_params(crystal::Rhombohedral, θ::AbstractVector) = [θ[1], θ[1], θ[1], θ[2], θ[2], θ[2], θ[3], θ[4]]
get_eight_params(crystal::Hexagonal, θ::AbstractVector) = [θ[1], θ[1], θ[2], pi/2, pi/2, 2*pi/3, θ[3], θ[4]]
get_eight_params(crystal::Monoclinic, θ::AbstractVector) = [θ[1], θ[2], θ[3], pi/2, θ[4], pi/2, θ[4], θ[5]]
get_eight_params(crystal::Triclinic, θ::AbstractVector) = θ

function get_peaks(lines)
    peaks = Vector{Peak}(undef, size(lines))
    for i in eachindex(lines)
        peaks[i] = Peak(String(lines[i]))
    end
    peaks
end

function get_parameters(CP::CrystalPhase)
    return [get_free_params(CP.cl)..., CP.act, CP.σ]
end

function get_parameters(CPs::AbstractVector{<:CrystalPhase})
    p = Vector{Float64}()
    for cp in CPs
        push!(p, get_parameters(cp)...) # Preallocation?
    end
    p
end

# Preallocating
function (CP::CrystalPhase)(x::AbstractVector, y::AbstractVector)
    @simd for i in eachindex(CP.peaks)
        q = (CP.cl)(CP.peaks[i]) * 10 # account for unit difference
        @. y += CP.act * CP.peaks[i].I * CP.profile((x-q)/CP.σ) # Main bottle neck
    end
    y
end

function (CPs::AbstractVector{<:CrystalPhase})(x::AbstractVector,
                                               y::AbstractVector)
    @simd for i in eachindex(CPs)
        CPs[i](x, y)
    end
    y
end

function reconstruct!(CP::CrystalPhase, θ::AbstractVector,
                      x::AbstractVector, y::AbstractVector)
    CrystalPhase(CP, θ)(x, y)
end

function reconstruct!(CPs::AbstractVector{<:CrystalPhase},
    θ::AbstractVector, x::AbstractVector, y::AbstractVector)
    s = 1
    for i in eachindex(CPs)
        num_of_param = get_param_nums(CPs[i])
        θ_temp = @view θ[s:s+num_of_param-1]
        reconstruct!(CPs[i], θ_temp, x, y)
        s += num_of_param
    end
    y
end

function res!(CP::CrystalPhase, θ::AbstractVector,
              x::AbstractVector, r::AbstractVector)
    res!(CrystalPhase(CP, θ), x, r)
end

function res!(CPs::AbstractVector{<:CrystalPhase},
              x::AbstractVector, r::AbstractVector)
    @simd for i in eachindex(CPs)
        res!(CPs[i], x, r)
    end
    r
end

function res!(CP::CrystalPhase, x::AbstractVector, r::AbstractVector)
    @simd for i in eachindex(CP.peaks)
        q = (CP.cl)(CP.peaks[i]) * 10 # account for unit difference
        @. r -= CP.act * CP.peaks[i].I * CP.profile((x-q)/CP.σ) # Main bottle neck
    end
    r
end

function res!(CPs::AbstractVector{<:CrystalPhase},
             θ::AbstractVector, x::AbstractVector, r::AbstractVector)
    s = 1
    for i in eachindex(CPs)
        num_of_param = get_param_nums(CPs[i])
        θ_temp = @view θ[s:s+num_of_param-1]
        res!(CPs[i], θ_temp, x, r)
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
    @simd for i in eachindex(CPs)
        y += CPs[i](x)
    end
    y
end

function reconstruct!(CP::CrystalPhase, θ::AbstractVector,
                      x::AbstractVector)
    CrystalPhase(CP, θ).(x)
end

# Index
function reconstruct!(CPs::AbstractVector{<:CrystalPhase},
                      θ::AbstractVector, x::AbstractVector)
    y = zeros(size(x))
    s = 1
    for i in eachindex(CPs)
        num_of_param = get_param_nums(CPs[i])
        θ_temp = @view θ[s:s+num_of_param-1]
        y += reconstruct!(CPs[i], θ_temp, x)
        s += num_of_param
    end
    y
end
