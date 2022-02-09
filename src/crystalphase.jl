struct CrystalPhase{T, V<:AbstractVector{T}, C, CL, P, K, M, N}
    cl::C # crystal object
    origin_cl::CL # save for later comparison
    peaks::V # Vector of peak object

    id::Int # Just index
    name::String # For printing

    act::K # Activation
    σ::M # Width of peaks
    profile::P # Peak profile
               # Do import from PhaseMapping if
               # other peak profiles are needed
    norm_constant::N
end

function Base.show(io::IO, CP::CrystalPhase)
    println("Phase name: $(CP.name), ID: $(CP.id)")
    println("Optimization parameters:")
    println("Activation: $(CP.act), Peak width: $(CP.σ)")
    println("Normalization: $(CP.norm_constant)")
    println("Lattice information:")
    println(CP.cl)
end

Base.Bool(CP::CrystalPhase) = true
Base.Bool(CPs::AbstractVector{<:CrystalPhase}) = true
get_param_nums(CP::CrystalPhase) = CP.cl.free_param + 2 + get_param_nums(CP.profile)
get_param_nums(CPs::AbstractVector{<:CrystalPhase}) = sum(get_param_nums.(CPs))

function CrystalPhase(_stn::String, wid_init::Real=.1,
                      profile::PeakProfile=PseudoVoigt(0.5))
    f = split(_stn, '\n')
    lattice_info = split(f[1], ',')
    id = parse(Int64, lattice_info[1])
    crystal = get_crystal(cast(lattice_info[4:end], Float64))
    peaks, norm_constant = get_peaks(f[2:end])
    name = String(lattice_info[2])
    act = 1.0

    CrystalPhase(crystal, crystal, peaks, id, name, act, wid_init, profile, norm_constant)
end

function CrystalPhase(CP::CrystalPhase, θ::AbstractVector)
    fp = CP.cl.free_param
    cl = get_intrinsic_crystal_type(typeof(CP.cl))
    t = eltype(θ)
    c = CrystalPhase(cl{t}(θ[1:fp]...), CP.origin_cl, CP.peaks, CP.id, CP.name,
                     θ[fp+1], θ[fp+2], CP.profile, CP.norm_constant)
    return c
end

function reconstruct_CPs!(θ::AbstractVector, CPs::AbstractVector{<:CrystalPhase})
    start = 1
    new_CPs = Vector{CrystalPhase}(undef, length(CPs))
    for i in eachindex(CPs)
		new_CPs[i] = CrystalPhase(CPs[i], θ[start:start + get_param_nums(CPs[i])-1])
		start += get_param_nums(CPs[i])
	end
    θ[start:end], new_CPs
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

function get_free_params(CP::CrystalPhase)
    return vcat(get_free_lattice_params(CP.cl), [CP.act, CP.σ], get_free_params(CP.profile))
end

function get_free_params(CPs::AbstractVector{<:CrystalPhase})
    p = Vector{Float64}()
    for cp in CPs
        push!(p, get_free_params(cp)...) # Preallocation?
    end
    p
end

function get_eight_params(CP::CrystalPhase)
    vcat([CP.cl.a, CP.cl.b, CP.cl.c, CP.cl.α, CP.cl.β, CP.cl.γ, CP.act, CP.σ], get_free_params(CP.profile))
end

get_eight_params(CP::CrystalPhase, θ::AbstractVector) = get_eight_parmams(CrystalPhase(CP, θ))
get_eight_params(crystal::Cubic, θ::AbstractVector) = [θ[1], θ[1], θ[1], pi/2, pi/2, pi/2, θ[2], θ[3]]
get_eight_params(crystal::Tetragonal, θ::AbstractVector) = [θ[1], θ[1], θ[2], pi/2, pi/2, pi/2, θ[3], θ[4]]
get_eight_params(crystal::Orthorhombic, θ::AbstractVector) = [θ[1], θ[2], θ[3], pi/2, pi/2, pi/2, θ[4], θ[5]]
get_eight_params(crystal::Rhombohedral, θ::AbstractVector) = [θ[1], θ[1], θ[1], θ[2], θ[2], θ[2], θ[3], θ[4]]
get_eight_params(crystal::Hexagonal, θ::AbstractVector) = [θ[1], θ[1], θ[2], pi/2, pi/2, 2*pi/3, θ[3], θ[4]]
get_eight_params(crystal::Monoclinic, θ::AbstractVector) = [θ[1], θ[2], θ[3], pi/2, θ[4], pi/2, θ[5], θ[6]]
get_eight_params(crystal::Triclinic, θ::AbstractVector) = θ

function get_free_lattice_params(CPs::AbstractVector{<:CrystalPhase})
    p = Vector{Float64}()
    for cp in CPs
        push!(p, get_free_lattice_params(cp)...) # Preallocation?
    end
    p
end
get_free_lattice_params(CP::CrystalPhase) = get_free_lattice_params(CP.cl)
get_free_lattice_params(CP::CrystalPhase, θ::AbstractVector) = get_free_lattice_params(CP.cl, θ)
get_free_lattice_params(cl::Cubic, θ::AbstractVector) = [θ[1]]
get_free_lattice_params(cl::Tetragonal, θ::AbstractVector) = [θ[1], θ[3]]
get_free_lattice_params(cl::Hexagonal, θ::AbstractVector) = [θ[1], θ[3]]
get_free_lattice_params(cl::Orthorhombic, θ::AbstractVector) = [θ[1], θ[2], θ[3]]
get_free_lattice_params(cl::Rhombohedral, θ::AbstractVector) = [θ[1], θ[4]]
get_free_lattice_params(cl::Monoclinic, θ::AbstractVector) = [θ[1], θ[2], θ[3], θ[5]]
get_free_lattice_params(cl::Triclinic, θ::AbstractVector) = θ

function get_eight_params(CP::CrystalPhase, θ::AbstractVector)
    get_eight_params(CP.cl, θ)
end

collect_crystals(CPs::AbstractVector{<:CrystalPhase}) = [CP.cl for CP in CPs]

# Preallocating
# Functor comes in handy but use evaluate! when you can to be clear
function (CP::CrystalPhase)(x::AbstractVector, y::AbstractVector)
    evaluate!(y, CP, x)
end

function (CPs::AbstractVector{<:CrystalPhase})(x::AbstractVector,
                                               y::AbstractVector)
    @simd for i in eachindex(CPs)
        CPs[i](x, y)
    end
    y
end

function evaluate!(y::AbstractVector, CP::CrystalPhase, x::AbstractVector)
    @simd for i in eachindex(CP.peaks)
        q = (CP.cl)(CP.peaks[i]) * 10 # account for unit difference
        @. y += CP.act * CP.peaks[i].I * CP.profile((x-q)/CP.σ) # Main bottle neck
    end
    y
end

function evaluate!(y::AbstractVector, CPs::AbstractVector{<:CrystalPhase}, x::AbstractVector)
    for CP in CPs
        evaluate!(y, CP, x)
    end
    y
end

function evaluate!(CP::CrystalPhase, θ::AbstractVector,
                   x::AbstractVector, y::AbstractVector)
    CrystalPhase(CP, θ)(x, y)
end

function evaluate!(y::AbstractVector, CPs::AbstractVector{<:Crystal},
                   θ::AbstractVector, x::AbstractVector)
    s = 1
    for i in eachindex(CPs)
        num_of_param = get_param_nums(CPs[i])
        θ_temp = @view θ[s : s+num_of_param-1]
        evalute!(y, CPs[i], θ_temp, x)
        s += num_of_param
    end
    y
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
