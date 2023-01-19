cast(f::AbstractVector, type::Type) = map(x->parse(type, x), f)

# "chess_x_y_tpeak_dwell"
function parse_cond(cond::AbstractString, T::Type)
    conds = split(cond, "_")[4:end]
    return cast(conds, T)
end

# TODO Should have a better way to identify center
function get_weighted_center(H::AbstractMatrix)
    avg = []
    for act in eachrow(H)
        ind = collect(1:size(act, 1))
        s = sum(act)
        if s>0
            push!(avg, Int(floor( sum(act.*ind)/sum(act) )))
        else
            print(size(act)[1]/2)
            push!(avg, Int(floor(size(act)[1]/2)))
        end
    end
    sum(avg)/length(avg)
end

# Helper functions for creating CrystalPhase object from CIFs
_get_key(cif_dict) = collect(keys(cif_dict))[1]
function _get_phase_name(info_dict)
    phase_name = ""
    space_group = ""
    if "_chemical_formula_structural" in keys(info_dict)
        phase_name = remove_blank!(info_dict["_chemical_formula_structural"])
    else
        phase_name = remove_blank!(info_dict["_chemical_formula_moiety"])
    end
    try
        space_group = remove_blank!(info_dict["_space_group_name_H-M_alt"])
    catch KeyError
        println("No _space_group_name_H-M_alt for $(phase_name)")
    end
    return phase_name * "_" * space_group
end

function _get_crystal_system(info_dict)
    try
        global sg_num = int(info_dict["_space_group_IT_number"])
    catch KeyError
        println("No space group info in cif. Default to triclinic")
        return "triclinic"
    end
    if sg_num in [1,2]
        return "triclinic"
    elseif 3 <= sg_num <= 15
        return "monoclinic"
    elseif 16 <= sg_num <= 74
        return "orthohombic"
    elseif 75 <= sg_num <= 142
        return "tetragonal"
    elseif 143 <= sg_num <= 167
        return "trigonal"
    elseif 168 <= sg_num <= 194
        return "hexagonal"
    elseif 195 <= sg_num <= 230
        return "cubic"
    end
end

function _get_lattice_parameters(info_dict)
    a, b, c = _get_cell_length(info_dict)
    alpha, beta, gamma = _get_cell_angle(info_dict)
    map(x->parse(Float64, x), [a, b, c, alpha, beta, gamma])
end

function _get_cell_length(info_dict)
    return (remove_parentheses(info_dict["_cell_length_a"]),
            remove_parentheses(info_dict["_cell_length_b"]),
            remove_parentheses(info_dict["_cell_length_c"]))
end

function _get_cell_angle(info_dict)
    return (remove_parentheses(info_dict["_cell_angle_alpha"]),
            remove_parentheses(info_dict["_cell_angle_beta"]),
            remove_parentheses(info_dict["_cell_angle_gamma"]))
end

remove_blank!(a::AbstractString) = replace(a, " "=>"")

function remove_parentheses(a::AbstractString)
    if occursin("(", a)
        return a[1:collect(findfirst("(", a))[1]-1]
    end
    a
end

## convenient macros
macro exported_enum(name, args...)
    esc(quote
        @enum($name, $(args...))
        export $name
        $([:(export $arg) for arg in args]...)
        end)
end

function volume(a::Real, b::Real, c::Real,
                α::Real, β::Real, γ::Real)
    try
        return (a * b * c
        * sqrt( 1+2*cos(α)*cos(β)*cos(γ)
                - cos(α)^2 - cos(β)^2 - cos(γ)^2 ) )
    catch DomainError
        return Inf
    end
end

function special_inv(x::Real)
    if x > 100.
        return 1E-12
    else
        return inv(x)
    end
end

function special_exp(x::Real)
    if x < -10.
        return 1E-12
    else
        return @fastmath exp(x)
    end
end



########################## objective functions #################################
# Kullback-Leibler divergence
function kl(p::Real, q::Real)
    !(isnan(p) || isinf(p) || isnan(q) || isinf(q)) || throw(DomainError("$(p) Nan or inf input"))
    (p > 0 && q > 0) || throw(DomainError("kl divergence undefined for negative inputs"))
    iszero(p) ? zero(p) : (p*log(p) - p*log(q)) # -p*log(q/p) # p*log(p/q) # is equal to
end

function kl(P::AbstractArray, Q::AbstractArray)
    length(P) == length(Q) || throw(DimensionMismatch("Lengths of inputs arrays do not match"))
    val = zero(promote_type(eltype(P), eltype(Q)))
    @inbounds @simd for i in eachindex(P)
        val += kl(P[i], Q[i])
    end
    return val
end


# negative Poisson likelihood, generalized for real (not just integer k)
function negative_log_poisson(λ::Real, k::Real)
    λ + loggamma(k) - k*log(λ)
end

# if we are optimizating w.r.t. λ, gamma(k) is irrelevant, so we remove it for this function
# NOTE: this is equivalent to optimizing kl w.r.t. 2nd argument and having l1 penalty
function negative_log_poisson_λ(λ::Real, k::Real)
    λ - k*log(λ)
end

poisson(λ::Real, k::Real) = λ^k * exp(-λ) / gamma(k)

var_lognormal(μ::Real, σ::Real) = (exp(σ^2) - 1) * exp(2μ+σ^2)
