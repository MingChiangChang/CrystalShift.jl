# Hand-written Jacobian for the CrystalPhase blocks of BlockJacobianLM (see lmsolver.jl).
#
# A phase block maps ℓ = log.([free lattice params; act; σ; profile params]) to
#   out[1:n]     = -w .* model,  model_j = Σ_i c_i (φ(u_ij) - F)·[j in window i] + (Σ_i c_i) F
#   out[n+1:end] = (ℓ - mean_log_θ) ./ (√2 std_θ)                   (if regularization)
# with c_i = act·I_i, u_ij = (x_j - q_i)/σ, q_i the peak positions and F = φ(2k) the profile
# floor, exactly as add_peaks! computes it. Writing the derivatives out evaluates each profile
# once per grid point instead of carrying one Dual partial per parameter through it:
#   ∂out/∂ℓ_act     = -w .* model
#   ∂out/∂ℓ_σ       =  w .* Σ_i c_i u_ij φ'(u_ij)
#   ∂out/∂ℓ_lattice =  w ./ σ .* Σ_i c_i φ'(u_ij) ∂q_i/∂ℓ_lattice
#   ∂out/∂ℓ_a       = -w .* a .* Σ_i c_i ∂φ/∂a(u_ij)            (PseudoVoigt mixing a)
# F is the same constant (1e-12) for every supported profile, so it has no derivative.
#
# Peak positions: every cl(peak) is 2π/V·√(quadratic form in h, k, l), so
#   q² = A h² + B k² + C l² + 2D hk + 2E kl + 2F hl.
# The 6 coefficients come from cl at 6 basis reflections, differentiated with ForwardDiff,
# and each peak is a fixed row [h² k² l² 2hk 2kl 2hl] of M: q² = M·coef and
# ∂q/∂ℓ = (M·∂coef/∂ℓ) / 2q, in plain Float64 however many peaks the phase has.
#
# Anything else (other profiles or phase types, unsorted x, an impossible cell) uses the
# block's ForwardDiff path. Set USE_ANALYTIC_JACOBIAN[] = false to always use ForwardDiff.

const USE_ANALYTIC_JACOBIAN = Ref(true)

const AnalyticProfile = Union{LorentzianProfile, GaussianProfile, PseudoVoigtProfile, FixedPseudoVoigtProfile}

# (φ, ∂φ/∂u) with the same floors as special_inv and special_exp
@inline function lorentz_derivative(u::Real)
    t = 1 + u^2
    t > 100 ? (1e-12, zero(u)) : (inv(t), -2u * inv(t)^2)
end
@inline function gauss_derivative(u::Real)
    t = -u^2 / 2
    t < -10 ? (1e-12, zero(u)) : (g = exp(t); (g, -u * g))
end

# (φ, ∂φ/∂u, ∂φ/∂a) where a is the profile's own parameter (zero if it has none)
@inline function profile_derivatives(::LorentzianProfile, u::Real)
    φ, d = lorentz_derivative(u)
    (φ, d, zero(u))
end
@inline function profile_derivatives(::GaussianProfile, u::Real)
    φ, d = gauss_derivative(u)
    (φ, d, zero(u))
end
@inline function profile_derivatives(P::FixedPseudoVoigtProfile, u::Real)
    (l, dl), (g, dg) = lorentz_derivative(u), gauss_derivative(u)
    (P.α * l + (1 - P.α) * g, P.α * dl + (1 - P.α) * dg, zero(u))
end
@inline function profile_derivatives(P::PseudoVoigtProfile, u::Real)
    (l, dl), (g, dg) = lorentz_derivative(u), gauss_derivative(u)
    s = P.sig_α
    ((s - 0.5) * l + (1.5 - s) * g, (s - 0.5) * dl + (1.5 - s) * dg, s * (1 - s) * (l - g))
end

################################ Peak positions ################################
const BASIS_HKL = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (0, 1, 1), (1, 0, 1))

struct PeakPositions{F, C}
    coef!::F               # coef!(coef, ℓ_lattice) writes (A, B, C, D, E, F)
    cfg::C
    coef::Vector{Float64}
    Dcoef::Matrix{Float64} # ∂coef/∂ℓ_lattice
    M::Matrix{Float64}     # one row [h² k² l² 2hk 2kl 2hl] per peak
    q::Vector{Float64}
    Dq::Matrix{Float64}    # ∂q/∂ℓ_lattice
end

function PeakPositions(cl::Crystal, peaks::AbstractVector{<:Peak})
    fp = Int(cl.free_param)
    cltype = get_intrinsic_crystal_type(typeof(cl))
    basis = [Peak(Int16(h), Int16(k), Int16(l), 0.0, 1.0) for (h, k, l) in BASIS_HKL]
    coef! = (coef, ℓ) -> begin
        c = cltype{eltype(ℓ)}(exp.(ℓ)...)
        s = ntuple(i -> (10 * c(basis[i]))^2, 6) # 10: same unit factor as add_peaks!
        coef[1], coef[2], coef[3] = s[1], s[2], s[3]
        coef[4] = (s[4] - s[1] - s[2]) / 2
        coef[5] = (s[5] - s[2] - s[3]) / 2
        coef[6] = (s[6] - s[1] - s[3]) / 2
        coef
    end
    M = zeros(length(peaks), 6)
    for (i, p) in enumerate(peaks)
        h, k, l = Float64(p.h), Float64(p.k), Float64(p.l) # Float64 first: h, k, l are Int16
        M[i, :] .= (h^2, k^2, l^2, 2h*k, 2k*l, 2h*l)
    end
    coef = zeros(6)
    PeakPositions(coef!, ForwardDiff.JacobianConfig(coef!, coef, zeros(fp)), coef, zeros(6, fp),
                  M, zeros(length(peaks)), zeros(length(peaks), fp))
end

# Fills P.q and P.Dq at ℓ_lattice; false if the cell is impossible (some q² ≤ 0)
function peak_positions!(P::PeakPositions, ℓ::AbstractVector)
    ForwardDiff.jacobian!(P.Dcoef, P.coef!, P.coef, ℓ, P.cfg)
    all(isfinite, P.coef) || return false
    mul!(P.q, P.M, P.coef)    # q²
    all(>(0), P.q) || return false
    @. P.q = sqrt(P.q)
    mul!(P.Dq, P.M, P.Dcoef)  # ∂q²/∂ℓ
    P.Dq ./= 2 .* P.q
    true
end

################################ Phase blocks ################################
# Jacobian block of one CrystalPhase. cols, out and J are the arrays of the ForwardDiff block
# `ad`, so add_block! and the fallback write to the same place.
struct PhaseJacobianBlock{A<:ADBlock, CPT<:CrystalPhase, X, W, P<:PeakPositions} <: JacobianBlock
    ad::A
    cols::UnitRange{Int}
    out::Vector{Float64}
    J::Matrix{Float64}
    CP::CPT
    x::X
    w::W
    n::Int
    reg::Bool
    mean_log_θ::Vector{Float64}
    std_θ::Vector{Float64}
    pos::P
    t::Vector{Float64} # per-grid-point buffer for c·φ'
end

# The analytic block when it applies to CP and x, else the ForwardDiff block itself
function phase_block(ad::ADBlock, CP::AbstractPhase, x::AbstractVector, w::AbstractVector,
                     n::Int, reg::Bool, mean_log_θ::AbstractVector, std_θ::AbstractVector)
    CP isa CrystalPhase && CP.profile isa AnalyticProfile && eltype(x) == Float64 && issorted(x) || return ad
    PhaseJacobianBlock(ad, ad.cols, ad.out, ad.J, CP, x, w, n, reg, Vector{Float64}(mean_log_θ),
                       Vector{Float64}(std_θ), PeakPositions(CP.cl, CP.peaks), zeros(n))
end

block_jacobian!(b::ADBlock, θ::AbstractVector) = ForwardDiff.jacobian!(b.J, b.g!, b.out, θ, b.cfg)

function block_jacobian!(b::PhaseJacobianBlock, θ::AbstractVector)
    USE_ANALYTIC_JACOBIAN[] && analytic_jacobian!(b, θ) && return b.J
    block_jacobian!(b.ad, θ)
end

# Writes b.out and b.J at θ (log space) like the ForwardDiff block would; false (writing
# nothing usable) when the peak positions are invalid, so the caller falls back.
function analytic_jacobian!(b::PhaseJacobianBlock, θ::AbstractVector)
    CP0, x, w, n, J, out, t = b.CP, b.x, b.w, b.n, b.J, b.out, b.t
    fp = Int(CP0.cl.free_param)
    peak_positions!(b.pos, view(θ, 1:fp)) || return false
    q, Dq = b.pos.q, b.pos.Dq
    CP = CrystalPhase(CP0, exp.(θ))
    act, σ, P = CP.act, CP.σ, CP.profile
    ia, iσ, ip = fp + 1, fp + 2, fp + 3
    hasprof = size(J, 2) > fp + 2
    k = support(P)
    F = P(2k)
    half = 1.01k * abs(σ) # same window as add_peaks!
    model = view(out, 1:n)
    model .= 0
    J .= 0
    total = zero(act)
    @inbounds for i in eachindex(CP.peaks)
        c = act * CP.peaks[i].I
        total += c
        qi = q[i]
        lo = searchsortedfirst(x, qi - half)
        hi = searchsortedlast(x, qi + half)
        lo > hi && continue
        @simd for j in lo:hi
            u = (x[j] - qi) / σ
            φ, dφ, dφa = profile_derivatives(P, u)
            model[j] += c * (φ - F)
            cd = c * dφ
            t[j] = cd
            J[j, iσ] += cd * u
            hasprof && (J[j, ip] += c * dφa)
        end
        for l in 1:fp # column-contiguous: J[lo:hi, l] += Dq[i, l] .* t[lo:hi]
            d = Dq[i, l]
            @simd for j in lo:hi
                J[j, l] += d * t[j]
            end
        end
    end
    @inbounds for l in 1:fp
        @simd for j in 1:n
            J[j, l] *= w[j] / σ
        end
    end
    @inbounds @simd for j in 1:n
        mj = model[j] + total * F
        J[j, ia] = -w[j] * mj
        J[j, iσ] *= w[j]
        out[j] = -w[j] * mj
    end
    if hasprof
        a = P.α
        @inbounds @simd for j in 1:n
            J[j, ip] *= -w[j] * a
        end
    end
    if b.reg
        @inbounds for l in eachindex(θ)
            s = sqrt(2) * b.std_θ[l]
            out[n+l] = (θ[l] - b.mean_log_θ[l]) / s
            J[n+l, l] = 1 / s
        end
    end
    true
end

############################ Second-order term ############################
# For the exact Hessian of Σⱼ rⱼ² (get_probabilities): H = 2(JᵀJ + Σⱼ rⱼ ∇²rⱼ). Phases add, so
# ∇²rⱼ has no cross-phase entries and Σⱼ rⱼ ∇²rⱼ is block diagonal; residual_curvature!
# writes one phase's block S = Σⱼ rⱼ ∂²rⱼ/∂ℓ² from a few sums per peak window, with
# r̃ = r .* w, the same u, φ windows as analytic_jacobian! and, for one peak i with c = act·Iᵢ,
#   lattice m, n: c (∂q/∂ℓₘ ∂q/∂ℓₙ Σr̃φ''/σ² - ∂²q/∂ℓₘ∂ℓₙ Σr̃φ'/σ)
#   lattice m, σ: c ∂q/∂ℓₘ (Σr̃uφ'' + Σr̃φ')/σ
#   σ, σ:         c (Σr̃u²φ'' + Σr̃uφ')
#   with a PseudoVoigt mixing a = exp(ℓ_a): lattice m: -c a ∂q/∂ℓₘ Σr̃∂²φ/∂u∂a / σ,
#   σ: -c a Σr̃u∂²φ/∂u∂a, a: c (a Σr̃∂φ/∂a + a² Σr̃∂²φ/∂a²)
# (all with a minus sign in S, as rⱼ = wⱼ(yⱼ - model)). The model is proportional to act, so
# the act row is Σⱼ rⱼ J[j, :].

# (φ, φ', φ'') with the same floors as special_inv and special_exp
@inline function lorentz_second_derivative(u::Real)
    t = 1 + u^2
    t > 100 ? (1e-12, zero(u), zero(u)) : (s = inv(t); (s, -2u * s^2, (6u^2 - 2) * s^3))
end
@inline function gauss_second_derivative(u::Real)
    t = -u^2 / 2
    t < -10 ? (1e-12, zero(u), zero(u)) : (g = exp(t); (g, -u * g, (u^2 - 1) * g))
end

# (∂φ/∂u, ∂²φ/∂u², ∂φ/∂a, ∂²φ/∂u∂a, ∂²φ/∂a²), a the profile's own parameter
@inline function profile_second_derivatives(::LorentzianProfile, u::Real)
    _, d, dd = lorentz_second_derivative(u)
    (d, dd, zero(u), zero(u), zero(u))
end
@inline function profile_second_derivatives(::GaussianProfile, u::Real)
    _, d, dd = gauss_second_derivative(u)
    (d, dd, zero(u), zero(u), zero(u))
end
@inline function profile_second_derivatives(P::FixedPseudoVoigtProfile, u::Real)
    (_, dl, ddl), (_, dg, ddg) = lorentz_second_derivative(u), gauss_second_derivative(u)
    (P.α * dl + (1 - P.α) * dg, P.α * ddl + (1 - P.α) * ddg, zero(u), zero(u), zero(u))
end
@inline function profile_second_derivatives(P::PseudoVoigtProfile, u::Real)
    (l, dl, ddl), (g, dg, ddg) = lorentz_second_derivative(u), gauss_second_derivative(u)
    s = P.sig_α
    ds = s * (1 - s) # d sigmoid / da
    ((s - 0.5) * dl + (1.5 - s) * dg, (s - 0.5) * ddl + (1.5 - s) * ddg,
     ds * (l - g), ds * (dl - dg), ds * (1 - 2s) * (l - g))
end

# ∂²q/∂ℓ∂ℓ of every peak (npeaks × fp × fp), from ∂²coef/∂ℓ∂ℓ at the 6 basis reflections;
# call after peak_positions!(P, ℓ), which fills P.q and P.Dq
function peak_position_hessians(P::PeakPositions, ℓ::AbstractVector)
    fp = length(ℓ)
    coef(t) = P.coef!(zeros(eltype(t), 6), t)
    D2coef = ForwardDiff.jacobian(t -> ForwardDiff.jacobian(coef, t), ℓ) # (6·fp) × fp
    position_hessians(P.M, Matrix{Float64}(reshape(D2coef, 6, fp * fp)), P.q, P.Dq)
end

# ∂²q/∂ℓ∂ℓ from q² = M·coef: (M·∂²coef)/2q - ∂q ∂qᵀ/q, as npeaks × fp × fp
function position_hessians(M::Matrix{Float64}, D2coef::Matrix{Float64}, q::Vector{Float64}, Dq::Matrix{Float64})
    fp = size(Dq, 2)
    D2q = reshape(M * D2coef, length(q), fp, fp) # ∂²q²/∂ℓₘ∂ℓₙ
    @inbounds for n in 1:fp, m in 1:fp, i in eachindex(q)
        D2q[i, m, n] = D2q[i, m, n] / (2q[i]) - Dq[i, m] * Dq[i, n] / q[i]
    end
    D2q
end

# Writes S = Σⱼ rⱼ ∂²rⱼ/∂θ² of the phase block b at θ (log space) for the full weighted
# residual r (data rows). Uses b.J, so call block_jacobian!(b, θ) first. Falls back to a
# ForwardDiff Hessian over this block's parameters only when the analytic path does not apply.
function residual_curvature!(S::AbstractMatrix, b::JacobianBlock, θ::AbstractVector, r::AbstractVector)
    g! = b isa PhaseJacobianBlock ? b.ad.g! : b.g!
    n = length(r)
    S .= ForwardDiff.hessian(t -> (o = zeros(eltype(t), length(b.out)); g!(o, t); dot(r, view(o, 1:n))), θ)
end

function residual_curvature!(S::AbstractMatrix, b::PhaseJacobianBlock, θ::AbstractVector, r::AbstractVector)
    USE_ANALYTIC_JACOBIAN[] || return invoke(residual_curvature!, Tuple{AbstractMatrix, JacobianBlock, AbstractVector, AbstractVector}, S, b, θ, r)
    CP0, x, w, n, J = b.CP, b.x, b.w, b.n, b.J
    fp = Int(CP0.cl.free_param)
    ℓ = θ[1:fp]
    peak_positions!(b.pos, ℓ) ||
        return invoke(residual_curvature!, Tuple{AbstractMatrix, JacobianBlock, AbstractVector, AbstractVector}, S, b, θ, r)
    q, Dq = b.pos.q, b.pos.Dq
    D2q = peak_position_hessians(b.pos, ℓ)
    CP = CrystalPhase(CP0, exp.(θ))
    act, σ, P = CP.act, CP.σ, CP.profile
    ia, iσ, ip = fp + 1, fp + 2, fp + 3
    hasprof = size(S, 2) > fp + 2
    a = hasprof ? P.α : 0.0
    half = 1.01support(P) * abs(σ) # same window as add_peaks!
    D = zeros(size(S)) # Σⱼ r̃ⱼ ∂²model/∂θ² without the act row
    @inbounds for i in eachindex(CP.peaks)
        c = act * CP.peaks[i].I
        qi = q[i]
        lo = searchsortedfirst(x, qi - half)
        hi = searchsortedlast(x, qi + half)
        lo > hi && continue
        A1 = A2 = B1 = B2 = B3 = C0 = C1 = C2 = C3 = 0.0
        @simd for j in lo:hi
            u = (x[j] - qi) / σ
            d, dd, da, dua, daa = profile_second_derivatives(P, u)
            rt = r[j] * w[j]
            A1 += rt * d; A2 += rt * dd
            B1 += rt * u * d; B2 += rt * u * dd; B3 += rt * u^2 * dd
            C0 += rt * da; C1 += rt * dua; C2 += rt * u * dua; C3 += rt * daa
        end
        for m in 1:fp
            gm = Dq[i, m]
            for l in 1:m
                D[m, l] += c * (gm * Dq[i, l] * A2 / σ^2 - D2q[i, m, l] * A1 / σ)
            end
            D[iσ, m] += c * gm * (B2 + A1) / σ
            hasprof && (D[ip, m] -= c * a * gm * C1 / σ)
        end
        D[iσ, iσ] += c * (B3 + B1)
        if hasprof
            D[ip, iσ] -= c * a * C2
            D[ip, ip] += c * (a * C0 + a^2 * C3)
        end
    end
    @inbounds for m in axes(S, 1), l in 1:m
        S[m, l] = S[l, m] = -D[m, l]
    end
    @inbounds for l in axes(S, 2) # act row and column: Σⱼ rⱼ ∂rⱼ/∂θ
        v = zero(eltype(S))
        @simd for j in 1:n
            v += r[j] * J[j, l]
        end
        S[ia, l] = S[l, ia] = v
    end
    S
end
