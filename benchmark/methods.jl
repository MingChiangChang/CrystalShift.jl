# Experimental optimization strategies, built on the wrappers in solvers.jl.
# Every strategy returns a full log-space θ for the original problem, so it is scored by
# exactly the same cost and pass criteria as the plain solvers.
using CrystalShift: reconstruct_BG!
using Random

const DOGLEG = lso(LSO.Dogleg())

# --- index helpers ----------------------------------------------------------------------

# Positions in θ of each phase's lattice parameters and activation
function phase_indices(p::Problem)
    lattice, act = Int[], Int[]
    p.pm.CPs isa AbstractVector{<:CrystalPhase} || return lattice, act
    off = 0
    for cp in p.pm.CPs
        fp = cp.cl.free_param
        append!(lattice, off+1:off+fp)
        push!(act, off + fp + 1)
        off += get_param_nums(cp)
    end
    lattice, act
end
bg_indices(p::Problem) = (n = get_param_nums(p.pm) ; (n - get_param_nums(p.pm.background) + 1):n)

with_y(p::Problem, y) = Problem(p.pm, p.x, y, p.stn)

# --- multi-start ------------------------------------------------------------------------
# Run Dogleg from the given start plus n-1 starts with lattice parameters jittered by up to
# ±δ (relative), keep the lowest cost. Deterministic per problem.
function multistart(n; δ = 0.0125, local_solver = DOGLEG)
    function (p, θ0)
        lattice, _ = phase_indices(p)
        rng = Xoshiro(hash(p.y))
        best, best_cost, total_iters = θ0, Inf, 0
        for k in 1:n
            θ = copy(θ0)
            k > 1 && (θ[lattice] .+= log1p.(δ .* (2 .* rand(rng, length(lattice)) .- 1)))
            θ, it = local_solver(p, θ)
            total_iters += it
            c = ls_cost(p, θ)
            c < best_cost && ((best, best_cost) = (θ, c))
        end
        best, total_iters
    end
end

# --- continuation -----------------------------------------------------------------------
# Fit progressively less-smoothed versions of the data, then the data itself. Smoothing
# widens the basin around each peak so lattice parameters can't lock onto the wrong one
# early; peak width σ is free and absorbs the extra broadening.
function gaussian_smooth(x, y, s)
    dx = x[2] - x[1]
    h = ceil(Int, 4s / dx)
    w = exp.(-0.5 .* ((-h:h) .* dx ./ s) .^ 2)
    out = similar(y)
    for i in eachindex(y)
        lo, hi = max(1, i - h), min(length(y), i + h)
        ww = @view w[(lo-i+h+1):(hi-i+h+1)]
        out[i] = sum(ww .* @view(y[lo:hi])) / sum(ww)
    end
    out
end

function continuation(widths = (0.4, 0.2); local_solver = DOGLEG)
    function (p, θ)
        iters = 0
        for s in widths
            θ, it = local_solver(with_y(p, gaussian_smooth(p.x, p.y, s)), θ)
            iters += it
        end
        θ, it = local_solver(p, θ)
        θ, iters + it
    end
end

# --- variable projection ----------------------------------------------------------------
# Activations and background coefficients enter the model linearly. For fixed nonlinear
# parameters φ (lattice, width, profile, wildcard), solve for them by linear least squares
# (data residual + the background's ridge prior) and let the outer solver see only φ.
# The outer objective is the full LM residual at (φ, linear solution), so priors on the
# activations still count, they're just not used in the inner solve. Activations are
# clamped positive (they live in log space). A final Dogleg pass on the full problem
# removes the small bias from leaving the activation prior out of the inner solve.
function linear_columns(p::Problem)
    # background columns don't depend on φ: precompute once
    nb = get_param_nums(p.pm.background)
    B = zeros(length(p.x), nb)
    for k in 1:nb
        e = zeros(nb); e[k] = 1
        _, bg = reconstruct_BG!(e, p.pm.background)
        evaluate!(view(B, :, k), bg, p.x)
    end
    D = zeros(nb, nb)  # ridge prior rows, from lm_prior!
    for k in 1:nb
        e = zeros(nb); e[k] = 1
        col = zeros(nb)
        CrystalShift.lm_prior!(col, p.pm.background, e)
        D[:, k] = col
    end
    B, D
end

function varpro(; polish = DOGLEG, outer = DOGLEG)
    function (p, θ0)
        lattice, act = phase_indices(p)
        bg = collect(bg_indices(p))
        lin = vcat(act, bg)
        isempty(lin) && return outer(p, θ0)
        nonlin = setdiff(eachindex(θ0), lin)
        nlogp = nlog(p)
        B, D = linear_columns(p)
        s = sqrt(2) * p.stn.priors.std_noise
        f_full = lm_residual(p)
        cps = p.pm.CPs isa AbstractVector{<:CrystalPhase} ? p.pm.CPs : CrystalPhase[]

        # assemble full log θ from φ, solving for the linear parameters
        function assemble(φ::AbstractVector{T}) where T
            θ = Vector{T}(undef, length(θ0))
            θ[nonlin] .= φ
            θ[act] .= 0  # act = 1 while building columns
            θr = copy(θ)
            θr[1:nlogp] .= exp.(θ[1:nlogp])
            A = zeros(T, length(p.x), length(act))
            off = 0
            for (j, cp) in enumerate(cps)
                m = get_param_nums(cp)
                evaluate!(view(A, :, j), cp, θr[off+1:off+m], p.x)
                off += m
            end
            # data minus the non-linear, non-background parts (wildcard) at current φ
            rest = T.(p.y)
            wn = get_param_nums(p.pm.wildcard)
            if wn > 0
                w0 = get_param_nums(p.pm.CPs)
                wpm = PhaseModel(nothing, p.pm.wildcard, nothing)
                rest .-= evaluate!(zeros(T, length(p.x)), wpm, θr[w0+1:w0+wn], p.x)
            end
            M = [hcat(A, B) ./ s; hcat(zeros(size(D, 1), length(act)), D)]
            rhs = [rest ./ s; zeros(size(D, 1))]
            sol = M \ rhs
            θ[act] .= log.(max.(sol[1:length(act)], 1e-8))
            θ[bg] .= sol[length(act)+1:end]
            θ
        end

        function f_outer!(r, φ)
            f_full(r, assemble(φ))
        end
        prob = LSO.LeastSquaresProblem(x = θ0[nonlin], f! = f_outer!, output_length = nresid(p), autodiff = :forward)
        res = LSO.optimize!(prob, LSO.Dogleg(); iterations = p.stn.maxiter)
        θ = assemble(res.minimizer)
        polish === nothing && return θ, res.iterations
        θ, it = polish(p, θ)
        θ, res.iterations + it
    end
end

# --- solver sets --------------------------------------------------------------------------
const NEW_SOLVERS = [
    "current: CrystalShift LM"             => current_lm,
    "LeastSquaresOptim Dogleg"             => DOGLEG,
    "NonlinearSolve LM, no geodesic"       => nls(NLS.LevenbergMarquardt(; autodiff = AutoForwardDiff(), disable_geodesic = Val(true))),
    "NonlinearSolve TrustRegion Hei"       => nls(NLS.TrustRegion(; autodiff = AutoForwardDiff(), radius_update_scheme = NLS.RadiusUpdateSchemes.Hei)),
    "NonlinearSolve TrustRegion Yuan"      => nls(NLS.TrustRegion(; autodiff = AutoForwardDiff(), radius_update_scheme = NLS.RadiusUpdateSchemes.Yuan)),
    "NonlinearSolve TrustRegion Fan"       => nls(NLS.TrustRegion(; autodiff = AutoForwardDiff(), radius_update_scheme = NLS.RadiusUpdateSchemes.Fan)),
    # TrustRegion with RadiusUpdateSchemes.Bastin is left out: it errors inside NonlinearSolve
    # ("No matching function wrapper was found!") on every problem.
    "NonlinearSolve TrustRegion NocedalWright" => nls(NLS.TrustRegion(; autodiff = AutoForwardDiff(), radius_update_scheme = NLS.RadiusUpdateSchemes.NocedalWright)),
    "Dogleg multi-start ×8"                => multistart(8),
    "Dogleg continuation (0.4, 0.2)"       => continuation(),
    "VarPro (Dogleg) + Dogleg polish"      => varpro(),
    "VarPro (Dogleg), no polish"           => varpro(polish = nothing),
]
