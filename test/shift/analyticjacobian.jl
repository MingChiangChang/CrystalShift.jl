module TestAnalyticJacobian
using CrystalShift
using CrystalShift: CrystalPhase, PhaseModel, OptimizationSettings, LM, Peak, optimize!
using CrystalShift: Lorentz, Gauss, PseudoVoigt, FixedPseudoVoigt, ApproxGaussian
using CrystalShift: Cubic, Tetragonal, Hexagonal, Rhombohedral, Orthorhombic, Monoclinic, Triclinic
using CrystalShift: block_lm, get_lm_objective_func, get_free_params, get_free_lattice_params,
                    get_intrinsic_crystal_type, ADBlock, PhaseJacobianBlock, PeakPositions,
                    peak_positions!, analytic_jacobian!, block_jacobian!, USE_ANALYTIC_JACOBIAN
using CrystalShift.OptimizationAlgorithms: update_jacobian!
using ForwardDiff
using LinearAlgebra
using Test

relerr(a, b) = maximum(abs.(a .- b)) / maximum(abs.(b))

f = open("../data/tree_test_sticks.csv", "r")
s = Sys.iswindows() ? split(read(f, String), "#\r\n") : split(read(f, String), "#\n")
close(f)
s = filter(!isempty ∘ strip, s)

x = collect(LinRange(8, 45, 700))
y_uncer = @. 0.005 * (1 + sin(x)) # non-constant weights
std_noise, mean_θ, std_θ = .01, [1., 1., .2], [.5, .5, 1.]

# The LM blocks of a phase model, evaluated away from the reference parameters
function lm_blocks(pm, y, reg; y_uncer = zero(y))
    opt = OptimizationSettings{Float64}(std_noise, mean_θ, std_θ, 64, reg, LM)
    θ = get_free_params(pm)
    log_θ = log.(θ) .+ 0.003 .* sin.(1:length(θ))
    f = get_lm_objective_func(pm, x, y, y_uncer, opt)
    r = zeros(reg ? length(y) + length(log_θ) : length(y))
    block_lm(f, pm, x, y, y_uncer, log_θ, r, opt), log_θ
end

function forwarddiff_block(b::PhaseJacobianBlock, θ)
    J, out = similar(b.J), similar(b.out)
    ForwardDiff.jacobian!(J, b.ad.g!, out, θ, b.ad.cfg)
    J, out
end

@testset "Peak positions from 6 basis reflections" begin
    hkl = [(1, 0, 0), (1, 1, 0), (1, 1, 1), (2, 1, 0), (2, -1, 1), (3, 1, -2), (0, 0, 2), (1, -2, 3), (4, 2, 1)]
    peaks = [Peak(Int16(h), Int16(k), Int16(l), 0.0, 1.0) for (h, k, l) in hkl]
    cells = (Cubic{Float64}(4.1), Tetragonal{Float64}(4.1, 6.3), Hexagonal{Float64}(3.2, 5.1),
             Rhombohedral{Float64}(5.4, 1.2), Orthorhombic{Float64}(4.1, 5.2, 6.3),
             Monoclinic{Float64}(4.1, 5.2, 6.3, 1.75), Triclinic{Float64}(4.1, 5.2, 6.3, 1.5, 1.7, 1.9))
    for cl in cells
        P = PeakPositions(cl, peaks)
        ℓ = log.(get_free_lattice_params(cl))
        @test peak_positions!(P, ℓ)
        T = get_intrinsic_crystal_type(typeof(cl))
        q(ℓ) = (c = T{eltype(ℓ)}(exp.(ℓ)...); [10 * c(p) for p in peaks])
        @test relerr(P.q, q(ℓ)) < 1e-12
        @test relerr(P.Dq, ForwardDiff.jacobian(q, ℓ)) < 1e-10
    end
    # impossible cell: volume is Inf and every q is 0, so the analytic path must decline
    bad = Triclinic{Float64}(4., 5., 6., 2.5, 2.5, 2.5)
    @test !peak_positions!(PeakPositions(bad, peaks), log.(get_free_lattice_params(bad)))
end

@testset "Phase blocks match ForwardDiff" begin
    for profile in (Lorentz(), Gauss(), PseudoVoigt(0.5), FixedPseudoVoigt(0.5)), reg in (true, false)
        cs = CrystalPhase.(String.(s), (0.1,), (profile,))
        for cp in cs
            y = cp.(x); y ./= maximum(y)
            LMb, log_θ = lm_blocks(PhaseModel([cp]), y, reg; y_uncer)
            b = only(LMb.blocks)
            @test b isa PhaseJacobianBlock
            θ = log_θ[b.cols]
            @test analytic_jacobian!(b, θ)
            J_fd, out_fd = forwarddiff_block(b, θ)
            @test relerr(b.J, J_fd) < 1e-10
            @test relerr(b.out, out_fd) < 1e-10
        end
    end
end

@testset "Assembled Jacobian of a two-phase model" begin
    cs = CrystalPhase.(String.(s), (0.1,), (PseudoVoigt(0.5),))
    pm = PhaseModel(cs[[1, 6]])
    y = cs[1].(x) .+ cs[6].(x); y ./= maximum(y)
    LMb, log_θ = lm_blocks(pm, y, true)
    @test all(b -> b isa PhaseJacobianBlock, LMb.blocks)
    y_an, J_an = copy.(update_jacobian!(LMb, log_θ))
    try
        USE_ANALYTIC_JACOBIAN[] = false
        y_fd, J_fd = update_jacobian!(LMb, log_θ)
        @test relerr(J_an, J_fd) < 1e-10
        @test relerr(y_an, y_fd) < 1e-10
    finally
        USE_ANALYTIC_JACOBIAN[] = true
    end
end

@testset "Fallback to ForwardDiff" begin
    # unsupported profile (no finite support)
    cp = CrystalPhase(String(s[1]), 0.1, ApproxGaussian())
    y = cp.(x); y ./= maximum(y)
    @test only(lm_blocks(PhaseModel([cp]), y, true)[1].blocks) isa ADBlock
    # unsorted grid
    cp = CrystalPhase(String(s[1]), 0.1, Lorentz())
    opt = OptimizationSettings{Float64}(std_noise, mean_θ, std_θ, 64, true, LM)
    pm, xr = PhaseModel([cp]), reverse(x)
    yr = cp.(xr); yr ./= maximum(yr)
    log_θ = log.(get_free_params(pm))
    LMb = block_lm(get_lm_objective_func(pm, xr, yr, zero(yr), opt), pm, xr, yr, zero(yr), log_θ,
                   zeros(length(yr) + length(log_θ)), opt)
    @test only(LMb.blocks) isa ADBlock
    # switched off: a phase block gives exactly the ForwardDiff result
    y = cp.(x); y ./= maximum(y)
    LMb, log_θ = lm_blocks(PhaseModel([cp]), y, true)
    b = only(LMb.blocks)
    θ = log_θ[b.cols]
    J_fd, out_fd = forwarddiff_block(b, θ)
    try
        USE_ANALYTIC_JACOBIAN[] = false
        block_jacobian!(b, θ)
        @test b.J == J_fd
        @test b.out == out_fd
    finally
        USE_ANALYTIC_JACOBIAN[] = true
    end
end

@testset "Fits match ForwardDiff fits" begin
    for profile in (Lorentz(), PseudoVoigt(0.5))
        cs = CrystalPhase.(String.(s), (0.1,), (profile,))
        strained(cp, f) = (θ = copy(get_free_params(cp)); θ[1:cp.cl.free_param] .*= f; CrystalPhase(cp, θ))
        y = strained(cs[1], 1.005).(x) .+ strained(cs[6], 0.996).(x); y ./= maximum(y)
        opt = OptimizationSettings{Float64}(std_noise, mean_θ, std_θ, 128, true, LM)
        θ_an = get_free_params(optimize!(PhaseModel(cs[[1, 6]]), x, y, opt))
        θ_fd = try
            USE_ANALYTIC_JACOBIAN[] = false
            get_free_params(optimize!(PhaseModel(cs[[1, 6]]), x, y, opt))
        finally
            USE_ANALYTIC_JACOBIAN[] = true
        end
        @test θ_an ≈ θ_fd rtol = 1e-8
    end
end

println("End of analyticjacobian.jl test")
end # module
