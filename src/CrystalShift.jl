module CrystalShift
using Base.Threads

using OptimizationAlgorithms: LevenbergMarquart, LevenbergMarquartSettings
using OptimizationAlgorithms: SaddleFreeNewton, DecreasingStep
using OptimizationAlgorithms: StoppingCriterion, fixedpoint!
using OptimizationAlgorithms: BFGS, LBFGS
using LinearAlgebra
using OptimizationAlgorithms
using SpecialFunctions
using ForwardDiff
using LazyInverses
using FastPow
using StatsBase
using Printf
using PrecompileTools
# using PyCall

# TODO: Update export list
export CrystalPhase, AbstractPhase, PeakModCP, Wildcard, Peak, PhaseModel, BackgroundModel
export Triclinic, Monoclinic, Orthorhombic, Tetragonal, Rhombohedral, Hexagonal, Cubic
export isCubic, isTetragonal, isHexagonal, isRhombohedral, isOrthohombic, isMonoclinic
export OptimizationMethods, OptimizationMode, OptimizationSettings # enums
export evaluate!, evaluate_residual!, optimize!, full_optimize!, fit_amorphous
export get_free_params, CifParser, CIFFile, Xtal, PowderDiffraction
export Gauss, Lorentz

# Python imports
# Note: Deprecated for ease of python wrapper installation
# Could be fixed by wrapping into docker image
# try
#     global github = ENV["GITHUB_WORKFLOW"]
# catch KeyError
#     global github = "false"
# end

# const CifParser = PyNULL()
# const CIFFile = PyNULL()
# const Xtal = PyNULL()
# const PowderDiffraction = PyNULL()

# function __init__()
#     if github == "false"
#         copy!(CifParser, pyimport("pymatgen.io.cif")."CifParser")
#         copy!(CIFFile, pyimport("xrayutilities.materials.cif")."CIFFile")
#         copy!(Xtal, pyimport("xrayutilities.materials.material")."Crystal")
#         copy!(PowderDiffraction, pyimport("xrayutilities.simpack")."PowderDiffraction")
#     end
# end

# Can be implemented to use AppleAccelerate
# function AppleAccelerate.exp(d::Vector{<:Dual{T}}) where T 
#     c = AppleAccelerate.exp(value.(d))
#     [Dual{T}(ci, ci*partials(di)) for (ci, di) in zip(c, d)]
# end


include("util.jl")
include("peakprofile.jl")
include("peak.jl")
include("crystal.jl")
include("crystalphase.jl")
include("wildcard.jl")
include("peakmodCP.jl")
include("background.jl")
include("fixedbackground.jl")
include("phasemodel.jl")
include("phaseresult.jl")
include("optimizationsettings.jl")
include("optimize.jl")

@setup_workload begin
    std_noise = 1.
    mean_θ = [1., .5, .2]
    std_θ = [.05, 2., 1.]
    test_path = "data/Ta-Sn-O/sticks.csv"
    x = collect(8:.1:60)
    y = zero(x)
    K = CovarianceFunctions.EQ()

    @compile_workload begin
        Gauss()
        Lorentz()
        FixedPseudoVoigt(0.5)
        PseudoVoigt(0.5)
        bg = BackgroundModel(x, K, 10)

        f = open(test_path, "r")
        CrystalPhase(f, 0.1, FixedPseudoVoigt(0.5) )

        f = open(test_path, "r")
        if Sys.iswindows()
            s = split(read(f, String), "#\r\n")
        else
            s = split(read(f, String), "#\n")
        end
        cs = CrystalPhase.(String.(s[1:end-1]), (0.1,), (FixedPseudoVoigt(0.5),));
        
        evaluate!(y, cs[1], [16.9, 4.86, 5.55, 1.5929519726520256, 0.687, 0.2, 0.5], x)
        optimize!(PhaseModel(cs[1:3]), x, y, std_noise, mean_θ, std_θ;
                  method = LM, maxiter = 32,
                  optimize_mode=Simple, em_loop_num=1,
                  regularization = true, verbose = false)
        optimize!(PhaseModel(cs[1:3]), x, y, std_noise, mean_θ, std_θ;
                  method = LM, maxiter = 32,
                  optimize_mode=EM, em_loop_num=1,
                  regularization = true, verbose = false)
        optimize!(cs[1:3], x, y, std_noise, mean_θ, std_θ;
                  method = LM, maxiter = 32,
                  optimize_mode=Simple, em_loop_num=1,
                  regularization = true, verbose = false)
        optimize!(PhaseModel(cs[1:3], nothing, bg), x, y, std_noise, mean_θ, std_θ;
                  method = LM, maxiter = 32,
                  optimize_mode=Simple, em_loop_num=1,
                  regularization = true, verbose = false)
        # Not doing this because it takes a long time, uncomment if you need to use Newton
        # optimize!(PhaseModel(cs[1:3], nothing, bg), x, y, std_noise, mean_θ, std_θ; 
        #           method = Newton, maxiter = 32,
        #           optimize_mode=Simple, em_loop_num=1,
        #           regularization = true, verbose = false)
    end
end


end
