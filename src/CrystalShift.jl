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
using Random
# using PyCall
using CovarianceFunctions: EQ
using Einsum
using LogExpFunctions: logsumexp

# TODO: Update export list
export CrystalPhase, AbstractPhase, PeakModCP, LinearPeakMod, Wildcard, Peak, PhaseModel, BackgroundModel
export FixedBackground, BasisBackground, PolynomialBackground, ChebyshevBackground, CosineBackground, InterpolationBackground
export Triclinic, Monoclinic, Orthorhombic, Tetragonal, Rhombohedral, Hexagonal, Cubic
export isCubic, isTetragonal, isHexagonal, isRhombohedral, isOrthohombic, isMonoclinic
export OptimizationMethod, OptimizationMode, OptimizationSettings # enums
export evaluate!, evaluate_residual!, optimize!, full_optimize!, fit_amorphous
export get_free_params, Gauss, Lorentz, FixedPseudoVoigt, PseudoVoigt
export LM, Newton, bfgs, l_bfgs
export Simple, EM, WithUncer

# Tree search exports
export Node, Tree, Lazytree, MPTree
export search!, search_k2n!
export TreeSearchSettings, MPTreeSearchSettings
export LeastSquares, KullbackLeibler, get_probabilities

# Python imports
# Note: Deprecated for ease of python wrapper installation
# Could be fixed by wrapping into docker image
# try
#     global github = ENV["GITHUB_WORKFLOW"]
# catch KeyError
#     global github = "false"
# end


# export CifParser, CIFFile, Xtal, PowderDiffraction # Only when using PyCall
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


# Forward model and optimization (formerly the whole of CrystalShift.jl)
include("shift/util.jl")
include("shift/peakprofile.jl")
include("shift/peak.jl")
include("shift/crystal.jl")
include("shift/crystalphase.jl")
include("shift/wildcard.jl")
include("shift/peakmodCP.jl")
include("shift/background.jl")
include("shift/fixedbackground.jl")
include("shift/basisbackground.jl")
include("shift/phasemodel.jl")
include("shift/phaseresult.jl")
include("shift/optimizationsettings.jl")
include("shift/lmsolver.jl")
include("shift/optimize.jl")

# Tree search and probabilistic labeling (merged from CrystalTree.jl)
include("tree/objective.jl")
include("tree/util.jl")
include("tree/node.jl")
include("tree/tree.jl")
include("tree/treesearchsettings.jl")
include("tree/lazytree.jl")
include("tree/mptree.jl")
include("tree/search.jl")
include("tree/probabilistic.jl")

# @setup_workload begin
#     std_noise = 1.0
#     mean_θ = [1.0, 0.5, 0.2]
#     std_θ = [0.05, 2.0, 1.0]
#     x = collect(8:0.1:60)
#     y = zero(x)
#     K = CovarianceFunctions.EQ()
#     lps = [[1., 1., 1., pi/2, pi/2, pi/2],
#             [1., 1., 2., pi/2, pi/2, pi/2],
#             [1., 2., 3., pi/2, pi/2, pi/2],
#             [1., 2., 3., pi/2, pi/2*0.98, pi/2],
#             [1., 1., 3., pi/2, pi/2, 2*pi/3],
#             [4., 5., 6., 2pi/3, 3pi/4, 5pi/6]]
#     hkls = [[[1, 1, 1]],
#             [[2, 3, 4]],
#             [[5, 6, 7]],
#             [[8, 9, 10]],
#             [[11, 12, 13]],
#             [[1, 3, 4]]]
#     peak_heights = [rand(1) for _ in 1:length(hkls) ]
#     _names = ["Cubic", "Tetragonal", "Orthorhombic", "Monoclinic", "Hexagonal", "Triclinic"]

#     @compile_workload begin
#         Gauss()
#         Lorentz()
#         FixedPseudoVoigt(0.5)
#         PseudoVoigt(0.5)
#         bg = BackgroundModel(x, K, 10)

#         cs = CrystalPhase.(lps, hkls, peak_heights, _names, collect(1:length(lps)), (0.1, ), (FixedApproxPseudoVoigt(0.5),))

#         evaluate!(y, cs[1], [16.9, 4.86, 5.55, 1.5929519726520256, 0.687, 0.2, 0.5], x)
#         optimize!(PhaseModel(cs[1:3]), x, y, std_noise, mean_θ, std_θ;
#             method=LM, maxiter=32,
#             optimize_mode=Simple, em_loop_num=1,
#             regularization=true, verbose=false)
#         optimize!(PhaseModel(cs[1:3]), x, y, std_noise, mean_θ, std_θ;
#             method=LM, maxiter=32,
#             optimize_mode=EM, em_loop_num=1,
#             regularization=true, verbose=false)
#         optimize!(cs[1:3], x, y, std_noise, mean_θ, std_θ;
#             method=LM, maxiter=32,
#             optimize_mode=Simple, em_loop_num=1,
#             regularization=true, verbose=false)
#         optimize!(PhaseModel(cs[1:3], nothing, bg), x, y, std_noise, mean_θ, std_θ;
#             method=LM, maxiter=32,
#             optimize_mode=Simple, em_loop_num=1,
#             regularization=true, verbose=false)
#         # Not doing this because it takes a long time, uncomment if you need to use Newton
#         # optimize!(PhaseModel(cs[1:3], nothing, bg), x, y, std_noise, mean_θ, std_θ; 
#         #           method = Newton, maxiter = 32,
#         #           optimize_mode=Simple, em_loop_num=1,
#         #           regularization = true, verbose = false)
#     end
# end


end
