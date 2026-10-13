module CrystalShift
using Base.Threads

using OptimizationAlgorithms: LevenbergMarquart, LevenbergMarquartSettings
using OptimizationAlgorithms: SaddleFreeNewton, DecreasingStep
using OptimizationAlgorithms: StoppingCriterion, fixedpoint!
using OptimizationAlgorithms: BFGS, LBFGS
using LinearAlgebra
using OptimizationAlgorithms
import LeastSquaresOptim
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
export LM, Newton, bfgs, l_bfgs, dogleg
export Simple, EM, WithUncer

# Tree search exports
export Node, Tree, Lazytree, MPTree
export search!, search_k2n!
export TreeSearchSettings, MPTreeSearchSettings, FullOptimizeSettings
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

include("precompile.jl")

end
