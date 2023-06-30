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

end
