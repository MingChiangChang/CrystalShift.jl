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
using StatsBase
using PyCall

# TODO: Update export list
export CrystalPhase, AbstractPhase, PeakModCP, Wildcard, Peak, PhaseModel, BackgroundModel
export Triclinic, Monoclinic, Orthorhombic, Tetragonal, Rhombohedral, Hexagonal, Cubic
export isCubic, isTetragonal, isHexagonal, isRhombohedral, isOrthohombic, isMonoclinic
export OptimizationMethods, OptimizationMode, OptimizationSettings # enums
export evaluate!, evaluate_residual!, optimize!, full_optimize!, fit_amorphous
export get_free_params

# Python imports
# IDEA: Should be a github CI flag that can decern this
ci = true # Change this to false once PyCall is properly install in the right env
if !ci
    global CifParser = pyimport("pymatgen.io.cif")["CifParser"]
    global CIFFile = pyimport("xrayutilities.materials.cif")["CIFFile"]
    global Xtal = pyimport("xrayutilities.materials.material")["Crystal"]
    global PowderDiffraction = pyimport("xrayutilities.simpack")["PowderDiffraction"]
end

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
