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
export get_free_params, CifParser, CIFFile, Xtal, PowderDiffraction

# Python imports
# IDEA: Should be a github CI flag that can decern this
try
    global github = ENV["GITHUB_WORKFLOW"]
catch KeyError
    global github = "false"
end
# if github == "false"

# CifParser = pyimport("pymatgen.io.cif")["CifParser"]
# CIFFile = pyimport("xrayutilities.materials.cif")["CIFFile"]
# Xtal = pyimport("xrayutilities.materials.material")["Crystal"]
# PowderDiffraction = pyimport("xrayutilities.simpack")["PowderDiffraction"]

const CifParser = PyNULL()
const CIFFile = PyNULL()
const Xtal = PyNULL()
const PowderDiffraction = PyNULL()

function __init__()
    if github == "false"
        copy!(CifParser, pyimport("pymatgen.io.cif")["CifParser"])
        copy!(CIFFile, pyimport("xrayutilities.materials.cif")["CIFFile"])
        copy!(Xtal, pyimport("xrayutilities.materials.material")["Crystal"])
        copy!(PowderDiffraction, pyimport("xrayutilities.simpack")["PowderDiffraction"])
    end
end

# CifParser = pyimport("pymatgen.io.cif")["CifParser"]
# CIFFile = pyimport("xrayutilities.materials.cif")["CIFFile"]
# Xtal = pyimport("xrayutilities.materials.material")["Crystal"]
# PowderDiffraction = pyimport("xrayutilities.simpack")["PowderDiffraction"]
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
