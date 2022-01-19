module CrystalShift
using Base.Threads

using OptimizationAlgorithms: LevenbergMarquart, LevenbergMarquartSettings
using OptimizationAlgorithms: update_jacobian!
using LinearAlgebra
using OptimizationAlgorithms

export CrystalPhase
export Triclinic, Monoclinic, Orthorhombic, Tetragonal, Rhombohedral, Hexagonal, Cubic
export isCubic, isTetragonal, isHexagonal, isRhombohedral, isOrthohombic, isMonoclinic
export OptimizationMethods, Peak

include("util.jl")
include("peakprofile.jl")
include("peak.jl")
include("crystal.jl")
include("crystalphase.jl")
include("phaseresult.jl")
include("optimize.jl")

end
