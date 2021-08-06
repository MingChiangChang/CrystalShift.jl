# test peak construction
using Plots

include("../src/Crystal.jl")
include("../src/peak.jl")

cl = Cubic{Float64}(1.0, 1.0, 1.0, pi/2, pi/2, pi/2)

# Peak list creation
# Spec:
p = Peak(1, 0, 0, 10)


# Reconstruction using cl and compare with Phasemapping?
cl(p)
