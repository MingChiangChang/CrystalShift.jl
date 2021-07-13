using Plots
using DelimitedFiles
using PhaseMapping
using PhaseMapping: readsticks, Lorentz, Phase, StickPattern
using PhaseMapping: pmp_path!, Library
using Test

# Backdoor
function Phase(c, μ, id::Int, a::Real, α::Real, σ::Real; profile = Lorentz(), width_init::Real = 1.)
	length(c) == length(μ) || throw(DimensionMismatch())
	c, μ = promote(c, μ)
	dc = zero(c)
	T = eltype(c)
	a, α, σ = a, α, σ
    Phase(c, μ, id, dc, a, α, σ, profile)
end

function Phase(S::StickPattern, a::Float64, α::Float64; profile = Lorentz(), width_init::Real = 1.)
    Phase(S.c, S.μ, S.id, a, α, width_init, profile = profile)
end

function get_ids(lib::Library)
    [p.id for p in lib.phases]
end

path = "C:/Users/r2121/Desktop/Code/crystal_shift/test/"
path = path*"sticks.txt"

sticks = readsticks(path, Float64)
println(typeof(sticks[1]))
phase_1 = Phase(sticks[1], 1.0, 1.0, profile=Lorentz(), width_init=.1)
phase_2 = Phase(sticks[2], 0.5, 1.0, profile=Lorentz(), width_init=.2)
#println(size(phase_arr))
phases = Phase.(sticks, profile=Lorentz(), width_init=.2)

x = LinRange(8, 45, 1024)
y = phase_1.(x)+phase_2.(x)

libraries, residuals = pmp_path!(phases, x, y, 2)
ids = get_ids(libraries[2])
@test 1 in ids
@test 2 in ids
@test 0.95<=libraries[2].phases[1].a<=1
@test 0.45<=libraries[2].phases[2].a<=0.55
