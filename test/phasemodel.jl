module Testphasemodel
using Test
using CrystalShift
using CrystalShift: evaluate!, get_param_nums, get_eight_params, get_free_lattice_params
using CrystalShift: collect_crystals, Lorentz, PseudoVoigt, get_moles, Gauss
using CrystalShift: get_fraction, get_free_params, evaluate_residual!, get_intrinsic_profile_type
using CrystalShift: get_phase_ids
using LinearAlgebra
using CovarianceFunctions: EQ

path = "../data/"
test_path = path * "Ta-Sn-O/sticks.csv" # when ]test is executed pwd() = /test
f = open(test_path, "r")


if Sys.iswindows()
    s = split(read(f, String), "#\r\n")
else
    s = split(read(f, String), "#\n")
end

cs = CrystalPhase.(String.(s[1:end-1]), (0.1,), (Lorentz(),))
gauss_cs = CrystalPhase(String(s[1]), 0.1, Gauss())
x = collect(8:.1:60)
y = zero(x)

empty_pm = PhaseModel()
@test isempty(empty_pm.CPs)
@test isnothing(empty_pm.background)

pm2 = PhaseModel(cs[1:2])
@test CrystalPhase(pm2) == cs[1:2]

pm3 = PhaseModel(cs[1])
bg = BackgroundModel(x, EQ(), 20)
pm = PhaseModel(cs[1], bg)
@test size(pm, 1) == 1
@test size(pm2, 1) == 2
@test pm3==pm

for i in pm2
    println(i.name)
end

@test isnothing(Nothing(nothing, 3))
@test evaluate!([1,2,3], nothing, [1.,.3,2.], x) == [1,2,3]
@test evaluate!([1,2,3], nothing, x) == [1,2,3]
@test evaluate_residual!(nothing, x, [1.,.3,2.], [1,2,3]) == [1,2,3]
@test evaluate_residual!(nothing, x, [1,2,3]) == [1,2,3]

@test pm(x) == evaluate!(zero(x), cs[1], x)
@test get_param_nums(pm3) == get_param_nums(cs[1])
@test get_param_nums(pm) == get_param_nums(cs[1]) + get_param_nums(bg)
@test get_phase_ids(pm) == [0]
@test get_phase_ids(pm2) == [0,1]
end