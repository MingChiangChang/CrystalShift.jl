module TestWildcard
using CrystalShift
using CrystalShift: Wildcard, Lorentz, evaluate, optimize!
using CovarianceFunctions: EQ

using NPZ 
using Test
using LinearAlgebra
# using Plots 

x = collect(10:.1:42)
w = Wildcard([20., 35.], [1., 0.2],  [2., 3.], "Amorphous", Lorentz(), [2., 2., 1., 1., .2, .5])
#w = Wildcard([20., ], [1., ], "Amorphous", [2.,], Lorentz(), [1., 1., .5,])
# plot(x, evaluate(w, x))

q = npzread("../data/test_q.npy")
y = npzread("../data/test_int.npy")
y ./= maximum(y)*2
# plt = plot(q, y)
bg = BackgroundModel(q, EQ(), 20, rank_tol=1e-3)
# plot!(q, evaluate!(zero(q), w, [22., 36., 1., 0.5, 2., 2.], q))

pm = PhaseModel(w, bg)

new_pm = optimize!(pm, q, y, 1e-2, [1.,1., 1.], [1., 1., 1.], method=LM, objective="LS",
                         maxiter=512, regularization=true, verbose=false)
t = zero(q)
# plot!(q, evaluate!(t, new_pm, q))
# plot!(q, evaluate!(zero(q), new_pm.wildcard, q))
# plot!(q, evaluate!(zero(q), new_pm.background, q))
# display(plt)
@test norm(y-evaluate!(t, new_pm, q)) < 0.3

end # module