using LinearAlgebra
using OptimizationAlgorithms: LevenbergMarquart
using Plots

coeff = randn(5)

function fe(coeff, x)
    a = zeros(size(x))
    for i in eachindex(coeff)
        @. a += coeff[i]*x^i
    end
    a
end

x = collect(-1:.001:1)

plot(x, fe(coeff, x))

function optimize!(θ::AbstractVector,
                   x::AbstractVector, y::AbstractVector; maxiter::Real=32)
    function residual!(r::AbstractVector, θ::AbstractVector)
        params = exp.(θ) # make a copy
		println(typeof(r), typeof(y), typeof(fe(θ, x)))
        @. r = y
		println(size(θ), size(x), size(r), size(fe(θ, x)))
		@. r -= fe(θ, x)
        return r
	end

    #@. θ = log(θ) # tramsform to log space for better conditioning
	r = zeros(eltype(θ), size(y))
    LM = LevenbergMarquart(residual!, θ, r)
    stn = LevenbergMarquartSettings(min_resnorm = 1e-2, min_res = 1e-3,
						min_decrease = 1e-8, max_iter = maxiter,
						decrease_factor = 7, increase_factor = 10, max_step = 1.0)
	λ = 1e-6
	OptimizationAlgorithms.optimize!(LM, θ, copy(r), stn, λ, Val(true))
	#@. θ = exp(θ) # transform back
	println(θ)
	return θ
end

θ = randn(5)

θ = optimize!(θ, x, fe(coeff, x))

plot!(x, fe(θ, x))
