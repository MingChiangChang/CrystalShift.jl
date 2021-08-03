"""
    optimize!

    This is the function that each of the node would call on.
    Try to fit a vector of CrystalPhase to the given spectrum
"""
function optimize!(phases::AbstractVector{:<CrystalPhase},
                   x::AbstractVector, y::AbstractVector,
                   std_noise::Reak = .01, mean_θ::AbstractVector = [1., 1.,.2],
                   std_θ::AbstractVector = [3., .01, 1.];
                   max_iter::Int = 32, regularization::Bool = true)

end

function optimize!(phase::CrystalPhase, x::AbstractVector, y::AbstractVector,
                   std_noise::Reak = .01, mean_θ::AbstractVector = [1., 1.,.2],
                   std_θ::AbstractVector = [3., .01, 1.],
                   max_iter::Int = 32, regularization::Bool = true)
    optimize!([phase], x, y, std_noise, mean_θ, std_θ, max)
end

function optimize!(θ::AbstracVector, x::AbstractVector, y::AbstractVector,
                   std_noise::Real, mean::AbstractVector,
                   std_θ::AbstractVector; maxiter::Int = 32,
                   regularization::Bool = true)
    # θ are the parameters
end


function optimize!(θ::AbstracVector, phases::AbstractVector{:<CrystalPhase},
                   x::AbstractVector, y::AbstractVector,
                   std_noise::Real, mean::AbstractVector,
                   std_θ::AbstractVector; maxiter::Int = 32,
                   regularization::Bool = true)
    function residual(θ::AbstractVector, )
end
