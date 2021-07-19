function optimize!(θ::AbstracVector, x::AbstractVector, y::AbstractVector,
                   std_noise::Real, mean::AbstractVector,
                   std_θ::AbstractVector; maxiter::Int = 32,
                   regularization::Bool = true)
    # θ are the parameters
end


function optimize!(θ::AbstracVector, phases::AbstractVector{:<Crystal},
                   x::AbstractVector, y::AbstractVector,
                   std_noise::Real, mean::AbstractVector,
                   std_θ::AbstractVector; maxiter::Int = 32,
                   regularization::Bool = true)
    function residual(θ::AbstractVector, )
end
