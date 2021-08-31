struct PhaseResult
    phase::AbstractVector{<:CrystalPhase}
    activation::AbstractVector
    raw_spectrum::AbstractVector
    isCenter::Bool
end

function PhaseResult(CP::CrystalPhase, H::AbstractVector,
                     y::AbstractVector, isCenter::Bool)
    PhaseResult([CP], H, y, isCenter)
end

(p::PhaseResult)(x::AbstractVector) = p.CP(x) # Quick Reconstruction

struct StripeResult
    phase_results::AbstractVector{<:PhaseResult}
    cation_ratio::AbstractVector
    x::Float64
    y::Float64
    tpeak::Float64
    dwell::Float64
end

function get_center(SR::StripeResult)
    for pr in SR.phase_results
        if pr.isCenter
            return pr
        end
    end
end

# TODO some helper function to quickly access data for plotting or processing
# TODO Constructor of StripeResult using an array of phaseresult
