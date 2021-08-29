struct PhaseResult
    phase::AbstractVector{<:CrystalPhase}
    acttivation::AbstractVector
    raw_spectrum::AbstractVector
    isCenter::Bool
end

struct StripeResult
    phase_results::AbstractVector{<:PhaseResult}
    x::Float64
    y::Float64
    cation_ratio::AbstractVector
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
