# Objects and functions for storing result and post analysis
struct PhaseResult{T, V<:AbstractVector{T}, K, M}
    phase::V
    phase_name::AbstractString
    activation::K
    raw_spectrum::M
    isCenter::Bool
end

function PhaseResult(cl::CrystalPhase, name::AbstractString,
                     H::AbstractVector, y::AbstractVector, isCenter::Bool)
    PhaseResult([cl], name, H, y, isCenter)
end

#(p::PhaseResult)(x::AbstractVector) = p.CP(x) # Quick Reconstruction

struct StripeResult{T, V<:AbstractVector{T}, K, M, L}
    phase_results::V
    cation_ratio::K
    x::M
    y::M
    tpeak::L
    dwell::L
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

# Collect parameter

# Generate plots
