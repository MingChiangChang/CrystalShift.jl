struct Peak
    h::Int16 # Stored separately to be explicit
    k::Int16 # Enforce h, k, l to be integer
    l::Int16

    q::Float64 # Reference for now, probably don't need in final version
    I::Float64 # Intensity
end

function Peak(s::String)
    info = split(s, ',')
    length(info) == 5 || throw("info must has length of 5")
    h, k, l = parse.(Int, info[1:3])
    q, I = parse.(Float64, info[4:5])
    Peak(h, k, l, q, I)
end

function get_peaks(lines)
    peaks = Vector{Peak}(undef, size(lines))
    for i in eachindex(lines)
        peaks[i] = Peak(String(lines[i]))
    end
    peaks, norm_constant = normalize_peaks!(peaks)
    return peaks, norm_constant
end

function normalize_peaks!(peaks::AbstractVector{Peak})
    intensities = [peaks[i].I for i in eachindex(peaks)]
    norm_constant = maximum(intensities)
    intensities ./= norm_constant
    for i in eachindex(peaks)
        peaks[i] = Peak(peaks[i].h, peaks[i].k, peaks[i].l, peaks[i].q, intensities[i])
    end
    return peaks, norm_constant
end

twoθ2q(twoθ::Float64, λ::Float64=1.5406) = 4*pi*sin(deg2rad(twoθ/2)) / λ