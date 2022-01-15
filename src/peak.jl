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

twoθ2q(twoθ::Float64, λ::Float64=1.5406) = 4*pi*sin(deg2rad(twoθ/2)) / λ