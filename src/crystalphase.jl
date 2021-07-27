using Phasemapping: Lorentz

struct CrystalPhase{T, V<:AbstarctVector{T}, C, P}
    cl::C
    peaks::V

    id::Int
    name::String

    profile::P
end


function CrystalPhase(path::String, profile=Lorentz())
    f = readdlm(path)
    # File Spec:
    #    First line: Crystal info


    CrystalPhase(crystal, peaks, id, name, profile)
end
