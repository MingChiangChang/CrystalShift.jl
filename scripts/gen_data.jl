using Plots
using NPZ
using CubicSplines

using CrystalShift
using CrystalShift: Gauss, q2twoθ

using Base.Threads

open("data/csw.txt") do f
    global txt = split(read(f, String), '\n')
end

function parsing(str::AbstractString)
    num, vec = split(str, '=')
    num = parse(Int64, num[2:end])
    vec = map(x->parse(Float64, x), split(vec, ","))
    return num, vec
end
composition = zeros(Float64, (231, 6))
shift = zero(composition)
width = zero(composition)

for i in eachindex(txt)
    if startswith(txt[i], "C")
        idx, vec = parsing(txt[i])
        composition[idx, :] = vec
    elseif startswith(txt[i], "S")
        idx, vec = parsing(txt[i])
        shift[idx, :] = vec
    elseif startswith(txt[i], "W")
        idx, vec = parsing(txt[i])
        width[idx, :] = vec
    end
end

test_path = "data/AlLiFeO/sticks_2.csv" # when ]test is executed pwd() = /test
# test_path = "/Users/ming/Downloads/AlLiFeO_new/sticks.csv"
f = open(test_path, "r")

if Sys.iswindows()
    s = split(read(f, String), "#\r\n")
else
    s = split(read(f, String), "#\n")
end

cs = CrystalPhase.(String.(s[1:end-1]), (0.1,), (Gauss(),))
# cs = cs[end-5:end]
x = collect(15:.1:79.9)
twoθ = q2twoθ.(x./10)
reg_itval_twoθ = collect(LinRange(minimum(twoθ), maximum(twoθ), 650))
ys = zeros(Float64, (231, 650))

# for i in eachindex(cs)
#     plt = plot(x, evaluate!(zero(x), cs[i], x), title=cs[i].name)
#     display(plt)
# end


p = Gauss()
for i in 1:231
    y = zero(x)
    for j in 1:6
        peak_locs = cs[j].cl.(cs[j].peaks) * 10 * shift[i, j]
        for idx in eachindex(peak_locs)
            y += cs[j].peaks[idx].I .* composition[i, j] .* p.((x .- peak_locs[idx])./(width[i, j]/4))
        end
    end
    ys[i, :] = y
    plt = plot(x, y, title="$(i)")
    display(plt)
end

npzwrite("alfeli_narrow.npy", ys)