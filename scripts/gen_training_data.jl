using NPZ
using Glob
using ProgressBars

using CrystalShift
using CrystalShift: Gauss, Lorentz

width_interval = 0.02
strain_interval = 0.01
noise_level = 0.1
num_spectra = 4000

cifs = glob("data/calibration/References/*.cif")

q_range = (5, 50)
q = collect(LinRange(8., 40., 4501))
twoθ = q2twoθ.(q./10)

cs = CrystalPhase.(cifs, (q_range, ), 1:length(cifs), (.1,), (Lorentz(), ))

data = Array{Float64, 4}(undef, (length(cs), num_spectra, 4501, 1))

for i in tqdm(eachindex(cs))
    for j in 1:num_spectra
        p = get_free_params(cs[i])
        #scaling = (interval_size.*rand(size(params, 1)).-interval_size/2).+1
        #@. params = params*scaling
        #params = vcat(params, 0.5.+0.1randn(1), 0.1.+0.02(randn(1)))
        scaling = (strain_interval.*rand(length(p[1:end-2])) .- strain_interval/2) .+ 1
        p[end:end] = 0.1 .+ width_interval.*(randn(1))
        p[1:end-2] .*= scaling
        # TODO: optional add noise for bcnn training data
        t = zero(q)
        evaluate!(t, CrystalPhase(cs[i], p), q)
        t ./= maximum(t)
        noise = noise_level.*rand(length(q))
        t += noise
        t ./= maximum(t)
        data[i, j, :, :] = transpose(t)
    end
end

npzwrite("calib_train_n=20k_w=0.02_s=0.01_noise=0.1.npy", data)
