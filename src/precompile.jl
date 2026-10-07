# Precompile workload: runs a small version of the usual workflow so that its methods
# are compiled into the package image instead of on first use. ForwardDiff compiles
# separate code per chunk size, which by default is the number of parameters (up to
# 12), so the workload fits single phases of every lattice system and all their pairs
# to cover the chunk sizes of one- and two-phase models.
@setup_workload begin
    # one small phase per lattice system: a, b, c, α, β, γ (degrees)
    lattices = [[4.1, 4.1, 4.1, 90, 90, 90],     # cubic
                [4.2, 4.2, 6.1, 90, 90, 90],     # tetragonal
                [3.3, 3.3, 5.2, 90, 90, 120],    # hexagonal
                [5.4, 5.4, 5.4, 80, 80, 80],     # rhombohedral
                [4.3, 5.1, 6.2, 90, 90, 90],     # orthorhombic
                [4.4, 5.2, 6.3, 90, 101, 90],    # monoclinic
                [4.5, 5.3, 6.4, 82, 95, 103]]    # triclinic
    hkls = [(1, 0, 0), (1, 1, 0), (1, 1, 1), (2, 0, 0), (2, 1, 0), (2, 1, 1), (0, 1, 2), (1, 0, 2)]
    # CrystalPhase(io) splits phases on "#\r\n" on Windows, so use the platform's line ending
    nl = Sys.iswindows() ? "\r\n" : "\n"
    blocks = map(enumerate(lattices)) do (i, lp)
        header = "$(i-1),phase$(i-1),system," * join(lp, ",")
        peaks = ["$h,$k,$l,0.0,$(100 / j)" for (j, (h, k, l)) in enumerate(hkls)]
        join(vcat(header, peaks), nl)
    end
    csv = join(blocks .* "#" .* nl) # stick-pattern CSV format: each phase ends with "#"
    x = collect(range(8., 45., length = 256))
    std_noise, mean_θ, std_θ = .01, [1., 1., .2], [.5, .5, 1.]
    opt_stn = OptimizationSettings{Float64}(std_noise, mean_θ, std_θ, 4)

    @compile_workload begin
        path, io = mktemp()
        write(io, csv); close(io)
        cs = open(f -> CrystalPhase(f, 0.1, FixedPseudoVoigt(0.5)), path)
        rm(path)
        y = evaluate!(zero(x), cs[1], x) .+ evaluate!(zero(x), cs[2], x)
        y ./= maximum(y)

        # all one- and two-phase models: covers the ForwardDiff chunk sizes 3 to 12, for
        # the fit (Jacobian) and for the model probabilities (Hessian)
        nodes = Node[]
        for i in eachindex(cs), j in i:length(cs)
            phases = i == j ? cs[i:i] : cs[[i, j]]
            pm = optimize!(PhaseModel(phases), x, y, std_noise, mean_θ, std_θ; method = LM, maxiter = 2)
            push!(nodes, Node(pm, Node[], length(nodes) + 1, x, y))
        end
        get_probabilities(nodes, x, y, std_noise, mean_θ, std_θ)

        # EM mode and a kernel background, for every lattice system
        bg = BackgroundModel(x, EQ(), 5.)
        for i in eachindex(cs)
            j = mod1(i + 1, length(cs))
            optimize!(PhaseModel(cs[[i, j]]), x, y, std_noise, mean_θ, std_θ; method = LM, maxiter = 2,
                      optimize_mode = EM, em_loop_num = 2) # 2: also runs the noise re-estimation step
            optimize!(PhaseModel(cs[[i, j]], nothing, bg), x, y, std_noise, mean_θ, std_θ; method = LM, maxiter = 2)
        end
        full_optimize!(PhaseModel(cs[1:2]), x, y, std_noise, mean_θ, std_θ; method = LM, loop_num = 1,
                       peak_shift_iter = 2, peak_mod_iter = 2)

        # tree search, with and without peak-height refinement, and probabilities
        for full_opt_stn in (nothing, FullOptimizeSettings(loop_num = 1, peak_mod_iter = 2))
            ts_stn = TreeSearchSettings{Float64}(2, 2, false, false, 5., opt_stn; full_opt_stn)
            res = reduce(vcat, search!(Lazytree(cs, x), x, y, ts_stn)[2:end])
            get_probabilities(res, x, y, std_noise, mean_θ, std_θ)
            get_probabilities(res, x, y, mean_θ, std_θ)
        end
        search_k2n!(Lazytree(cs, x), x, y, TreeSearchSettings{Float64}(2, 2, false, false, 5., opt_stn))

        # other peak profiles (the string constructor defaults to PseudoVoigt), for every
        # lattice system
        for profile in (Lorentz(), Gauss(), PseudoVoigt(0.5))
            p = CrystalPhase.(String.(blocks), 0.1, (profile,))
            for i in eachindex(p)
                optimize!(PhaseModel(p[[i, mod1(i + 1, length(p))]]), x, y, std_noise, mean_θ, std_θ;
                          method = LM, maxiter = 2)
            end
        end
    end
end
