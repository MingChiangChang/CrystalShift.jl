using Plots

correct = [7728, 438, 311, 273, 240, 241, 295, 418, 832, 39208]
totl = [796356, 1526, 902, 738, 604, 575, 698, 828, 1411, 46380]

plt = plot([0., 1.], [0., 1.],
            linestyle=:dash, color=:black,
            legend=false, figsize=(10,10), dpi=300,
            xlims=(0, 1), ylims=(0, 1), xtickfontsize=10, ytickfontsize=10,
            xlabelfontsize=12, ylabelfontsize=12, xlabel="Predicted Probabilities", ylabel="Frequency of Correct Matches", markersize=5)

plot!(collect(0.05:.1: 0.95), calibration)
scatter!(collect(0.05:.1: 0.95), calibration)

plt = plot(legend=false, figsize=(10,10), dpi=300,
           xlims=(0.0, 1.0), ylims=(1.0, 10^6), xtickfontsize=10, ytickfontsize=10,
           xlabelfontsize=12, ylabelfontsize=12, xlabel="Predicted Probabilities", ylabel="Total Counts", markersize=5, yscale=:log10)

x = [0.05 + i*0.1 for i in 0:9]


for i in 1:10
    plot!([x[i], x[i]], [1.0, totl[i]*1.05], width=53, color=:black)
end

for i in 1:10
    plot!([x[i], x[i]], [1.0, totl[i]], width=51.5, color=:red)
end

display(plt)

