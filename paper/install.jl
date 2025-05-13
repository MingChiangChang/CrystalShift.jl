using Pkg

packages = ["BackgroundSubtraction", "CSV", "Combinatorics",
            "CovarianceFunctions", "DataFrames", "ForwardDiff",
            "ProgressBars", "Plots", "JSON", "LaTeXStrings", "LazyInverses",
            "Measurements", "NPZ", "Plots", "ProgressBars", "StatsBase", "YAML"]

for package in packages
    Pkg.add(package)
end