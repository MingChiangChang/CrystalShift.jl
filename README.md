# CrystalShift.jl
[![CI](https://github.com/MingChiangChang/CrystalShift.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/MingChiangChang/CrystalShift.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/MingChiangChang/CrystalShift.jl/branch/main/graph/badge.svg?token=3A8XI43H0C)](https://codecov.io/gh/MingChiangChang/CrystalShift.jl)

`CrystalShift.jl` is a package for probabilistic multi-phase labeling of X-ray diffraction(XRD) patterns. It is an efficient and flexible way to probabilistically label phases in XRD patterns based on pseudo-refinement optimzation, best-first tree search and Bayesian model comparison techniques.

## Installation
```julia
using Pkg
Pkg.add(url="https://github.com/MingChiangChang/CrystalShift.jl")
Pkg.add(url="https://github.com/MingChiangChang/CrystalTree.jl")
```
There are plans to merge these two module and add it to the julia general repo in the future.

## Usage
### Generate precalculate input files
```console
python cif_to_input_file.py -c /path/to/cif/folder -o /path/to/output/location -qmin 10.0 -qmax 80.0 -w 1.5406
```
This script will grab all of the cifs in /path/to/cif/folder, precalculate the XRD stick patterns, and store it as `/path/to/output/location.csv`. 

### Pseudo-refinement of lattices
Once you have the input files generated, it can be used to generate an array of `CrystalPhase` objects
```julia
using CrystalShift
cs = CrystalPhase("/path/to/output/location.csv", 0.1, CrystalShift.FixedPseudoVoigt(0.5))
```
The array of `CrystalPhase` object can then be passed to the pseudo-refinement optimization process
```julia
std_noise = .1
mean_θ = [1., .5, .2]
std_θ = [.05, 2., 1.]
# y is the target XRD, which has to be normalized to have maximum of 1
phasemodel = CrystalShift.optimize!(cs, q, y, std_noise, mean_θ, std_θ) 
```
This is the minimal setup for the optimization. User can also pass in the following named arguments:
* `y_uncer::AbstractVector`: uncertainty of y that will be take into account when optimizing
* `method::OptimizationMethod`: `LM`(Levenberg-Marquart), `Newton`, `bfgs`, `l_bfgs`
* `optimize_mode::OptimizationMode`: `Simple`, `EM`(expectation-maximization mode for jointly optimize `std_noise`), `WithUncer` (returns uncertainty of lattice parameters)
* `ovejctive::String`: `"LS"` (least square), `"KL"` (KL-divergence)
* `maxiter::Int`: number of maximum iterations
* `regularization::Bool`: whether to include regularization on lattice strain, peak width and activation of phase activations
* `em_loop_num::Int`: how many iterations of EM performed
* `λ::Float64`: For regularizing background

One can also pass in `PhaseModel` array, of which each contains a list of `CrystalPhase`, a `BaackgroundModel` and a `Wildcard` object.
`BackgroundModel` is a modified kernel regression model that construct a Gram matrix based on the kernel function you selected, perform SVD to obtain the most significant basis and use that to form a linear model to simulate the background.
`Wildcard` is an object that tries to optimize signals that are not explained by the above models with multiple manually construct but modifiable gaussian distributions.
All three objects are jointly optimized if `PhaseModel` is passed in to `optimize!` function.

### Tree search and probability estimates
```julia
using CrystalTree: LazyTree, search!, get_probabilities

max_depth = 3         # maximum number of allowed phases
k = 3                 # degree of expansion for each top node
amorphous = false     # whether trying to fit the data with an smooth signal
background = false    # whether
background_length = 8.

lt = LazyTree(cs, q)
results = search!(lt, q, y, max_depth, k, amorphous, background, background_length)

results = results[2:end]         # the result is a vector with (k+1) length, one for each level.
                                 # If you did not try to optimize for amorphous phase, you will have to exclude the root node.
results = reduce(vcat, results)  # flatten the result

# std_noise is optional, can be calculated instead of a preset number
probibilities = get_probabilities(results, q, y, std_noise, mean_θ, std_θ, renormalize=true, normalization_constant=2.5)
```
To obtain a proper normalization constant, one will have to go through the optimization procedure described in the paper.

## Cite this package
When using this package for your work, please cite this package using the following Bibtex citation:
```
@article{chang2025probabilistic,
  title={Probabilistic phase labeling and lattice refinement for autonomous materials research},
  author={Chang, Ming-Chiang and Ament, Sebastian and Amsler, Maximilian and Sutherland, Duncan R and Zhou, Lan and Gregoire, John M and Gomes, Carla P and van Dover, R Bruce and Thompson, Michael O},
  journal={npj Computational Materials},
  volume={11},
  number={1},
  pages={148},
  year={2025},
  publisher={Nature Publishing Group UK London}
}```

## Other links
* [CrystalTree.jl](https://github.com/MingChiangChang/crystaltree.jl) is a package based on `CrystalShift` that builds all of the tree search and probibilistic capabilities.
* [pyPhaseLabel](https://github.com/MingChiangChang/pyPhaseLabel) is a python wrapper of `CrystalShift`.
* [phiddle](https://github.com/MingChiangChang/phiddle) is a GUI that builds on top of `CrysatlShift` for rapid labeling of large XRD dataset and build-in visualization tools.
Currently is designed for the need for our group but it is easily modifiable for other usage.
