# Optimizer benchmark: findings

## Update 2026-10-06: Dogleg added as `method = dogleg`

Rerun on current `main` (files in `src/shift/`, tree search merged; the current LM now uses
`SplitJacobianLM`, with the exact background Jacobian, for linear backgrounds). Apple M3 Pro,
Julia 1.10.9.

**Decision:** `optimize!(...; method = dogleg)` runs LeastSquaresOptim's Dogleg on the same
problem as `lm_optimize!` (`dogleg_optimize!` in `src/shift/optimize.jl`). It is opt-in;
the default stays `LM`.

**Tree search needed a lattice bound.** The scenarios below only fit the *true* phases. In a tree
search most fitted models contain wrong phases, and there the unbounded Dogleg was worse than
the current LM. Its trust region takes large steps, so wrong phases strained 5–9% (once 40× on one
axis) to fit part of a pattern and pushed the true phases out of the top k. The current LM caps
every step at 0.1 in log space (`max_step`), which had been acting as protective regularization.
`dogleg_optimize!` therefore bounds each lattice parameter to ±`DOGLEG_MAX_STRAIN` (5%) of its
starting value. Phase assignment over 28 patterns (`tree.jl`; 7 phase pairs × perturbed/reference
peak heights × with/without noise; search depth 2, k 3, then `get_probabilities`):

| Method | correct (plain search) | search time | correct (with `FullOptimizeSettings` refinement) | search time |
|---|---|---|---|---|
| Current LM | 24/28 | 18.4 s | 28/28 | 44.2 s |
| **Dogleg, ±5% bound** | **24/28** | **9.3 s** | **28/28** | **26.4 s** |
| Dogleg, ±10% bound | 18/28 | – | 24/28 | – |
| Dogleg, unbounded | 15/28 | – | crashed (non-physical lattice in refinement) | – |

(4 threads. The ±10% and unbounded rows come from an earlier run of the same check with Julia 1.11;
their times aren't comparable, so they're left out.)

**Single fits (scenarios below, rerun with the bounded Dogleg).** The bound never binds there
(strains ≤ 1.25%), and Dogleg is still faster with equal or better lattice recovery:

| Scenario (trials) | Current LM: pass / lat ok / total time | Dogleg: pass / lat ok / total time |
|---|---|---|
| 1 phase, ±1.25% strain (60) | 100% / 100% / 368 ms | 100% / 100% / 127 ms |
| 2 phases, ±0.5% strain (40) | 92% / 88% / 627 ms (lat err p90 0.365%) | **100% / 95%** / 268 ms (p90 0.026%) |
| 1 phase + smooth background (15) | 100% / 100% / 307 ms | 100% / 100% / 81 ms |
| 1 phase + polynomial bg + noise (15) | 100% / 93% / 1.27 s | 100% / **100%** / 66 ms |
| wildcard + background, measured (1) | pass / 19.6 ms | pass / 3.6 ms |

`SplitJacobianLM` made the current LM faster on the smooth-background scenario than in the
original study (474 → about 270–310 ms), but Dogleg remains several times faster there.

Logs: `results/2026-10-06-scenarios-dogleg.log` (bounded Dogleg),
`results/2026-10-06-scenarios-packages.log` (unbounded, plus LeastSquaresOptim's LM),
`results/2026-10-06-tree-search.log`.

Notes:
- Dogleg runs to LeastSquaresOptim's default tolerances (1e-8), while the current LM stops once
  the residual norm is below 1e-2. On fits that start at the answer, the LM can therefore
  be faster.
- `WithUncer` uses the same least-squares Hessian for both methods; their uncertainties agree.
- New dependency: LeastSquaresOptim (pulls in Optim, LineSearches, NLSolversBase,
  DifferentiationInterface, FiniteDiff, ADTypes).

**Not done:** making Dogleg the default. It changes results for every user, so it needs more
evidence first, especially real (not synthetic) patterns and the strain priors used in practice.

---

# Original study (2026-10-03)

Based on commit `37155e7` (before the `src/shift/` reorganization). Machine: Apple M3 Pro,
Julia 1.10.9, 1 thread. File paths below are updated to the current layout; the numbers are
from that commit.

## Question

Can a well-maintained optimization package replace the `OptimizationAlgorithms.jl`
Levenberg–Marquardt used in `lm_optimize!` (`src/shift/optimize.jl`)? The replacement should be
faster without hurting fit quality, and above all without hurting **lattice-parameter
recovery**.

Tree search (`Lazytree`, `search!`) was out of scope; see the update above.

## Summary

- **Recommendation: replace the current LM with `LeastSquaresOptim.Dogleg()`.** Across every
  scenario it was 2.5–25× faster in total time. It passed every trial and recovered the
  lattice at least as well as the current LM: in two-phase fits, 95% of trials were within
  0.1% of the true lattice, vs 88% for the current LM.
- The current LM's two-phase failures are real lattice errors. The worst 10% of its trials
  are off by about 0.4%. Dogleg's worst 10% are within 0.026%.
- **VarPro + Dogleg polish** and **Dogleg multi-start ×8** reach a lower final cost on two-phase
  fits (median 0.176 vs 0.187). That only improved lattice recovery by one trial in 40
  (98% vs 95%). This isn't enough evidence yet to justify the extra complexity (VarPro) or
  the 8× cost (multi-start). A larger two-phase sample is the next step.
- **NonlinearSolve.jl** was slower than both the current LM and LeastSquaresOptim in almost
  every case. One trial of its LM without geodesic acceleration ran for about 5 minutes.
  The TrustRegion `Bastin` update rule errors on every problem.
- **Optim.jl** (BFGS/LBFGS) has no consistent advantage over the current BFGS. Both
  LBFGS implementations fail on the wildcard test.
- **Continuation** (fitting smoothed data first) and **VarPro without the polish step** didn't help.

## Method

Every solver gets exactly the problem that `simple_optimize!` builds:

- the same `PhaseModel`
- the same starting point, including `initialize_activation!`
- log-space parameters for phases and wildcards
- the same objective closure (`get_lm_objective_func` for least-squares solvers,
  `get_newton_objective_func` for quasi-Newton solvers)
- the same `maxiter`

Each package keeps its own default stopping rule, so time is always reported next to final
cost and accuracy.

Timing is a second pass over all trials, after an untimed first pass. Each crystal system is
its own Julia type, and two-phase models mix types, so a warm-up on a single problem leaves
compilation inside the timings. Before this fix the totals were inflated up to 10×.

### Scenarios (`benchmark/scenarios.jl`)

Each scenario copies the data generation, priors, `maxiter` and pass criterion of one test
file, with a fixed seed and more trials:

| Scenario | Source | Trials | Params | Pass criterion |
|---|---|---|---|---|
| 1 phase, ±1.25% lattice strain | `test/shift/optimize.jl` | 60 (each of the 15 Ta-Sn-O phases ×4) | 4–6 | ‖fit − y‖ < 0.1 |
| 2 random phases, ±0.5% strain | `test/shift/optimize.jl` | 40 | 8 | ‖fit − y‖ < 0.1 |
| 1 phase + smooth background (`BackgroundModel`) | `test/shift/background.jl` | 15 | 19 | ‖fit − y‖ < 0.1 |
| 1 phase + polynomial `FixedBackground` + noise | `test/shift/fixedbackground.jl` | 15 | 8 | MSE < 0.01 |
| Wildcard + background, measured data | `test/shift/wildcard.jl` | 1 | 13 | ‖fit − y‖ < 0.3 |

Metrics:

- **pass**: fraction of trials meeting the test's own pass criterion.
- **lat ok**: fraction of trials where every fitted free lattice parameter is within 0.1% of
  the true value used to generate the data. When the same phase is drawn twice, each fit is
  matched to its closest true lattice.
- **lat err med / p90**: median and 90th percentile of the largest relative lattice error
  in each trial.
- **≤ cur**: fraction of trials whose final cost is no worse than the current LM's (+1%).

## Results

### Lattice recovery and speed: two phases (40 trials)

This is the only scenario where the methods differ meaningfully.

| Method | pass | lat ok | lat err med / p90 | median cost | total time |
|---|---|---|---|---|---|
| Current LM (`OptimizationAlgorithms`) | 92% | 88% | 0.001% / **0.365%** | 0.1975 | 539 ms |
| **LSO Dogleg** | 100% | 95% | 0.001% / 0.026% | 0.1867 | **213 ms** |
| VarPro + Dogleg polish | 100% | **98%** | 0.001% / 0.025% | **0.1763** | 566 ms |
| Dogleg multi-start ×8 | 100% | **98%** | 0.001% / 0.025% | **0.1763** | 2.09 s |
| NLS TrustRegion Fan / NocedalWright | 100% | 95% | 0.001% / 0.026% | 0.1867 | 1.1–1.3 s |
| NLS TrustRegion Hei / Yuan | 98% | 92% | 0.001% / 0.037% | 0.1867 | 1.4 s / 10.1 s |
| Dogleg continuation (0.4, 0.2) | 98% | 92% | 0.001% / 0.037% | 0.1867 | 1.20 s |
| NLS LM, no geodesic | 92% | 88% | 0.001% / 0.139% | 0.1916 | 344 s (one trial ran ~5 min) |
| VarPro, no polish | 88% | 88% | 0.001% / 0.412% | 0.2099 | 416 ms |

### Total time per scenario (all trials, compiled)

| Method | 1 phase (60) | 2 phases (40) | smooth bg (15) | poly bg (15) | wildcard (1) |
|---|---|---|---|---|---|
| Current LM | 282 ms | 539 ms | 474 ms | 1.20 s | 19 ms |
| **LSO Dogleg** | **88 ms** | **213 ms** | **59 ms** | **47 ms** | **2.4 ms** |
| LSO LevenbergMarquardt | 86 ms | 241 ms | 78 ms | 78 ms ¹ | 3.3 ms |
| VarPro + Dogleg polish | 298 ms | 566 ms | 109 ms | 92 ms | 12.6 ms |
| Dogleg multi-start ×8 | 851 ms | 2.09 s | 593 ms | 459 ms | 20 ms |
| Dogleg continuation | 503 ms | 1.20 s | 216 ms | 173 ms | 4.6 ms |
| NLS TrustRegion (default) | 504 ms | 1.23 s | 268 ms | 253 ms | 14 ms |
| NLS TrustRegion Hei | 420 ms | 1.39 s | 512 ms | 223 ms | 10.5 ms |
| NLS LevenbergMarquardt (geodesic, default) | 2.32 s | 4.96 s | 1.87 s | 936 ms | 104 ms |
| NLS GaussNewton | 1.15 s | 5.82 s | 872 ms | 512 ms | 38 ms |

¹ LSO LevenbergMarquardt still passed every trial, but in 53% of them it stopped at a worse
minimum than the current LM (cost 40115 vs 38499). Dogleg didn't have this problem.

### Lattice recovery: single-phase scenarios

Every method recovers the lattice almost exactly in the single-phase scenarios, with the
worst 10% of trials within 0.002% (no background) or 0.042% (polynomial background). The
only misses:

- Polynomial background + noise: the **current LM** missed in 1 of 15 trials (p90 0.061%).
  Every other method got all 15 (p90 0.042%).
- Smooth background: **continuation** missed in 1 of 15.

### Quasi-Newton methods (wildcard test, measured data)

| Method | pass | final cost | time |
|---|---|---|---|
| Current BFGS | yes | 351.64 | 17–33 ms |
| Current LBFGS | **no** | 600.88 | 48–53 ms |
| Optim BFGS | yes | 351.64 | 10–26 ms |
| Optim LBFGS | **no** | 712.15 | 31–45 ms |

On the easy fixed problems in `benchmark/optimizers.jl`, Optim BFGS was faster than the
current BFGS for 1 phase (2.8 vs 4.6 ms) but much slower for 3 phases (1.46 s vs 12.6 ms).
It hit the iteration cap there. Both LBFGS implementations fail on 3 phases + background.

## Notes on the methods tried

- **Dogleg** (trust region): at each iteration it picks between the Gauss–Newton step and a
  steepest-descent step, limited by a trust radius that adapts to how well the local model
  predicted the last step. A rejected step only shortens the path and needs no new linear
  solve. The current LM can retry up to 32 times per iteration (`lm_backtrack!`), each with
  a new solve and fixed ×10 / ÷7 changes to λ, plus a 0.1 cap on each step in log space.
  Dogleg needed 6–8 iterations where the current LM needed 18–184.
- **VarPro** (prototype in `benchmark/methods.jl`): activations and background coefficients
  enter the model linearly, so they are solved by an inner linear least-squares solve. That
  solve covers the data residual plus the background ridge prior, and clamps activations
  positive. The outer Dogleg sees only the nonlinear parameters. The inner solve leaves out
  the activation prior, which is why the final full-problem Dogleg "polish" step is needed;
  without it, 12% of two-phase trials fail.
- **Multi-start**: 8 Dogleg runs from starting lattices jittered by up to ±1.25%; the best
  final cost is kept.
- **Continuation**: Dogleg on data smoothed with Gaussians of width 0.4 and then 0.2, then
  on the raw data.
- **NonlinearSolve**: `LevenbergMarquardt` has geodesic acceleration *on* by default. The
  TrustRegion rules tried were Hei, Yuan, Fan, Bastin and NocedalWright. Bastin errors inside
  NonlinearSolve on every problem ("No matching function wrapper was found!").

## Caveats

- Apart from the wildcard case, the scenarios use synthetic data from the test files, with
  small strains (±0.5–1.25%). Larger strains and real multi-phase patterns aren't covered yet.
- 40 two-phase trials are too few to separate Dogleg (38/40) from VarPro or multi-start (39/40).
- The wildcard scenario is a single run, so its timing is noisy: 19–47 ms for the current LM
  across runs.
- Single thread only. The tree search runs fits in parallel, which wasn't measured.
- Lower final cost didn't reliably mean better lattice recovery. Most of VarPro's cost gain
  came from activations and peak widths.

## Next steps

1. ~~Add `Dogleg` as a `method` option~~ Done (2026-10-06, with a lattice bound; see the
   update above). Switching the default is still open.
2. Rerun the two-phase scenario with more trials (`scenarios.jl 5 --set new`, 200 trials),
   limited to Dogleg, VarPro + polish and multi-start, to see if their edge holds up.
3. Add harder scenarios: larger strains, 3+ phases, real patterns from `data/`.
4. Independent of the solver: a hand-written Jacobian. One 3-phase ForwardDiff Jacobian
   takes 185 µs vs 59 µs for the residual alone (`benchmark/results/baseline.json`).

## Reproducing

```bash
# from the repo root; the first run installs the benchmark environment
julia --project=benchmark benchmark/run.jl [name] [--compare baseline] [--filter substring]  # general suite → results/<name>.json
julia --project=benchmark benchmark/optimizers.jl                    # fixed easy problems, BenchmarkTools timing
julia --project=benchmark benchmark/scenarios.jl [mult]              # test scenarios, package solvers
julia --project=benchmark benchmark/scenarios.jl [mult] --set new    # test scenarios, new strategies + lattice error (~15 min)
julia --project=benchmark benchmark/scenarios.jl --filter dogleg     # only the current LM and the package's dogleg
julia --project=benchmark -t 4 benchmark/tree.jl                     # phase assignment in tree search, LM vs dogleg
```

| File | Contents |
|---|---|
| `benchmarks.jl` | `SUITE`: forward model, objective/Jacobian, full `optimize!` |
| `run.jl` | Runs `SUITE`, saves JSON, compares against a saved baseline |
| `solvers.jl` | Shared `Problem` type and wrappers for the current solvers, LeastSquaresOptim, NonlinearSolve, Optim |
| `methods.jl` | Multi-start, continuation, VarPro prototype, NonlinearSolve variants |
| `optimizers.jl` | Package comparison on three fixed problems |
| `scenarios.jl` | Randomized scenarios from `test/shift/`, pass rate and lattice error |
| `tree.jl` | Phase assignment through tree search + `get_probabilities`, LM vs dogleg |
| `results/` | Raw logs of the runs above (`run.jl` writes JSON results here; they are not committed) |
