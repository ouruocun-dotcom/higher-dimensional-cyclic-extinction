# Final pair statistics in higher dimensional cyclic systems

This repository contains the simulator, processed stochastic outputs, and plotting code accompanying a short numerical note on extinction cascades in cyclic populations with a switching carrying capacity.

The repository is intended to make the displayed numerical observations directly reproducible. It does not assign a broader interpretation to the observations.

## Contents

- `docs/Ou_extinction_cascade_observation.tex` is the LaTeX source for the two-page numerical note. Compile it from the `docs/` directory, or upload the repository contents to Overleaf and select this file as the main document.
- `figures/fragmentation_core_result.pdf` is the vector version of the main figure.
- `figures/fragmentation_core_result.png` is the screen-resolution version.
- `code/fast_switching_sim.mjs` runs the event-driven stochastic simulation.
- `code/plot_fragmentation_core.py` generates the main figure and summary statistics.
- `data/` contains the processed simulation outputs used by the plotting script.

## Model and recorded quantities

For species `i`, the simulation uses birth rate

```text
n_i [1 + s (A x)_i]
```

and density-dependent death rate

```text
n_i M / K(t).
```

`K(t)` switches between `K(1-delta)` and `K(1+delta)`. The interaction matrix is a weighted antisymmetric nearest-neighbour cycle. Each trajectory begins at the positive interior equilibrium and is simulated until a single species remains.

The output records:

- the complete extinction order;
- the time of every extinction;
- residence time at every surviving-species count;
- the identity of the final pair;
- the interaction magnitude within the final pair.

## Baseline design

- odd dimensions `q = 5, 7, 9`;
- selection strength `s = 0.5`;
- environmental amplitude `delta = 0.7`;
- switching rate `nu = 20`;
- effective population per initial species `K_eff/q = 20`;
- independently seeded exact event-driven trajectories.

Every JSON file stores its full configuration, derived parameters, equilibrium composition, trajectory counts, timing summaries, and random seed.

## Added interaction calculation

The `leaky_ring` topology retains the original nearest-neighbour cycle and adds antisymmetric interactions between previously unconnected pairs. The command-line option `--leakage` controls their common multiplier.

The displayed values are:

- `q = 7`: `0, 0.02, 0.05, 0.10, 0.25, 0.50`;
- `q = 9`: `0, 0.10, 0.50`.

## Reproduce the main figure

Python requirements:

```text
numpy
scipy
matplotlib
```

From the repository root:

```bash
python code/plot_fragmentation_core.py
```

The plotting script expects the JSON files in `data/` and writes the figure to `figures/` and the numerical summary to `data/`.

## Run a representative simulation

Node.js 18 or later is recommended.

```bash
node code/fast_switching_sim.mjs \
  --q 7 \
  --trials 2500 \
  --k-mean 275 \
  --delta-k 0.7 \
  --s 0.5 \
  --nu 20 \
  --heterogeneity 0.55 \
  --phase 0 \
  --full-fixation true \
  --topology ring \
  --design switching \
  --seed 20261023 \
  --output data/example_q7.json
```

To add the previously missing interactions, replace `--topology ring` with `--topology leaky_ring` and add, for example, `--leakage 0.25`.

## Reproduce the numerical checks

Environmental perturbations at `delta = 0.5` and `nu = 40` are summarized in `data/robustness_summary.json`.

A separate NumPy implementation can be compared with the optimized Node implementation using:

```bash
python code/validate_independent_implementations.py
```

The stored comparison is `data/independent_implementation_validation.json`. The script recomputes the homogeneity tests from stored outputs produced by independently written NumPy and Node implementations.

## Compile the note

From the `docs/` directory:

```bash
pdflatex Ou_extinction_cascade_observation.tex
pdflatex Ou_extinction_cascade_observation.tex
```

The second pass resolves the figure reference. On Overleaf, upload the complete repository contents and set `docs/Ou_extinction_cascade_observation.tex` as the main document. If Overleaf flattens the uploaded folder structure, move `fragmentation_core_result.pdf` beside the TeX file; the source accepts either location.

## Author

Ruocun Ou  
Nanyang Normal University  
ouruocun@gmail.com
