# Code for the Brinkman local-ROM paper

This MATLAB package reproduces the numerical experiments in the manuscript:

> X. Li, L. Ju, and Z. Wang, *Residual-Screened Local Reduced-Order Modeling
> for Brinkman-Penalized Shape Optimization of Incompressible Flows*.

The four main cases are the Stokes cannula, the isolated body at Reynolds
numbers 40 and 200, and Pironneau's three-dimensional Stokes body. The main
comparison is with a **preconditioned full-order optimization**, whose state
and adjoint solves use GMRES and a retained sparse factorization. The compact
local ROM uses eight state and adjoint snapshots and screened defect
corrections before a full-order fallback.

## Requirements and setup

- MATLAB R2024b; no extra toolboxes are used.
- Several gigabytes of memory for the three-dimensional sparse LU factors.
- A long, otherwise idle session for comparable wall times. The complete
  preconditioned chain took 7.7 hours on the paper's workstation, including
  its idle intervals. Fresh-factorization reference runs take additional time.

Run from this directory in MATLAB:

```matlab
setup_paths
addpath('experiments/rom_efficiency')
```

Output is written to `results/`, which is created on demand. Set the
`BRK_RESULTS` environment variable to an absolute directory before starting
MATLAB to write elsewhere. Saved result files are generated outputs and are
not part of the source package.

## Reproduce the manuscript

Run the experiment scripts below in the listed order. Each saves `.mat`
results through `results_dir`. The old fresh-factorization runs supply the
geometry panels, the cannula error study, and the fresh-LU comparison in the
text. The preconditioned chain supplies the final Tables 1--3, the compact
same-design audit, the compact runs under their own stopping rule, and the
objective/indicator histories.

```matlab
run_tests                % Appendix C: manufactured solutions and adjoint checks
eff_check_maps           % exact convection-map checks for the compact NSE ROM
run_sec51_cannula        % fresh-LU cannula, local ROM, cannula error study
run_sec52_body2d         % fresh-LU Re=40 and Re=200 body runs
run_sec53_body3d         % fresh-LU 3D body and local-ROM own-rule run
run_fig1_svd            % Figure 1 snapshot spectra
run_appC_validation      % Appendix C refinement, ablation, sampled errors
run_pc_baseline_chain    % final preconditioned FOM/local/compact comparisons
```

The chain's default phases are `cannula`, `re40`, `re200`, `body3d`, `appC`,
`ownrule`, and `audit`. It runs sequentially, waits five minutes before each
timed main leg and one minute before each timed Appendix C leg, and saves
`pc_*.mat` files. Pass a cell array of phase names to run a subset. Its
second argument, `true`, selects small grids and iteration counts for a
workflow check; those results do not reproduce the paper.

Generate the plots and numerical summary after the runs:

```matlab
fig1_svd
fig_sec51
fig_sec52
fig_sec53
fig_pc_main             % run last: replaces the older two-curve history panels
tables_appC
export_pc_numbers(results_dir('numbers.json'), results_dir())
```

The final call writes the numbers behind the principal tables and audits.
Compare its scientific quantities with `reference/published_numbers.json`;
wall times can vary by workstation. The generated figure files are written
under `results/figures/`. The manuscript's geometry and validation panels
come from `fig_sec51`, `fig_sec52`, and `fig_sec53`; its final three-method
history panels come from `fig_pc_main`.

## Published settings

The experiment drivers explicitly set `vol_project = false` (clip after the
volume-multiplier update) and `lsq = 'normal'` (normal equations for reduced
least squares). The driver defaults differ: `vol_project = true` and
`lsq = 'qr'`. The latter settings are used in the Appendix C refinement and
ablation runs. Running a driver directly with its defaults will not
necessarily recover a published trajectory.

The main chain first takes each preconditioned FOM to the manuscript's
stopping rule, then runs both ROMs for the same number of iterations. It uses
window/rank caps of 40/30 for the local ROM and 8/8 for the compact ROM. The
compact POD cutoff is `1e-7` for Navier--Stokes and `1e-10` for Stokes. A
maximum of two defect corrections is allowed after a failed residual screen
and two more after an objective trigger. The compact NSE reduced solve uses
at most two Gauss--Newton steps. Full residual checks, objective verification,
and the preconditioned full-order backstop remain in use.

Audits replay the compact trajectories with a same-design full-order solve
on sampled provisional pairs. They do not alter the design path, and audit
solve times are excluded from the reported optimization wall times.

## Code layout

| Directory | Purpose |
|---|---|
| `src/model` | MAC Brinkman flow operators, objectives, sensitivities, direct and preconditioned full-order solvers |
| `src/rom` | Local POD, LSPG, residual indicators, and dual-weighted estimate |
| `src/design` | Volume-constrained Allen--Cahn update |
| `src/drivers` | Cannula, 2D body, and 3D body optimizers |
| `experiments` | Published run drivers and numerical-summary exporter |
| `experiments/rom_efficiency` | Compact-ROM implementation and its direct dependencies |
| `figures` | Manuscript plots and Appendix C tables |
| `tests` | Code verification used in Appendix C |

The compact NSE driver retains helper routines for optional predictor and
sampled-solve branches because the executed driver references them. The
published chain does not enable those branches. Exploratory study drivers,
pilot outputs, and earlier compact-run reports have been removed.

## Provenance limit

The saved main-run log and `.mat` configurations match
`run_pc_baseline_chain` and the values in `reference/published_numbers.json`.
The code files used by the chain were last modified before its recorded start,
and the earlier 60-file source checksum list still matched this tree before
pruning. The original timed run did not record a Git commit or an attached
source archive, so exact byte-for-byte provenance of that historical run
cannot be established from the saved output alone. `reference/source_sha256.txt`
records the retained source supplied here for future releases.
