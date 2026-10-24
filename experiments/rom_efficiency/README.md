# Compact local-ROM implementation

The published compact method is called from
`experiments/run_pc_baseline_chain.m`. Its case drivers are:

- `eff_compact_cannula.m` for the two-dimensional Stokes cannula;
- `eff_shapeopt_goal2d.m` for the Navier--Stokes bodies at `Re=40` and
  `Re=200`, and for their preconditioned full-order references;
- `eff_compact_body3d.m` for Pironneau's three-dimensional Stokes body.

The Stokes drivers call `eff_stokes_correct.m` and `eff_stokes_dwr.m`.
The NSE driver calls `eff_rom_step.m`, `eff_enrich.m`,
`eff_defect_candidate.m`, `eff_defect_finish.m`, and
`eff_convection_maps.m`. It contains optional predictor and row-sampling
branches, so their helpers remain here to keep the executed driver intact;
the published configuration does not activate them. `eff_check_maps.m`
verifies the exact nonlinear convection, Jacobian, and transposed actions
used by the compact NSE reduced solve.

For each failed state/adjoint screen, the method first tries up to two
retained-factor defect-correction pairs. The corrected candidate must pass
the same calibrated full-residual screen. After a dual-weighted objective
trigger, it can try up to two further corrections. A full-order solve remains
the backstop. Accepted corrected pairs enter the eight-pair snapshot window.
The objective estimator enriches its adjoint trial space with the corrected
adjoint direction. The correction, full residual checks, and design update
still scale with the full-order dimension.

See the package [README](../../README.md) for the exact run order, published
settings, figure commands, and provenance limits.
