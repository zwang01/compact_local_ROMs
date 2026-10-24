function [U, W, J, Jvisc, dJdU, info] = brinkman3d_stokes_solve_pc(par, S, K, au, av, aw, Kref, U0, W0, opts)
% BRINKMAN3D_STOKES_SOLVE_PC  3D analogue of brinkman_stokes_solve_pc.m: the
% full-order state and adjoint at a design whose matrix differs from the
% retained factorization Kref only in the Brinkman diagonal, by pc_gmres
% preconditioned with K_ref^{-1} (K_ref^{-T} for the adjoint); refactors once
% the state solve has needed more than opts.refresh_its iterations.
if nargin < 8, U0 = []; end
if nargin < 9, W0 = []; end
[U, fl, ~, it] = pc_gmres(K, S.F, U0, Kref.solve, opts.tol, opts.maxit);
info = struct('ok', fl == 0, 'its', [it, 0], 'refreshed', false, 'Fref', Kref);
W = []; J = NaN; Jvisc = NaN; dJdU = [];
if ~info.ok, return, end
[J, Jvisc, dJdU] = brinkman3d_J_eval(par, S, au, av, aw, U);
[W, flw, ~, itw] = pc_gmres(K.', dJdU, W0, Kref.solvet, opts.tol, opts.maxit);
info.its(2) = itw;
info.ok = (flw == 0);
if info.ok && it > opts.refresh_its
    brk_count_fact;
    info.Fref = brk_lu(K, true); info.refreshed = true;
end
end
