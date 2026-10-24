function [U, W, J, Jvisc, dJdU, info] = brinkman_stokes_solve_pc(par, S, K, au, av, Kref, U0, W0, opts)
% BRINKMAN_STOKES_SOLVE_PC  Full-order Stokes-Brinkman state and adjoint at a
% design whose matrix K = K_ref + diag(d) differs from the RETAINED
% factorization Kref only in the Brinkman diagonal, solved WITHOUT a new
% factorization by GMRES preconditioned with K_ref^{-1} (and K_ref^{-T} for
% the adjoint). K_ref^{-1}K = I + K_ref^{-1}diag(d) with d supported on the
% faces whose coefficient moved since K_ref was factored, so GMRES needs few
% iterations, each one product with K and one pair of triangular solves.
% This is the full-order solve of Algorithm 2 at a fallback or a verification.
%
%   U0, W0 : warm starts (the rejected or verified reduced pair); [] = zero
%   opts   : .tol         relative preconditioned residual (paper: 1e-10)
%            .maxit       GMRES iteration cap (unrestarted)
%            .refresh_its refactor once the state solve needed more (m_max = 8)
%   info   : .ok          both solves converged (otherwise the caller refactors)
%            .its         GMRES iterations [state, adjoint]
%            .refreshed   K was factored afresh; .Fref is the factorization to retain
if nargin < 7, U0 = []; end
if nargin < 8, W0 = []; end
[U, fl, ~, it] = pc_gmres(K, S.F, U0, Kref.solve, opts.tol, opts.maxit);
info = struct('ok', fl == 0, 'its', [it, 0], 'refreshed', false, 'Fref', Kref);
W = []; J = NaN; Jvisc = NaN; dJdU = [];
if ~info.ok, return, end

[J, Jvisc, dJdU] = brinkman_J_eval(par, S, au, av, U);
[W, flw, ~, itw] = pc_gmres(K.', dJdU, W0, Kref.solvet, opts.tol, opts.maxit);
info.its(2) = itw;
info.ok = (flw == 0);
% Refresh the retained factorization once the preconditioner has degraded
% (the design has drifted from the reference): the state solve needed more
% than m_max iterations.
if info.ok && it > opts.refresh_its
    brk_count_fact;
    info.Fref = brk_lu(K, true); info.refreshed = true;
end
end
