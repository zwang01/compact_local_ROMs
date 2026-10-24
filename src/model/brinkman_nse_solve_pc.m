function [J, sens, out, info] = brinkman_nse_solve_pc(par, S, phi, U0, W0, Kref, opts)
% BRINKMAN_NSE_SOLVE_PC  Full-order Navier-Stokes-Brinkman solve at a new
% design WITHOUT a new factorization (fallback and verification of
% Algorithm 2): Newton-Krylov in which every Newton system
%     (K(phi) + dN(U)) dU = R(U)
% is solved by GMRES preconditioned with the factorization Kref RETAINED from
% the last factored solve, and the adjoint likewise with the transposed pair.
% The preconditioned Jacobian differs from the identity by the Brinkman
% diagonal that moved since Kref was formed plus the change in the convection
% Jacobian. Returns the same triple as brinkman_nse_solve_J plus info;
% info.ok = false when a GMRES call missed opts.tol within opts.maxit
% iterations or Newton did not converge, and the caller then factors.
%
%   U0, W0 : warm starts (the rejected or verified reduced state and adjoint)
%   opts   : .tol, .maxit as in brinkman_stokes_solve_pc; .refresh_its (m_max = 20):
%            refactor at the converged state once the solve needed more GMRES
%            iterations in total (all Newton steps and the adjoint)
%   info   : .ok, .its (total GMRES iterations), .newton_its, .res, .refreshed
Ntot = S.Ntot;
[dvec, au, av, dalpha] = brinkman_coeffs(par, S, phi);
K = S.K0 + spdiags(dvec, 0, Ntot, Ntot);
nF = max(norm(S.F), 1e-14);

U = U0;
info = struct('ok', true, 'its', 0, 'newton_its', 0, 'res', NaN, 'refreshed', false);
J = NaN; sens = []; out = struct();
stale = true; newton_its = 0; Nc = []; dN = [];
for it = 1:25
    [Nc, dN] = brinkman_nse_convection(par, S, U);
    stale = false;
    R = K*U + Nc - S.F;
    nr = norm(R)/nF;
    if nr < 1e-11, break, end
    % inexact Newton (Eisenstat-Walker): early steps need not be solved tightly
    tol_k = max(opts.tol, min(1e-2, 0.1*nr));
    [dU, fl, ~, git] = pc_gmres(K + dN, R, [], Kref.solve, tol_k, opts.maxit);
    info.its = info.its + git;
    if fl ~= 0, info.ok = false; return, end
    % damped update: backtrack if the residual does not decrease
    lam = 1;
    for bt = 1:6
        Ut = U - lam*dU;
        Nt = brinkman_nse_convection(par, S, Ut);
        if norm(K*Ut + Nt - S.F)/nF < nr, break, end
        lam = lam/2;
    end
    U = U - lam*dU;
    stale = true;
    newton_its = it;
end
if stale
    [Nc, dN] = brinkman_nse_convection(par, S, U);
end
res_final = norm(K*U + Nc - S.F)/nF;
info.newton_its = newton_its; info.res = res_final;
if res_final > 1e-8, info.ok = false; return, end

Nx = par.Nx; Ny = par.Ny; Nu_ = S.Nu; Nv_ = S.Nv;
[J, Jvisc, dJdU, Jdrag] = brinkman_J_eval(par, S, au, av, U);
[W, flw, ~, gitw] = pc_gmres((K + dN).', dJdU, W0, Kref.solvet, opts.tol, opts.maxit);
info.its = info.its + gitw;
if flw ~= 0, info.ok = false; J = NaN; return, end
sens = brinkman_sens_gather(par, S, dalpha, U, W);
out = struct('U',U, 'u',reshape(U(1:Nu_),Nx+1,Ny), ...
    'v',reshape(U(Nu_+1:Nu_+Nv_),Nx,Ny+1), ...
    'p',reshape(U(Nu_+Nv_+1:end),Nx,Ny), ...
    'Jvisc',Jvisc, 'Jdrag',Jdrag, 'newton_its',newton_its, 'res',res_final, ...
    'W',W, 'dJac',Kref);    % the retained factorization stays the reference...
if info.its > opts.refresh_its
    brk_count_fact;          % ...unless the preconditioner has degraded: refactor at the converged state
    out.dJac = brk_lu(K + dN, true); info.refreshed = true;
end
end
