function [J, sens, out] = brinkman_nse_solve_J(par, S, phi, U0)
% BRINKMAN_NSE_SOLVE_J  Full-order Navier-Stokes-Brinkman solve at the design
% phi by damped Newton with sparse LU, the objective, the exact discrete
% adjoint and the sensitivity dJ/dphi (Algorithm 1, lines 4-6, of the paper):
%   R(U) = K(phi)*U + N(U) - F = 0,     (K(phi) + N'(U))' W = dJ/dU.
% U0 warm-starts Newton; without it the Stokes solution is used (one extra
% factorization).
%
% par.adj_reuse_jac = false (Algorithm 1): the Jacobian is factored once more
%   at the converged state and the adjoint is solved with its transpose.
% par.adj_reuse_jac = true (used by the reduced runs of Algorithm 2): the
%   factorization of the last Newton step is reused for the adjoint; that
%   Jacobian belongs to the iterate one step before convergence, so the
%   adjoint is perturbed at the level of the last Newton increment, and a
%   full-order solve costs one factorization fewer.
% par.lu_explicit = true returns explicit LU factors (brk_lu.m), as needed
%   when the factorization is retained as a GMRES preconditioner.
%
% out: U, u, v, p, W, Jvisc, Jdrag, newton_its, res, and dJac, the factorization
% of the Jacobian used for the adjoint (the retained K_ref of Algorithm 2).

Ntot = S.Ntot;
[dvec, au, av, dalpha] = brinkman_coeffs(par, S, phi);
K = S.K0 + spdiags(dvec, 0, Ntot, Ntot);
lx = isfield(par, 'lu_explicit') && par.lu_explicit;

if nargin < 4 || isempty(U0)
    brk_count_fact;
    U = brk_lu(K, lx).solve(S.F);            % Stokes initial guess
else
    U = U0;
end

nF = max(norm(S.F), 1e-14);
newton_its = 0; dJac = [];
stale = true;                  % are (Nc, dN) out of date with respect to U?
for it = 1:25
    [Nc, dN] = brinkman_nse_convection(par, S, U);
    stale = false;
    R = K*U + Nc - S.F;
    nr = norm(R)/nF;
    if nr < 1e-11, break, end
    brk_count_fact;
    dJac = brk_lu(K + dN, lx);
    dU = dJac.solve(R);
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
if stale        % Newton left on the step-count limit: reassemble at the final U
    [Nc, dN] = brinkman_nse_convection(par, S, U);
end
res_final = norm(K*U + Nc - S.F)/nF;
if res_final > 1e-8
    warning('brinkman_nse_solve_J: Newton residual %.2e after %d iterations', res_final, newton_its);
end

Nx = par.Nx; Ny = par.Ny; Nu_ = S.Nu; Nv_ = S.Nv;
[J, Jvisc, dJdU, Jdrag] = brinkman_J_eval(par, S, au, av, U);
out = struct('U',U, 'u',reshape(U(1:Nu_),Nx+1,Ny), ...
    'v',reshape(U(Nu_+1:Nu_+Nv_),Nx,Ny+1), ...
    'p',reshape(U(Nu_+Nv_+1:end),Nx,Ny), ...
    'Jvisc',Jvisc, 'Jdrag',Jdrag, 'newton_its',newton_its, 'res',res_final);
if nargout < 2, sens = []; return, end

% exact discrete adjoint (K + dN(U))' W = dJ/dU
if ~(isfield(par,'adj_reuse_jac') && par.adj_reuse_jac && ~isempty(dJac))
    brk_count_fact;               % factor at the converged state (also when Newton took no step)
    dJac = brk_lu(K + dN, lx);
end
W = dJac.solvet(dJdU);
sens = brinkman_sens_gather(par, S, dalpha, U, W);
out.W = W;
out.dJac = dJac;
end
