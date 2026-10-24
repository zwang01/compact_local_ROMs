function [U, W, J, Jvisc, sens, eta_s, eta_a, nwt, rs, ra, dJdU, etaJ] = ...
    brinkman_nse_rom_step(par, S, rom, dvec, au, av, dalpha, U0)
% BRINKMAN_NSE_ROM_STEP  One reduced iteration of Algorithm 2 for 2D
% Navier-Stokes: LSPG (Gauss-Newton) reduced state, LSPG reduced adjoint,
% full-order residuals and the dual-weighted estimate, as brinkman_rom_step.m.
%   state  : a = argmin || F - K(phi) Psi_s a - N(Psi_s a) ||, by Gauss-Newton
%            with the exact Jacobian J_N = K(phi) + N'(U):
%            J_r = K0 Psi_s + dvec.*Psi_s + N'(U) Psi_s,  a <- a + J_r \ R;
%   adjoint: b = argmin || dJ/dU - J_N' Psi_a b ||.
% U0 is the warm start (the last accepted state), projected onto Psi_s.
% The convection blocks depend on the state, not on the design, so N'(U) Psi_s
% is formed at every Gauss-Newton step; eta_s is the relative full-order
% NONLINEAR residual. nwt returns the number of Gauss-Newton steps. The
% least-squares problems are solved as set by rom.lsq (rom_lsq.m).
Phi = rom.Phi;
a = Phi'*U0;                 % warm start (Phi orthonormal)
U = Phi*a;
KPhi = rom.K0Phi + dvec.*Phi;
nF = max(norm(S.F), 1e-14);
eta_s = inf; nwt = 0;

stale = true;                % are (Nc, dN) out of date with respect to U?
for it = 1:15
    [Nc, dN] = brinkman_nse_convection(par, S, U);
    stale = false;
    R  = S.F - (S.K0*U + dvec.*U + Nc);
    eta_new = norm(R)/nF;
    if it > 1 && (eta_new > eta_s || eta_new < 1e-12)
        eta_s = min(eta_s, eta_new); break        % converged or stalled
    end
    eta_s = eta_new;
    Jr = KPhi + dN*Phi;
    da = rom_lsq(Jr, R, rom.lsq);          % Gauss-Newton increment
    a = a + da;
    U = Phi*a;
    stale = true;
    nwt = it;
    if norm(da) <= 1e-12*max(norm(a),1), break, end
end
if stale                                   % (Nc, dN) at the final reduced state
    [Nc, dN] = brinkman_nse_convection(par, S, U);
end
rs = S.F - (S.K0*U + dvec.*U + Nc);        % full-order nonlinear residual
eta_s = norm(rs)/nF;

[J, Jvisc, dJdU] = brinkman_J_eval(par, S, au, av, U);

PhiA  = rom.PhiA;
JtPhi = rom.K0tPhiA + dvec.*PhiA + (dN.')*PhiA;
b = rom_lsq(JtPhi, dJdU, rom.lsq);
W = PhiA*b;
ra = dJdU - JtPhi*b;                       % full-order adjoint residual
eta_a = norm(ra)/max(norm(dJdU), 1e-14);

% dual-weighted estimate: adjoint re-solved in span{Psi_a, rs}
if strcmp(rom.lsq, 'normal')               % as in the runs of Sections 5.1-5.3
    rdir = rs/max(norm(rs),1e-30); inspan = false;
else                                       % residual orthogonalized against Psi_a
    [rdir, inspan] = dwr_direction(rs, PhiA);
end
if inspan
    etaJ = W.'*rs;
else
    JtP2 = [JtPhi, S.K0.'*rdir + dvec.*rdir + dN.'*rdir];
    b2   = rom_lsq(JtP2, dJdU, rom.lsq);
    etaJ = ([PhiA, rdir]*b2).'*rs;
end

sens = brinkman_sens_gather(par, S, dalpha, U, W);
end
