function [U, W, J, Jvisc, sens, eta_s, eta_a, rs, ra, dJdU, etaJ] = brinkman_rom_step(par, S, rom, dvec, au, av, dalpha)
% BRINKMAN_ROM_STEP  One reduced iteration of Algorithm 2 for 2D Stokes:
% least-squares Petrov-Galerkin (LSPG) state and adjoint, their full-order
% residuals, the dual-weighted estimate of the objective error, and the
% sensitivity at the reduced pair.
%   state  : a = argmin || F - K(phi) Psi_s a ||,        U = Psi_s a,
%   adjoint: b = argmin || dJ/dU - K(phi)' Psi_a b ||,   W = Psi_a b,
% with K(phi) Psi_s = K0 Psi_s + dvec.*Psi_s formed in O(n r) from the products
% of brinkman_rom_enrich.m. The least-squares problems are solved by QR, or
% by the normal equations when rom.lsq = 'normal' (rom_lsq.m).
% Returns the residuals rs, ra and the Euclidean indicators [eq:indicator_flow]
%   eta_s = ||rs||/||F||,  eta_a = ||ra||/||dJ/dU||,
% and etaJ, the dual-weighted estimate of J(U) - J(U_r) [eq:dwr_flow]: the
% adjoint re-solved in span{Psi_a, rs} (one more product K'*r_dir and an
% (r_a+1)-column least squares), etaJ = w_+' rs.
KPhi = rom.K0Phi + dvec.*rom.Phi;
a = rom_lsq(KPhi, S.F, rom.lsq);
U = rom.Phi*a;
rs = S.F - KPhi*a;                         % full-order state residual
eta_s = norm(rs)/max(norm(S.F),1e-14);

[J, Jvisc, dJdU] = brinkman_J_eval(par, S, au, av, U);

KtPhi = rom.K0tPhiA + dvec.*rom.PhiA;
b = rom_lsq(KtPhi, dJdU, rom.lsq);
W = rom.PhiA*b;
ra = dJdU - KtPhi*b;                       % full-order adjoint residual
eta_a = norm(ra)/max(norm(dJdU),1e-14);

if strcmp(rom.lsq, 'normal')               % as in the runs of Sections 5.1-5.3
    rdir = rs/max(norm(rs),1e-30); inspan = false;
else                                       % residual orthogonalized against Psi_a
    [rdir, inspan] = dwr_direction(rs, rom.PhiA);
end
if inspan
    etaJ = W.'*rs;
else
    KtP2 = [KtPhi, S.K0.'*rdir + dvec.*rdir];
    b2   = rom_lsq(KtP2, dJdU, rom.lsq);
    etaJ = ([rom.PhiA, rdir]*b2).'*rs;
end

sens = brinkman_sens_gather(par, S, dalpha, U, W);
end
