function [U, W, J, Jvisc, sens, eta_s, eta_a, rs, ra, dJdU, etaJ] = brinkman3d_rom_step(par, S, rom, dvec, au, av, aw, dalpha)
% BRINKMAN3D_ROM_STEP  3D analogue of brinkman_rom_step.m: LSPG reduced state
% and adjoint of the Stokes-Brinkman system, full-order residuals, Euclidean
% indicators, and the dual-weighted estimate etaJ of J(U) - J(U_r).
KPhi = rom.K0Phi + dvec.*rom.Phi;
a = rom_lsq(KPhi, S.F, rom.lsq);
U = rom.Phi*a;
rs = S.F - KPhi*a;
eta_s = norm(rs)/max(norm(S.F),1e-14);

[J, Jvisc, dJdU] = brinkman3d_J_eval(par, S, au, av, aw, U);

KtPhi = rom.K0tPhiA + dvec.*rom.PhiA;
b = rom_lsq(KtPhi, dJdU, rom.lsq);
W = rom.PhiA*b;
ra = dJdU - KtPhi*b;
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

sens = brinkman3d_sens_gather(par, S, dalpha, U, W);
end
