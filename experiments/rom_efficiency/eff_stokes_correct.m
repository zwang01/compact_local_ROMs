function [U,W,J,Jvisc,sens,eS,eA,rs,ra,g,etaJ] = ...
    eff_stokes_correct(par,S,rom,dvec,au,av,aw,dalpha,U,W,rs,Kref)
% One retained-factor defect correction for 2D/3D Stokes. The corrected
% state and adjoint are screened with the original full residual norms.
U = U + Kref.solve(rs);
rs = S.F - S.K0*U - dvec.*U;
if isempty(aw)
    [J,Jvisc,g] = brinkman_J_eval(par,S,au,av,U);
else
    [J,Jvisc,g] = brinkman3d_J_eval(par,S,au,av,aw,U);
end
ra = g - S.K0.'*W - dvec.*W;
W = W + Kref.solvet(ra);
ra = g - S.K0.'*W - dvec.*W;
eS = norm(rs)/max(norm(S.F),1e-14);
eA = norm(ra)/max(norm(g),1e-14);
if isempty(aw)
    sens = brinkman_sens_gather(par,S,dalpha,U,W);
else
    sens = brinkman3d_sens_gather(par,S,dalpha,U,W);
end
etaJ = eff_stokes_dwr(S,rom,dvec,W,rs,g);
end
