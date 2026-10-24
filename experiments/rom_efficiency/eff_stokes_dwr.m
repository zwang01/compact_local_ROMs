function etaJ = eff_stokes_dwr(S,rom,dvec,W,rs,g)
% Enriched adjoint LSPG solve containing the corrected adjoint and the
% residual direction, consistent with the original normal-equation mode.
T = rom.PhiA;
[extra,inspan] = dwr_direction(W,T);
if ~inspan, T = [T,extra]; end
if strcmp(rom.lsq,'normal')
    rdir = rs/max(norm(rs),1e-30); inspan = false;
else
    [rdir,inspan] = dwr_direction(rs,T);
end
if ~inspan, T = [T,rdir]; end
KtT = S.K0.'*T + dvec.*T;
b = rom_lsq(KtT,g,rom.lsq);
etaJ = (T*b).'*rs;
end
