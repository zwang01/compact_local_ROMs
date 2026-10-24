function [sens,etaJ]=eff_defect_finish(par,S,rom,dvec,dalpha,U,W,rs,g,work)
% The DWR trial space explicitly contains the corrected adjoint, as well as
% the previous adjoint basis and the current residual. This is a full LSPG
% enriched-adjoint solve, not a substitute based only on W'*rs.
[extra,inspan]=dwr_direction(W,rom.PhiA);
if ~inspan
    rom.PhiA=[rom.PhiA,extra]; rom.K0tPhiA=[rom.K0tPhiA,S.K0'*extra];
end
[sens,etaJ]=eff_predict_finish(par,S,rom,dvec,dalpha,U,W,rs,g,work);
end
