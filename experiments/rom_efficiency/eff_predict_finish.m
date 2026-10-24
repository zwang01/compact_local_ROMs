function [sens,etaJ] = eff_predict_finish(par,S,rom,dvec,dalpha,U,W,rs,g,work)
% Accepted predictions retain the original full enriched-adjoint DWR solve.
% We save the state/adjoint reduced solves, not the objective safeguard.
C=rom.C; PA=rom.PhiA; x=work.x; y=work.y; dx=work.dx; dy=work.dy;
A=rom.K0tPhiA+dvec.*PA+jt(PA(C.rows,:));
if strcmp(rom.lsq,'normal')
    rdir=rs/max(norm(rs),1e-30); inspan=false;
else
    [rdir,inspan]=dwr_direction(rs,PA);
end
if inspan
    b=rom_lsq(A,g,rom.lsq); etaJ=(PA'*rs)'*b;
else
    z=S.K0'*rdir+dvec.*rdir+jt(rdir(C.rows));
    b=rom_lsq([A,z],g,rom.lsq);
    etaJ=(PA'*rs)'*b(1:end-1)+(rdir'*rs)*b(end);
end
sens=brinkman_sens_gather(par,S,dalpha,U,W);
    function out=jt(V)
        out=C.X'*(dx.*V)+C.Dx'*(x.*V)+C.Y'*(dy.*V)+C.Dy'*(y.*V);
    end
end
