function [J,Jvisc,eta_s,eta_a,rs,ra,g,work] = eff_predict(par,S,rom,dvec,au,av,U,W)
% Screen a linear extrapolation before spending a reduced nonlinear solve.
% Both residuals are exact; no convection matrix is assembled.
C=rom.C;
x=C.X*U; y=C.Y*U; dx=C.Dx*U+C.gx; dy=C.Dy*U+C.gy;
rs=S.F-S.K0*U-dvec.*U; rs(C.rows)=rs(C.rows)-x.*dx-y.*dy;
[J,Jvisc,g]=brinkman_J_eval(par,S,au,av,U);
wr=W(C.rows);
ra=g-S.K0'*W-dvec.*W-C.X'*(dx.*wr)-C.Dx'*(x.*wr)-C.Y'*(dy.*wr)-C.Dy'*(y.*wr);
eta_s=norm(rs)/max(norm(S.F),1e-14); eta_a=norm(ra)/max(norm(g),1e-14);
work=struct('x',x,'y',y,'dx',dx,'dy',dy);
end
