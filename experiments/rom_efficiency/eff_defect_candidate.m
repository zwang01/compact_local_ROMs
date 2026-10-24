function [U,W,J,Jvisc,eS,eA,rs,ra,g,work]=eff_defect_candidate(par,S,rom,dvec,au,av,U,W,rs,Kref)
% One inexact Newton state correction followed by an adjoint correction
% at the corrected state. Two retained-factor inverse applications total.
U=U+Kref.solve(rs);
[J,Jvisc,eS,~,rs,ra,g,work]=eff_predict(par,S,rom,dvec,au,av,U,W);
dW=Kref.solvet(ra); W=W+dW;
C=rom.C; V=dW(C.rows);
ra=ra-S.K0'*dW-dvec.*dW-C.X'*(work.dx.*V)-C.Dx'*(work.x.*V) ...
    -C.Y'*(work.dy.*V)-C.Dy'*(work.y.*V);
eA=norm(ra)/max(norm(g),1e-14);
end
