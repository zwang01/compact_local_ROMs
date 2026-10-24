function [U,W,J,Jvisc,sens,eta_s,eta_a,nwt,rs,ra,dJdU,etaJ] = ...
    eff_rom_step(par,S,rom,dvec,au,av,dalpha,U0)
% Exact reduced convection/Jacobian actions, without global sparse Jacobian
% assembly. Full residual screening and the original enriched-adjoint DWR
% space are retained. 'mapped' uses original stopping; 'capped' changes only
% the Gauss-Newton iteration budget and increment tolerance.
if strcmp(rom.eff_variant,'baseline')
    [U,W,J,Jvisc,sens,eta_s,eta_a,nwt,rs,ra,dJdU,etaJ]= ...
        brinkman_nse_rom_step(par,S,rom,dvec,au,av,dalpha,U0);
    return
end
C=rom.C; rows=C.rows; P=rom.Phi; PA=rom.PhiA;
a=P'*U0; KP=rom.K0Phi+dvec.*P; nF=max(norm(S.F),1e-14);
eta_s=inf; nwt=0;
if strcmp(rom.eff_variant,'sampled')
    [a,nwt]=eff_sample_state(rom,S,dvec,a);
else
for it=1:rom.eff_max_newton
    x=rom.X*a; y=rom.Y*a; dx=rom.Dx*a+C.gx; dy=rom.Dy*a+C.gy;
    R=S.F-KP*a; R(rows)=R(rows)-x.*dx-y.*dy;
    eta_new=norm(R)/nF;
    if it>1 && (eta_new>eta_s || eta_new<1e-12), break; end
    eta_s=eta_new;
    Jr=KP;
    Jr(rows,:)=Jr(rows,:)+dx.*rom.X+x.*rom.Dx+dy.*rom.Y+y.*rom.Dy;
    da=rom_lsq(Jr,R,rom.lsq); a=a+da; nwt=it;
    if norm(da)<=rom.eff_newton_tol*max(norm(a),1), break; end
end
end
U=P*a;
x=rom.X*a; y=rom.Y*a; dx=rom.Dx*a+C.gx; dy=rom.Dy*a+C.gy;
rs=S.F-KP*a; rs(rows)=rs(rows)-x.*dx-y.*dy; eta_s=norm(rs)/nF;
[J,Jvisc,dJdU]=brinkman_J_eval(par,S,au,av,U);
At=rom.K0tPhiA+dvec.*PA+jt(PA(rows,:));
if strcmp(rom.lsq,'normal')
    G=At'*At; rhs=At'*dJdU; b=G\rhs;
else
    b=rom_lsq(At,dJdU,rom.lsq);
end
W=PA*b; ra=dJdU-At*b; eta_a=norm(ra)/max(norm(dJdU),1e-14);
if strcmp(rom.lsq,'normal')
    rdir=rs/max(norm(rs),1e-30); inspan=false;
else
    [rdir,inspan]=dwr_direction(rs,PA);
end
if inspan
    etaJ=W'*rs;
else
    z=S.K0'*rdir+dvec.*rdir+jt(rdir(rows));
    if strcmp(rom.lsq,'normal')
        % Reuse G and rhs: exactly the same bordered normal system as
        % [At,z]'*[At,z], avoiding its second tall matrix multiplication.
        c=At'*z; b2=[G,c;c',z'*z]\[rhs;z'*dJdU];
    else
        b2=rom_lsq([At,z],dJdU,rom.lsq);
    end
    etaJ=(PA'*rs)'*b2(1:end-1)+(rdir'*rs)*b2(end);
end
sens=brinkman_sens_gather(par,S,dalpha,U,W);
    function out=jt(V)
        out=C.X'*(dx.*V)+C.Dx'*(x.*V)+C.Y'*(dy.*V)+C.Dy'*(y.*V);
    end
end
