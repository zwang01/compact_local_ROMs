function [a,nwt]=eff_sample_state(rom,S,dvec,a)
% All Gauss-Newton operations here use m sampled equations and r unknowns;
% complexity of this inner solve is independent of the FOM dimension.
si=rom.si; active=rom.sactive; w=rom.sw;
KP=rom.K0Phi(si,:)+dvec(si).*rom.Phi(si,:); F=S.F(si);
eta=inf; nwt=0;
for it=1:rom.eff_max_newton
    x=rom.sX*a; y=rom.sY*a; dx=rom.sDx*a+rom.sgx; dy=rom.sDy*a+rom.sgy;
    R=F-KP*a; R(active)=R(active)-x.*dx-y.*dy;
    eta_new=norm(w.*R)/max(norm(w.*F),1e-14);
    if it>1 && (eta_new>eta || eta_new<1e-12), break; end
    eta=eta_new;
    A=KP; A(active,:)=A(active,:)+dx.*rom.sX+x.*rom.sDx+dy.*rom.sY+y.*rom.sDy;
    da=rom_lsq(w.*A,w.*R,rom.lsq); a=a+da; nwt=it;
    if norm(da)<=rom.eff_newton_tol*max(norm(a),1), break; end
end
end
