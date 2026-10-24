function eff_check_maps()
% Exactness check on mixed open/wall boundaries and nonzero ghost values.
rng(71);
par=struct('Nx',12,'Ny',6,'h',0.1,'nu',0.025,'x0',-0.2,'y0',-0.3, ...
    'openRight',mod((1:6)',2)==0,'openBottom',mod((1:12)',3)==0, ...
    'ubc',@(x,y) 0.3+0.2*x-0.1*y,'vbc',@(x,y) -0.1+0.1*x+0.3*y);
S=brinkman_base(par); C=eff_convection_maps(par,S);
U=randn(S.Ntot,1); V=randn(S.Ntot,4); Z=randn(S.Ntot,3);
[N,A]=brinkman_nse_convection(par,S,U);
x=C.X*U; y=C.Y*U; dx=C.Dx*U+C.gx; dy=C.Dy*U+C.gy;
Nm=zeros(S.Ntot,1); Nm(C.rows)=x.*dx+y.*dy;
AV=zeros(size(V)); AV(C.rows,:)=dx.*(C.X*V)+x.*(C.Dx*V)+dy.*(C.Y*V)+y.*(C.Dy*V);
Zr=Z(C.rows,:);
ATZ=C.X'*(dx.*Zr)+C.Dx'*(x.*Zr)+C.Y'*(dy.*Zr)+C.Dy'*(y.*Zr);
errors=[norm(Nm-N)/norm(N),norm(AV-A*V,'fro')/norm(A*V,'fro'), ...
    norm(ATZ-A'*Z,'fro')/norm(A'*Z,'fro')];
fprintf('Convection / Jacobian action / transpose action errors: %.3e %.3e %.3e\n',errors);
assert(all(errors<1e-12),'Exact stencil factorization disagrees with original convection.');
end
