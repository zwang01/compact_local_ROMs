function rom = eff_sample_prepare(rom,S,dvec,U)
% Leverage-score sampling of the local state residual. The sampled Newton
% solve uses a fixed number of rows. Full residual screening is still exact.
% Local RandStream makes repeats identical without altering application RNG.
C=rom.C; P=rom.Phi; x=C.X*U; y=C.Y*U;
dx=C.Dx*U+C.gx; dy=C.Dy*U+C.gy;
A=rom.K0Phi+dvec.*P;
A(C.rows,:)=A(C.rows,:)+dx.*rom.X+x.*rom.Dx+dy.*rom.Y+y.*rom.Dy;
[Q,~]=qr([A,S.F],0);
p=0.98*sum(Q.^2,2)/size(Q,2)+0.02/S.Ntot;
p=p/sum(p);
if ~isfield(rom,'sample_generation'), rom.sample_generation=0; end
rom.sample_generation=rom.sample_generation+1;
stream=RandStream('mt19937ar','Seed',1729+rom.sample_generation);
m=min(S.Ntot,rom.eff_sample_size);
edges=[0;cumsum(p)]; edges(end)=1;
si=sort(discretize(rand(stream,m,1),edges));
rom.si=si; rom.sw=1./sqrt(m*p(si));
loc=zeros(S.Ntot,1); loc(C.rows)=1:numel(C.rows);
rom.sactive=find(loc(si)>0); cr=loc(si(rom.sactive));
rom.sX=rom.X(cr,:); rom.sY=rom.Y(cr,:);
rom.sDx=rom.Dx(cr,:); rom.sDy=rom.Dy(cr,:);
rom.sgx=C.gx(cr); rom.sgy=C.gy(cr);
end
