function C = eff_convection_maps(par, S)
% Exact stencil factorization N(U)(rows) = (X U).*(Dx U+gx) +
% (Y U).*(Dy U+gy). Includes the boundary ghost contributions. No sampling.
nx=par.Nx; ny=par.Ny; h=par.h; n=S.Ntot;
iu=@(i,j) i+(j-1)*(nx+1);
iv=@(i,j) S.Nu+i+(j-1)*nx;
[I,J]=ndgrid(2:nx,1:ny); ru=iu(I,J); mu=numel(I); a=(1:mu)';
[P,Q]=ndgrid(1:nx,2:ny); rv=iv(P,Q); mv=numel(P); b=mu+(1:mv)'; m=mu+mv;
C.rows=[ru(:);rv(:)];
cols=[ru(:);reshape(iu(P,Q-1),[],1);reshape(iu(P,Q),[],1); ...
    reshape(iu(P+1,Q-1),[],1);reshape(iu(P+1,Q),[],1)];
C.X=sparse([a;b;b;b;b],cols,[ones(mu,1);0.25*ones(4*mv,1)],m,n);
cols=[reshape(iv(I-1,J),[],1);reshape(iv(I-1,J+1),[],1); ...
    reshape(iv(I,J),[],1);reshape(iv(I,J+1),[],1);rv(:)];
C.Y=sparse([a;a;a;a;b],cols,[0.25*ones(4*mu,1);ones(mv,1)],m,n);

% x derivatives: ordinary u stencil; v left/right ghost conventions.
vl=-ones(size(P))/(2*h); vr=-vl;
vl(P==1)=1/(2*h);
openR=reshape(par.openRight(Q(:)-1)|par.openRight(Q(:)),size(Q));
vr(P==nx & ~openR)=-1/(2*h);
C.Dx=sparse([a;a;b;b], ...
    [reshape(iu(I+1,J),[],1);reshape(iu(I-1,J),[],1); ...
     reshape(iv(max(P-1,1),Q),[],1);reshape(iv(min(P+1,nx),Q),[],1)], ...
    [ones(mu,1)/(2*h);-ones(mu,1)/(2*h);vl(:);vr(:)],m,n);
gx=zeros(size(P));
left=P==1; right=P==nx & ~openR;
gx(left)=-par.vbc(par.x0*ones(nnz(left),1),par.y0+(Q(left)-1)*h)/h;
gx(right)=par.vbc((par.x0+nx*h)*ones(nnz(right),1),par.y0+(Q(right)-1)*h)/h;
C.gx=[zeros(mu,1);gx(:)];

% y derivatives: u bottom/top ghosts; ordinary v stencil.
openB=reshape(par.openBottom(I(:)-1)|par.openBottom(I(:)),size(I));
ud=-ones(size(I))/(2*h); uu=-ud;
ud(J==1 & ~openB)=1/(2*h); uu(J==ny)=-1/(2*h);
C.Dy=sparse([a;a;b;b], ...
    [reshape(iu(I,max(J-1,1)),[],1);reshape(iu(I,min(J+1,ny)),[],1); ...
     reshape(iv(P,Q+1),[],1);reshape(iv(P,Q-1),[],1)], ...
    [ud(:);uu(:);ones(mv,1)/(2*h);-ones(mv,1)/(2*h)],m,n);
gy=zeros(size(I)); bottom=J==1 & ~openB; top=J==ny;
gy(bottom)=-par.ubc(par.x0+(I(bottom)-1)*h,par.y0*ones(nnz(bottom),1))/h;
gy(top)=par.ubc(par.x0+(I(top)-1)*h,(par.y0+ny*h)*ones(nnz(top),1))/h;
C.gy=[gy(:);zeros(mv,1)];
end
