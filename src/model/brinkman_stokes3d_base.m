function S = brinkman_stokes3d_base(par)
% BRINKMAN_STOKES3D_BASE  Design-independent part of the 3D MAC discretization
% of the Brinkman-penalized Stokes equations on the box
% [x0, x0+Nx*h] x [y0, y0+Ny*h] x [z0, z0+Nz*h] (paper, Section 5.3):
%   -nu*Lap(u) + alpha(phi)*u + grad p = 0,   div u = 0,
% K(phi) = S.K0 + diag(dvec) as in 2D (brinkman_base.m).
% Boundary conditions of the duct:
%   x = x0          : Dirichlet inflow u = par.uin(y,z), v = w = 0
%   x = x0 + Nx*h   : do-nothing outflow, nu*du/dx - p = 0, dv/dx = dw/dx = 0
%   y, z walls      : no slip (zero wall values, so the energy shift e0 vanishes)
% par.sym = true: y = y0 and z = z0 are SYMMETRY planes (the quarter duct
%   y, z >= 0 of the paper): tangential components see an even ghost there and
%   carry no wall-edge energy row.
% Grid: p at cell centres; u at x-faces (Nx+1,Ny,Nz); v at y-faces (Nx,Ny+1,Nz);
% w at z-faces (Nx,Ny,Nz+1). U = [u(:); v(:); w(:); p(:)].

Nx=par.Nx; Ny=par.Ny; Nz=par.Nz; h=par.h; nu=par.nu;
y0=par.y0; z0=par.z0;
sym = isfield(par,'sym') && par.sym;
if sym, sg = -1; else, sg = +1; end
Un=(Nx+1)*Ny*Nz; Vn=Nx*(Ny+1)*Nz; Wn=Nx*Ny*(Nz+1); Pn=Nx*Ny*Nz;
Ntot=Un+Vn+Wn+Pn;
iu=@(i,j,k) i+(j-1)*(Nx+1)+(k-1)*(Nx+1)*Ny;
iv=@(i,j,k) Un+i+(j-1)*Nx+(k-1)*Nx*(Ny+1);
iw=@(i,j,k) Un+Vn+i+(j-1)*Nx+(k-1)*Nx*Ny;
ip=@(i,j,k) Un+Vn+Wn+i+(j-1)*Nx+(k-1)*Nx*Ny;
yc=@(j) y0+(j-0.5)*h; zc=@(k) z0+(k-0.5)*h;

II=zeros(10*Ntot,1); JJ=II; VV=II; nz_=0; F=zeros(Ntot,1);
momu=false(Nx+1,Ny,Nz); momv=false(Nx,Ny+1,Nz); momw=false(Nx,Ny,Nz+1);
    function add(r,c,val), nz_=nz_+1; II(nz_)=r; JJ(nz_)=c; VV(nz_)=val; end

vh=nu/h^2;
% ---- u rows ----
for k=1:Nz, for j=1:Ny, for i=1:Nx+1
    r=iu(i,j,k);
    if i==1
        add(r,r,1); F(r)=par.uin(yc(j),zc(k));
    elseif i==Nx+1
        add(r,r,nu/h); add(r,iu(Nx,j,k),-nu/h); add(r,ip(Nx,j,k),-1);
    else
        momu(i,j,k)=true; d=6*vh;
        add(r,iu(i-1,j,k),-vh); add(r,iu(i+1,j,k),-vh);
        if j>1,  add(r,iu(i,j-1,k),-vh); else, d=d+sg*vh; end
        if j<Ny, add(r,iu(i,j+1,k),-vh); else, d=d+vh; end
        if k>1,  add(r,iu(i,j,k-1),-vh); else, d=d+sg*vh; end
        if k<Nz, add(r,iu(i,j,k+1),-vh); else, d=d+vh; end
        add(r,r,d);
        add(r,ip(i,j,k),1/h); add(r,ip(i-1,j,k),-1/h);
    end
end, end, end
% ---- v rows ----
for k=1:Nz, for j=1:Ny+1, for i=1:Nx
    r=iv(i,j,k);
    if j==1 || j==Ny+1
        add(r,r,1);
    else
        momv(i,j,k)=true; d=6*vh;
        add(r,iv(i,j-1,k),-vh); add(r,iv(i,j+1,k),-vh);
        if i>1,  add(r,iv(i-1,j,k),-vh); else, d=d+vh; end
        if i<Nx, add(r,iv(i+1,j,k),-vh); else, d=d-vh; end   % open right
        if k>1,  add(r,iv(i,j,k-1),-vh); else, d=d+sg*vh; end
        if k<Nz, add(r,iv(i,j,k+1),-vh); else, d=d+vh; end
        add(r,r,d);
        add(r,ip(i,j,k),1/h); add(r,ip(i,j-1,k),-1/h);
    end
end, end, end
% ---- w rows ----
for k=1:Nz+1, for j=1:Ny, for i=1:Nx
    r=iw(i,j,k);
    if k==1 || k==Nz+1
        add(r,r,1);
    else
        momw(i,j,k)=true; d=6*vh;
        add(r,iw(i,j,k-1),-vh); add(r,iw(i,j,k+1),-vh);
        if i>1,  add(r,iw(i-1,j,k),-vh); else, d=d+vh; end
        if i<Nx, add(r,iw(i+1,j,k),-vh); else, d=d-vh; end   % open right
        if j>1,  add(r,iw(i,j-1,k),-vh); else, d=d+sg*vh; end
        if j<Ny, add(r,iw(i,j+1,k),-vh); else, d=d+vh; end
        add(r,r,d);
        add(r,ip(i,j,k),1/h); add(r,ip(i,j,k-1),-1/h);
    end
end, end, end
% ---- continuity ----
for k=1:Nz, for j=1:Ny, for i=1:Nx
    r=ip(i,j,k);
    add(r,iu(i+1,j,k),1/h); add(r,iu(i,j,k),-1/h);
    add(r,iv(i,j+1,k),1/h); add(r,iv(i,j,k),-1/h);
    add(r,iw(i,j,k+1),1/h); add(r,iw(i,j,k),-1/h);
end, end, end
K0=sparse(II(1:nz_),JJ(1:nz_),VV(1:nz_),Ntot,Ntot);

% ---- energy operators (all first differences; e0 = 0, see header) ----
EI=zeros(12*Ntot,1); EJ=EI; EV=EI; ez=0; erow=0;
    function eadd(cols,vals)
        erow=erow+1;
        for q=1:numel(cols), ez=ez+1; EI(ez)=erow; EJ(ez)=cols(q); EV(ez)=vals(q); end
    end
% u derivatives
for k=1:Nz, for j=1:Ny
    for i=1:Nx, eadd([iu(i+1,j,k),iu(i,j,k)],[1/h,-1/h]); end        % du/dx
end, end
for k=1:Nz, for i=1:Nx+1
    for j=1:Ny-1, eadd([iu(i,j+1,k),iu(i,j,k)],[1/h,-1/h]); end      % du/dy
    if ~sym, eadd(iu(i,1,k),2/h); end; eadd(iu(i,Ny,k),2/h);         % y walls (none at a symmetry plane)
    end, end
for j=1:Ny, for i=1:Nx+1
    for k=1:Nz-1, eadd([iu(i,j,k+1),iu(i,j,k)],[1/h,-1/h]); end      % du/dz
    if ~sym, eadd(iu(i,j,1),2/h); end; eadd(iu(i,j,Nz),2/h);         % z walls
end, end
% v derivatives
for k=1:Nz, for i=1:Nx
    for j=1:Ny, eadd([iv(i,j+1,k),iv(i,j,k)],[1/h,-1/h]); end        % dv/dy
end, end
for k=1:Nz, for j=1:Ny+1
    for i=1:Nx-1, eadd([iv(i+1,j,k),iv(i,j,k)],[1/h,-1/h]); end      % dv/dx
    eadd(iv(1,j,k),2/h);                                             % inlet wall
end, end
for j=1:Ny+1, for i=1:Nx
    for k=1:Nz-1, eadd([iv(i,j,k+1),iv(i,j,k)],[1/h,-1/h]); end      % dv/dz
    if ~sym, eadd(iv(i,j,1),2/h); end; eadd(iv(i,j,Nz),2/h);
end, end
% w derivatives
for j=1:Ny, for i=1:Nx
    for k=1:Nz, eadd([iw(i,j,k+1),iw(i,j,k)],[1/h,-1/h]); end        % dw/dz
end, end
for k=1:Nz+1, for j=1:Ny
    for i=1:Nx-1, eadd([iw(i+1,j,k),iw(i,j,k)],[1/h,-1/h]); end      % dw/dx
    eadd(iw(1,j,k),2/h);
end, end
for k=1:Nz+1, for i=1:Nx
    for j=1:Ny-1, eadd([iw(i,j+1,k),iw(i,j,k)],[1/h,-1/h]); end      % dw/dy
    if ~sym, eadd(iw(i,1,k),2/h); end; eadd(iw(i,Ny,k),2/h);
end, end
E=sparse(EI(1:ez),EJ(1:ez),EV(1:ez),erow,Un+Vn+Wn);

S=struct('K0',K0,'F',F,'E',E,'momu',momu,'momv',momv,'momw',momw, ...
    'Un',Un,'Vn',Vn,'Wn',Wn,'Pn',Pn,'Ntot',Ntot);
end
