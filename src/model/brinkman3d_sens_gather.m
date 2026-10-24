function sens = brinkman3d_sens_gather(par, S, dalpha, U, W)
% BRINKMAN3D_SENS_GATHER  3D analogue of brinkman_sens_gather.m: discrete-adjoint
% sensitivity dJ/dphi per cell from a state U and adjoint W.
Nx = par.Nx; Ny = par.Ny; Nz = par.Nz; h = par.h;
Un = S.Un; Vn = S.Vn; Wn = S.Wn;
u = reshape(U(1:Un),Nx+1,Ny,Nz); v = reshape(U(Un+1:Un+Vn),Nx,Ny+1,Nz); w = reshape(U(Un+Vn+1:Un+Vn+Wn),Nx,Ny,Nz+1);
wu = reshape(W(1:Un),Nx+1,Ny,Nz); wv = reshape(W(Un+1:Un+Vn),Nx,Ny+1,Nz); ww = reshape(W(Un+Vn+1:Un+Vn+Wn),Nx,Ny,Nz+1);
if isfield(par,'objective') && strcmp(par.objective,'drag')
    su = h^3*u    - (wu.*u).*S.momu;
    sv =          - (wv.*v).*S.momv;
    sw =          - (ww.*w).*S.momw;
else
    su = h^3*u.^2 - (wu.*u).*S.momu;
    sv = h^3*v.^2 - (wv.*v).*S.momv;
    sw = h^3*w.^2 - (ww.*w).*S.momw;
end
sens = zeros(Nx,Ny,Nz);
sens(1:Nx-1,:,:) = sens(1:Nx-1,:,:) + 0.5*su(2:Nx,:,:);
sens(2:Nx,:,:)   = sens(2:Nx,:,:)   + 0.5*su(2:Nx,:,:);
sens(1,:,:) = sens(1,:,:) + su(1,:,:);   sens(Nx,:,:) = sens(Nx,:,:) + su(Nx+1,:,:);
sens(:,1:Ny-1,:) = sens(:,1:Ny-1,:) + 0.5*sv(:,2:Ny,:);
sens(:,2:Ny,:)   = sens(:,2:Ny,:)   + 0.5*sv(:,2:Ny,:);
sens(:,1,:) = sens(:,1,:) + sv(:,1,:);   sens(:,Ny,:) = sens(:,Ny,:) + sv(:,Ny+1,:);
sens(:,:,1:Nz-1) = sens(:,:,1:Nz-1) + 0.5*sw(:,:,2:Nz);
sens(:,:,2:Nz)   = sens(:,:,2:Nz)   + 0.5*sw(:,:,2:Nz);
sens(:,:,1) = sens(:,:,1) + sw(:,:,1);   sens(:,:,Nz) = sens(:,:,Nz) + sw(:,:,Nz+1);
sens = dalpha.*sens;
end
