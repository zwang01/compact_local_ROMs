function [dvec, au, av, aw, dalpha] = brinkman3d_coeffs(par, S, phi)
% BRINKMAN3D_COEFFS  3D analogue of brinkman_coeffs.m: RAMP coefficient on the
% cells, its averages on the three face families, its cell derivative, and the
% Brinkman diagonal dvec on the momentum rows, K(phi) = S.K0 + diag(dvec).
Nx = par.Nx; Ny = par.Ny; Nz = par.Nz;
if isfield(par,'qramp'), q = par.qramp; else, q = 0.1; end
alpha  = par.alpha_max * q*(1-phi)./(q+phi);
dalpha = -par.alpha_max * q*(1+q)./(q+phi).^2;
au = zeros(Nx+1,Ny,Nz);
au(2:Nx,:,:) = 0.5*(alpha(1:Nx-1,:,:)+alpha(2:Nx,:,:));
au(1,:,:) = alpha(1,:,:); au(Nx+1,:,:) = alpha(Nx,:,:);
av = zeros(Nx,Ny+1,Nz);
av(:,2:Ny,:) = 0.5*(alpha(:,1:Ny-1,:)+alpha(:,2:Ny,:));
av(:,1,:) = alpha(:,1,:); av(:,Ny+1,:) = alpha(:,Ny,:);
aw = zeros(Nx,Ny,Nz+1);
aw(:,:,2:Nz) = 0.5*(alpha(:,:,1:Nz-1)+alpha(:,:,2:Nz));
aw(:,:,1) = alpha(:,:,1); aw(:,:,Nz+1) = alpha(:,:,Nz);
dvec = [au(:).*S.momu(:); av(:).*S.momv(:); aw(:).*S.momw(:); zeros(S.Pn,1)];
end
