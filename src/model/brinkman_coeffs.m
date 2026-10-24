function [dvec, au, av, dalpha] = brinkman_coeffs(par, S, phi)
% BRINKMAN_COEFFS  Design-dependent coefficients of the 2D Brinkman system:
% the RAMP coefficient [eq:ramp] on the cells,
%   alpha(phi) = alpha_max * q*(1 - phi)/(q + phi),   q = par.qramp (default 0.1),
% so that alpha(1) = 0 (fluid) and alpha(0) = alpha_max (solid); its face
% averages au, av; its cell derivative dalpha; and the Brinkman diagonal dvec
% on the momentum rows, K(phi) = S.K0 + diag(dvec).
Nx = par.Nx; Ny = par.Ny;
if isfield(par,'qramp'), q = par.qramp; else, q = 0.1; end
alpha  = par.alpha_max * q*(1 - phi)./(q + phi);
dalpha = -par.alpha_max * q*(1 + q)./(q + phi).^2;
au = zeros(Nx+1,Ny);
au(2:Nx,:) = 0.5*(alpha(1:Nx-1,:) + alpha(2:Nx,:));
au(1,:) = alpha(1,:); au(Nx+1,:) = alpha(Nx,:);
av = zeros(Nx,Ny+1);
av(:,2:Ny) = 0.5*(alpha(:,1:Ny-1) + alpha(:,2:Ny));
av(:,1) = alpha(:,1); av(:,Ny+1) = alpha(:,Ny);
dvec = [au(:).*S.momu(:); av(:).*S.momv(:); zeros(S.Np,1)];
end
