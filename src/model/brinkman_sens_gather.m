function sens = brinkman_sens_gather(par, S, dalpha, U, W)
% BRINKMAN_SENS_GATHER  Discrete-adjoint sensitivity dJ/dphi per cell
% [eq:sens_flow] from a state U and adjoint W: on each momentum face the
% objective term (h^2 u^2 for the dissipation, h^2 u_x for the drag) minus the
% constraint term w.*u, gathered to the cells with the weights that average
% alpha to the faces, times alpha'(phi). The convection does not depend on phi,
% so the same formula serves the Stokes and the Navier-Stokes state.
Nx = par.Nx; Ny = par.Ny; h = par.h;
Nu = S.Nu; Nv = S.Nv;
u  = reshape(U(1:Nu), Nx+1, Ny);
v  = reshape(U(Nu+1:Nu+Nv), Nx, Ny+1);
wu = reshape(W(1:Nu), Nx+1, Ny);
wv = reshape(W(Nu+1:Nu+Nv), Nx, Ny+1);
if isfield(par,'objective') && strcmp(par.objective,'drag')
    su = h^2*u      - (wu.*u).*S.momu;     % dJ/dalpha_face = h^2*u_x for the drag
    sv =            - (wv.*v).*S.momv;
else
    su = h^2*u.^2   - (wu.*u).*S.momu;
    sv = h^2*v.^2   - (wv.*v).*S.momv;
end
sens = zeros(Nx,Ny);
sens(1:Nx-1,:) = sens(1:Nx-1,:) + 0.5*su(2:Nx,:);
sens(2:Nx,:)   = sens(2:Nx,:)   + 0.5*su(2:Nx,:);
sens(1,:)      = sens(1,:)      + su(1,:);
sens(Nx,:)     = sens(Nx,:)     + su(Nx+1,:);
sens(:,1:Ny-1) = sens(:,1:Ny-1) + 0.5*sv(:,2:Ny);
sens(:,2:Ny)   = sens(:,2:Ny)   + 0.5*sv(:,2:Ny);
sens(:,1)      = sens(:,1)      + sv(:,1);
sens(:,Ny)     = sens(:,Ny)     + sv(:,Ny+1);
sens = dalpha .* sens;
end
