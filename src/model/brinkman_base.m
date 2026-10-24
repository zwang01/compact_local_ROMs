function S = brinkman_base(par)
% BRINKMAN_BASE  Design-independent part of the 2D MAC discretization of the
% Brinkman-penalized Stokes equations on the box [x0, x0+Nx*h] x [y0, y0+Ny*h],
%
%   -nu*Lap(u) + alpha(phi)*u + grad p = 0,   div u = 0,
%
% with alpha(phi) the RAMP coefficient of brinkman_coeffs.m (fluid phi = 1,
% solid phi = 0). The Brinkman term is the ONLY design dependence and it is
% diagonal, so the saddle-point matrix is
%   K(phi) = S.K0 + diag(dvec),   dvec = face-averaged alpha on the momentum rows,
% and the reduced operators K(phi)*Psi = K0*Psi + dvec.*Psi are formed in O(n r).
% The Navier-Stokes convection is added by brinkman_nse_convection.m.
%
% Grid: p(i,j) at cell centres (Nx x Ny); u(i,j) on x-faces ((Nx+1) x Ny);
% v(i,j) on y-faces (Nx x (Ny+1)). Unknowns U = [u(:); v(:); p(:)].
%
% Boundary conditions:
%  * Dirichlet velocity (values par.ubc(x,y), par.vbc(x,y)) except on the open
%    (do-nothing, traction-free) segments
%      par.openRight(j)  : cell rows j open at x = x0 + Nx*h,
%      par.openBottom(i) : cell columns i open at y = y0,
%    where nu*du_n/dn - p = 0 and du_t/dn = 0;
%  * one pressure is pinned to zero if no open segment exists.
%
% Also builds the energy operators of the dissipation objective,
%   J = h^2 * ( nu*||E*[u;v] + e0||^2 + sum(a_u.*u.^2) + sum(a_v.*v.^2) ),
% where E collects all first differences of u and v (wall values through the
% ghost convention of the momentum stencil; no rows at open walls).

Nx = par.Nx; Ny = par.Ny; h = par.h; nu = par.nu;
x0 = par.x0; y0 = par.y0;
openR = par.openRight(:);    % Ny x 1 logical
openB = par.openBottom(:);   % Nx x 1 logical

Nu = (Nx+1)*Ny; Nv = Nx*(Ny+1); Np = Nx*Ny; Ntot = Nu+Nv+Np;
iu = @(i,j) i + (j-1)*(Nx+1);
iv = @(i,j) Nu + i + (j-1)*Nx;
ip = @(i,j) Nu + Nv + i + (j-1)*Nx;

xu = @(i) x0 + (i-1)*h;   yu = @(j) y0 + (j-0.5)*h;
xv = @(i) x0 + (i-0.5)*h; yv = @(j) y0 + (j-1)*h;
xR = x0 + Nx*h; yT = y0 + Ny*h;

% open-ness seen from tangential unknowns
obu = @(i) (i>1 && openB(i-1)) || (i<=Nx && openB(min(i,Nx)));  % u-face i vs bottom
orv = @(j) (j>1 && openR(j-1)) || (j<=Ny && openR(min(j,Ny)));  % v-face j vs right

II = zeros(12*Ntot,1); JJ = II; VV = II; nz = 0;
F  = zeros(Ntot,1);
momu = false(Nx+1,Ny); momv = false(Nx,Ny+1);
    function add(r,c,val)
        nz = nz+1; II(nz)=r; JJ(nz)=c; VV(nz)=val;
    end

% ---------------- u rows ----------------
for j = 1:Ny
    for i = 1:Nx+1
        r = iu(i,j);
        if i == 1
            add(r,r,1); F(r) = par.ubc(x0, yu(j));
        elseif i == Nx+1 && openR(j)
            % do-nothing: nu*(u(Nx+1)-u(Nx))/h - p(Nx,j) = 0
            add(r, r, nu/h); add(r, iu(Nx,j), -nu/h); add(r, ip(Nx,j), -1);
        elseif i == Nx+1
            add(r,r,1); F(r) = par.ubc(xR, yu(j));
        else
            momu(i,j) = true;
            d = 4*nu/h^2;
            add(r, iu(i-1,j), -nu/h^2); add(r, iu(i+1,j), -nu/h^2);
            if j > 1, add(r, iu(i,j-1), -nu/h^2);
            elseif obu(i), d = d - nu/h^2;                       % open: du/dy = 0
            else, d = d + nu/h^2; F(r) = F(r) + 2*nu/h^2*par.ubc(xu(i), y0);
            end
            if j < Ny, add(r, iu(i,j+1), -nu/h^2);
            else, d = d + nu/h^2; F(r) = F(r) + 2*nu/h^2*par.ubc(xu(i), yT);
            end
            add(r, r, d);
            add(r, ip(i,j), 1/h); add(r, ip(i-1,j), -1/h);
        end
    end
end

% ---------------- v rows ----------------
for j = 1:Ny+1
    for i = 1:Nx
        r = iv(i,j);
        if j == 1 && openB(i)
            % do-nothing: nu*(v(i,2)-v(i,1))/h - p(i,1) = 0
            add(r, iv(i,2), nu/h); add(r, iv(i,1), -nu/h); add(r, ip(i,1), -1);
        elseif j == 1
            add(r,r,1); F(r) = par.vbc(xv(i), y0);
        elseif j == Ny+1
            add(r,r,1); F(r) = par.vbc(xv(i), yT);
        else
            momv(i,j) = true;
            d = 4*nu/h^2;
            add(r, iv(i,j-1), -nu/h^2); add(r, iv(i,j+1), -nu/h^2);
            if i > 1, add(r, iv(i-1,j), -nu/h^2);
            else, d = d + nu/h^2; F(r) = F(r) + 2*nu/h^2*par.vbc(x0, yv(j));
            end
            if i < Nx, add(r, iv(i+1,j), -nu/h^2);
            elseif orv(j), d = d - nu/h^2;                       % open: dv/dx = 0
            else, d = d + nu/h^2; F(r) = F(r) + 2*nu/h^2*par.vbc(xR, yv(j));
            end
            add(r, r, d);
            add(r, ip(i,j), 1/h); add(r, ip(i,j-1), -1/h);
        end
    end
end

% ---------------- continuity rows ----------------
for j = 1:Ny
    for i = 1:Nx
        r = ip(i,j);
        add(r, iu(i+1,j), 1/h); add(r, iu(i,j), -1/h);
        add(r, iv(i,j+1), 1/h); add(r, iv(i,j), -1/h);
    end
end

K0 = sparse(II(1:nz), JJ(1:nz), VV(1:nz), Ntot, Ntot);

% pin one pressure if the domain has no open boundary
if ~any(openR) && ~any(openB)
    r = ip(1,1);
    K0(r,:) = 0; K0(r,r) = 1; F(r) = 0;
end

% ---------------- energy (dissipation) operators ----------------
EI = zeros(8*(Nu+Nv),1); EJ = EI; EV = EI; ez = 0; e0v = []; erow = 0;
    function eadd(cols, vals, shift)
        erow = erow + 1;
        for kk = 1:numel(cols)
            ez = ez+1; EI(ez)=erow; EJ(ez)=cols(kk); EV(ez)=vals(kk);
        end
        e0v(erow,1) = shift; %#ok<AGROW>
    end
% du/dx at cell centres
for j = 1:Ny, for i = 1:Nx
    eadd([iu(i+1,j), iu(i,j)], [1/h, -1/h], 0);
end, end
% du/dy at interior horizontal edges and wall edges
for i = 1:Nx+1
    for j = 1:Ny-1
        eadd([iu(i,j+1), iu(i,j)], [1/h, -1/h], 0);
    end
    if ~obu(i)  % bottom wall edge: (u1 - (2ub-u1))/h
        eadd(iu(i,1), 2/h, -2*par.ubc(xu(i),y0)/h);
    end
    eadd(iu(i,Ny), 2/h, -2*par.ubc(xu(i),yT)/h);   % top wall edge
end
% dv/dy at cell centres
for j = 1:Ny, for i = 1:Nx
    eadd([iv(i,j+1), iv(i,j)], [1/h, -1/h], 0);
end, end
% dv/dx at interior vertical edges and wall edges
for j = 1:Ny+1
    for i = 1:Nx-1
        eadd([iv(i+1,j), iv(i,j)], [1/h, -1/h], 0);
    end
    eadd(iv(1,j), 2/h, -2*par.vbc(x0,yv(j))/h);    % left wall edge
    if ~orv(j)  % right wall edge
        eadd(iv(Nx,j), 2/h, -2*par.vbc(xR,yv(j))/h);
    end
end
E  = sparse(EI(1:ez), EJ(1:ez), EV(1:ez), erow, Nu+Nv);
e0 = e0v;

S = struct('K0',K0,'F',F,'E',E,'e0',e0,'momu',momu,'momv',momv, ...
    'Nu',Nu,'Nv',Nv,'Np',Np,'Ntot',Ntot);
end
