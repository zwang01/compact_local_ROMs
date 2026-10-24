function [N, dN] = brinkman_nse_convection(par, S, U)
% BRINKMAN_NSE_CONVECTION  Convection term N(U) = (u.grad)u of the steady
% Navier-Stokes equations (advective form, central differences on the MAC
% grid) and its exact Jacobian dN = N'(U). Contributions live only on the
% interior momentum rows (S.momu, S.momv); continuity and boundary rows carry
% none. Ghost values follow brinkman_base.m: tangential Dirichlet walls use
% ghost = 2*u_bc - u_in, open boundaries a zero normal derivative.
% N is exactly quadratic in U, and dN is skipped when only N is requested.
Nx = par.Nx; Ny = par.Ny; h = par.h;
x0 = par.x0; y0 = par.y0; xR = x0 + Nx*h; yT = y0 + Ny*h;
Nu_ = S.Nu; Nv_ = S.Nv;
u = reshape(U(1:Nu_), Nx+1, Ny);
v = reshape(U(Nu_+1:Nu_+Nv_), Nx, Ny+1);
openR = par.openRight(:); openB = par.openBottom(:);
iu = @(i,j) i + (j-1)*(Nx+1);
iv = @(i,j) Nu_ + i + (j-1)*Nx;

N = zeros(S.Ntot,1);
needJ = nargout > 1;          % the Jacobian is the bulk of the work: skip it
Ic = {}; Jc = {}; Vc = {};    % when only the residual vector is wanted
    function put(R, C, V)
        if ~needJ, return, end
        Ic{end+1} = R(:); Jc{end+1} = C(:); Vc{end+1} = V(:);
    end

% ---------------- u-momentum faces: i = 2..Nx, j = 1..Ny ----------------
ii = (2:Nx)';
xu = x0 + (ii-1)*h;
obu = false(Nx-1,1);                       % bottom open at this face?
for m = 1:Nx-1
    obu(m) = openB(ii(m)-1) || openB(ii(m));
end
uc   = u(2:Nx, :);
dudx = (u(3:Nx+1, :) - u(1:Nx-1, :)) / (2*h);
u_dn = zeros(Nx-1, Ny); u_dn(:, 2:Ny) = u(2:Nx, 1:Ny-1);
gb = 2*par.ubc(xu, y0*ones(Nx-1,1)) - u(2:Nx, 1);     % bottom wall ghost
ub1 = u(2:Nx, 1); gb(obu) = ub1(obu);                 % open: du/dy = 0
u_dn(:, 1) = gb;
u_up = zeros(Nx-1, Ny); u_up(:, 1:Ny-1) = u(2:Nx, 2:Ny);
u_up(:, Ny) = 2*par.ubc(xu, yT*ones(Nx-1,1)) - u(2:Nx, Ny);  % top wall
dudy = (u_up - u_dn) / (2*h);
vbar = 0.25*( v(1:Nx-1, 1:Ny) + v(1:Nx-1, 2:Ny+1) ...
            + v(2:Nx,   1:Ny) + v(2:Nx,   2:Ny+1) );
Nu_val = uc.*dudx + vbar.*dudy;
[IIg, JJg] = ndgrid(2:Nx, 1:Ny);
R = iu(IIg, JJg);
N(R(:)) = Nu_val(:);

% Jacobian, u-part
Dd = dudx;                                             % d/du(i,j) from uc
Dd(:,1)  = Dd(:,1)  + (~obu).*vbar(:,1)/(2*h) - obu.*vbar(:,1)/(2*h);
Dd(:,Ny) = Dd(:,Ny) - vbar(:,Ny)/(2*h);
put(R, R, Dd);
put(R, iu(IIg+1, JJg),  uc/(2*h));
put(R, iu(IIg-1, JJg), -uc/(2*h));
put(R(:,1:Ny-1), iu(IIg(:,1:Ny-1), JJg(:,1:Ny-1)+1),  vbar(:,1:Ny-1)/(2*h));
put(R(:,2:Ny),   iu(IIg(:,2:Ny),   JJg(:,2:Ny)-1),   -vbar(:,2:Ny)/(2*h));
put(R, iv(IIg-1, JJg),   0.25*dudy);
put(R, iv(IIg-1, JJg+1), 0.25*dudy);
put(R, iv(IIg,   JJg),   0.25*dudy);
put(R, iv(IIg,   JJg+1), 0.25*dudy);

% ---------------- v-momentum faces: i = 1..Nx, j = 2..Ny ----------------
jj = (2:Ny)';
yv = y0 + (jj-1)*h;
orv = false(Ny-1,1);                       % right open at this face?
for m = 1:Ny-1
    orv(m) = openR(jj(m)-1) || openR(jj(m));
end
vc   = v(:, 2:Ny);
dvdy = (v(:, 3:Ny+1) - v(:, 1:Ny-1)) / (2*h);
v_lt = zeros(Nx, Ny-1); v_lt(2:Nx, :) = v(1:Nx-1, 2:Ny);
v_lt(1, :) = 2*par.vbc(x0*ones(Ny-1,1), yv)' - v(1, 2:Ny);    % left wall
v_rt = zeros(Nx, Ny-1); v_rt(1:Nx-1, :) = v(2:Nx, 2:Ny);
grt = 2*par.vbc(xR*ones(Ny-1,1), yv)' - v(Nx, 2:Ny);          % right wall
vN = v(Nx, 2:Ny); grt(orv') = vN(orv');                       % open: dv/dx = 0
v_rt(Nx, :) = grt;
dvdx = (v_rt - v_lt) / (2*h);
ubar = 0.25*( u(1:Nx, 1:Ny-1) + u(1:Nx, 2:Ny) ...
            + u(2:Nx+1, 1:Ny-1) + u(2:Nx+1, 2:Ny) );
Nv_val = ubar.*dvdx + vc.*dvdy;
[IIg2, JJg2] = ndgrid(1:Nx, 2:Ny);
R2 = iv(IIg2, JJg2);
N(R2(:)) = Nv_val(:);

% Jacobian, v-part
Dd2 = dvdy;
Dd2(1,:)  = Dd2(1,:)  + ubar(1,:)/(2*h);               % left wall ghost
Dd2(Nx,:) = Dd2(Nx,:) - (~orv').*ubar(Nx,:)/(2*h) + orv'.*ubar(Nx,:)/(2*h);
put(R2, R2, Dd2);
put(R2, iv(IIg2, JJg2+1),  vc/(2*h));
put(R2, iv(IIg2, JJg2-1), -vc/(2*h));
put(R2(2:Nx,:),   iv(IIg2(2:Nx,:)-1,   JJg2(2:Nx,:)),   -ubar(2:Nx,:)/(2*h));
put(R2(1:Nx-1,:), iv(IIg2(1:Nx-1,:)+1, JJg2(1:Nx-1,:)),  ubar(1:Nx-1,:)/(2*h));
put(R2, iu(IIg2,   JJg2-1), 0.25*dvdx);
put(R2, iu(IIg2,   JJg2),   0.25*dvdx);
put(R2, iu(IIg2+1, JJg2-1), 0.25*dvdx);
put(R2, iu(IIg2+1, JJg2),   0.25*dvdx);

if needJ
    % The triplet list has ~820k entries collapsing onto ~456k nonzeros, and
    % sparse() spends most of its time accumulating those duplicates: 39 ms
    % against 8 ms from a deduplicated list at production size. The row and
    % column indices depend only on the grid and the open-boundary flags, not
    % on U, so the accumulation map is built once and reused.
    persistent ck cI cJ cmap
    key = {Nx, Ny, openR, openB, S.Ntot};
    if isempty(ck) || ~isequal(ck, key)
        Iall = vertcat(Ic{:}); Jall = vertcat(Jc{:});
        [IJ, ~, cmap] = unique([Iall, Jall], 'rows');
        cI = IJ(:,1); cJ = IJ(:,2); ck = key;
    end
    dN = sparse(cI, cJ, accumarray(cmap, vertcat(Vc{:})), S.Ntot, S.Ntot);
else
    dN = [];
end
end
