% TEST_ADJOINT_CONSISTENCY  Is the transposed discrete Jacobian J_N' a
% consistent discretization of the continuous adjoint operator (Appendix C of
% the paper, code verification)
%   A* w = -nu Lap w - (u.grad) w + (grad u)' w + alpha w - grad pi_h,   div w = 0 ?
% Smooth u, w, pi and the smooth band alpha of test_mms_order.m are sampled on
% the MAC grid, J_N = K(phi) + N'(U) is assembled at the sampled u, J_N' is
% applied to the sampled (w, pi), and the result is compared with A* w at
% interior faces (two cells away from the boundary) and the continuity rows
% with -div w at interior cells. Expected: second order.
%   matlab -batch "setup_paths; test_adjoint_consistency"
nu = 0.1; amax = 2500; q = 0.1; w0 = 0.1; r0 = 0.3;
uex = @(x,y) 2*pi*sin(pi*x).*cos(2*pi*y);  vex = @(x,y) -pi*cos(pi*x).*sin(2*pi*y);
ux = @(x,y) 2*pi^2*cos(pi*x).*cos(2*pi*y);  uy = @(x,y) -4*pi^2*sin(pi*x).*sin(2*pi*y);
vx = @(x,y) pi^2*sin(pi*x).*sin(2*pi*y);    vy = @(x,y) -2*pi^2*cos(pi*x).*cos(2*pi*y);
% adjoint test fields (smooth, not related to u)
wu = @(x,y) sin(2*pi*x).*cos(pi*y);  wv = @(x,y) cos(pi*x).*sin(3*pi*y);  piex = @(x,y) sin(pi*x).*sin(2*pi*y);
wux = @(x,y) 2*pi*cos(2*pi*x).*cos(pi*y);  wuy = @(x,y) -pi*sin(2*pi*x).*sin(pi*y);
wvx = @(x,y) -pi*sin(pi*x).*sin(3*pi*y);   wvy = @(x,y) 3*pi*cos(pi*x).*cos(3*pi*y);
lapwu = @(x,y) -5*pi^2*wu(x,y);  lapwv = @(x,y) -10*pi^2*wv(x,y);
pix = @(x,y) pi*cos(pi*x).*sin(2*pi*y);   piy = @(x,y) 2*pi*sin(pi*x).*cos(2*pi*y);
phif = @(x,y) 0.5*(1 + tanh((sqrt((x-1).^2 + y.^2) - r0)/w0));
alphaf = @(x,y) amax*q*(1 - phif(x,y))./(q + phif(x,y));
% continuous adjoint operator (8), momentum components, with -grad pi_h
Au = @(x,y) -nu*lapwu(x,y) - (uex(x,y).*wux(x,y) + vex(x,y).*wuy(x,y)) + (ux(x,y).*wu(x,y) + vx(x,y).*wv(x,y)) + alphaf(x,y).*wu(x,y) - pix(x,y);
Av = @(x,y) -nu*lapwv(x,y) - (uex(x,y).*wvx(x,y) + vex(x,y).*wvy(x,y)) + (uy(x,y).*wu(x,y) + vy(x,y).*wv(x,y)) + alphaf(x,y).*wv(x,y) - piy(x,y);
divw = @(x,y) wux(x,y) + wvy(x,y);

fprintf('%5s %10s | %12s %12s %12s | %6s %6s %6s\n', 'N', 'h', 'mom-u err', 'mom-v err', 'div err', 'ord', 'ord', 'ord');
Nl = [32 64 128 256]; E = zeros(numel(Nl),3); hs = zeros(numel(Nl),1);
for g = 1:numel(Nl)
    N = Nl(g); par = struct('Nx',N,'Ny',N/2,'x0',0,'y0',-0.5,'h',2/N,'nu',nu,'alpha_max',amax,'qramp',q);
    par.ubc = uex; par.vbc = vex; par.openRight = false(par.Ny,1); par.openBottom = false(par.Nx,1);
    S = brinkman_base(par); h = par.h; Nx = N; Ny = N/2;
    xu = (0:Nx)*h; yu = -0.5 + ((1:Ny)-0.5)*h; [XU,YU] = ndgrid(xu,yu);
    xv = ((1:Nx)-0.5)*h; yv = -0.5 + (0:Ny)*h; [XV,YV] = ndgrid(xv,yv);
    xc = ((1:Nx)-0.5)*h; yc = -0.5 + ((1:Ny)-0.5)*h; [XC,YC] = ndgrid(xc,yc);
    U = [reshape(uex(XU,YU),[],1); reshape(vex(XV,YV),[],1); zeros(S.Np,1)];
    phi = phif(XC,YC); dvec = brinkman_coeffs(par, S, phi);
    [~, dN] = brinkman_nse_convection(par, S, U);
    JN = S.K0 + spdiags(dvec, 0, S.Ntot, S.Ntot) + dN;
    W = [reshape(wu(XU,YU),[],1); reshape(wv(XV,YV),[],1); reshape(piex(XC,YC),[],1)];
    R = JN.' * W;
    Ru = reshape(R(1:S.Nu), Nx+1, Ny); Rv = reshape(R(S.Nu+1:S.Nu+S.Nv), Nx, Ny+1); Rp = reshape(R(S.Nu+S.Nv+1:end), Nx, Ny);
    iu = 3:Nx-1; ju = 3:Ny-2; iv = 3:Nx-2; jv = 3:Ny-1; ic = 3:Nx-2; jc = 3:Ny-2;   % interior, away from boundary rows
    E(g,1) = max(abs(Ru(iu,ju) - Au(XU(iu,ju),YU(iu,ju))), [], 'all');
    E(g,2) = max(abs(Rv(iv,jv) - Av(XV(iv,jv),YV(iv,jv))), [], 'all');
    E(g,3) = max(abs(Rp(ic,jc) + divw(XC(ic,jc),YC(ic,jc))), [], 'all');     % continuity rows of J_N' give -div w
    hs(g) = h;
    if g == 1, o = nan(1,3); else, o = log(E(g-1,:)./E(g,:))/log(hs(g-1)/hs(g)); end
    fprintf('%5d %10.5f | %12.4e %12.4e %12.4e | %6.2f %6.2f %6.2f\n', N, h, E(g,:), o);
end
