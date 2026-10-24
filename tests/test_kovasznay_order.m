% TEST_KOVASZNAY_ORDER  Code verification of the Navier-Stokes solver on the
% Kovasznay solution at Re = 40 (Table C.1 of the paper, first row): all fluid
% (alpha = 0), all-Dirichlet box (0,2) x (-1/2,1/2), h = 1/10 ... 1/160.
% Prints the maximum and L2 velocity errors, the L2 error of the mean-free
% pressure, the relative error of the dissipation nu*int|grad u|^2 (exact value
% in closed form), their observed orders, and the Newton steps from the Stokes
% initial guess.
%   matlab -batch "setup_paths; test_kovasznay_order"
Re = 40; nu = 1/Re;
lam = 1/(2*nu) - sqrt(1/(4*nu^2) + 4*pi^2);
uex = @(x,y) 1 - exp(lam*x).*cos(2*pi*y);
vex = @(x,y) (lam/(2*pi))*exp(lam*x).*sin(2*pi*y);
pex = @(x,y) 0.5*(1 - exp(2*lam*x));
Nlist = [20 40 80 160 320]; ng = numel(Nlist);
E = zeros(ng, 4); hs = zeros(ng,1); nwt = zeros(ng,1); res = zeros(ng,1);
for g = 1:ng
    par = struct(); par.Nx = Nlist(g); par.Ny = Nlist(g)/2;
    par.x0 = 0; par.y0 = -0.5; par.h = 2/par.Nx;
    par.nu = nu; par.alpha_max = 1;                  % phi = 1 everywhere: alpha = 0
    par.ubc = @(x,y) uex(x,y); par.vbc = @(x,y) vex(x,y);
    par.openRight = false(par.Ny,1); par.openBottom = false(par.Nx,1);
    S = brinkman_base(par); phi = ones(par.Nx, par.Ny);
    [~, ~, out] = brinkman_nse_solve_J(par, S, phi);
    h = par.h; hs(g) = h;
    xu = par.x0 + (0:par.Nx)*h;       yu = par.y0 + ((1:par.Ny)-0.5)*h;   [XU, YU] = ndgrid(xu, yu);
    xv = par.x0 + ((1:par.Nx)-0.5)*h; yv = par.y0 + (0:par.Ny)*h;         [XV, YV] = ndgrid(xv, yv);
    xp = par.x0 + ((1:par.Nx)-0.5)*h; yp = par.y0 + ((1:par.Ny)-0.5)*h;   [XP, YP] = ndgrid(xp, yp);
    du = out.u - uex(XU, YU); dv = out.v - vex(XV, YV);
    dp = out.p - pex(XP, YP); dp = dp - mean(dp(:));
    % dissipation nu*int |grad u|^2 over (0,2) x (-1/2,1/2): the y-integrals of cos^2 and sin^2 are 1/2
    Jex = nu*(exp(4*lam) - 1)/(2*lam) * 0.5*(2*lam^2 + 4*pi^2 + lam^4/(4*pi^2));
    E(g,:) = [max(abs([du(:); dv(:)])), sqrt(h^2*(sum(du(:).^2) + sum(dv(:).^2))), sqrt(h^2*sum(dp(:).^2)), abs(out.Jvisc - Jex)/Jex];
    nwt(g) = out.newton_its; res(g) = out.res;
end
fprintf('Kovasznay flow, Re = %g (lambda = %.4f), domain (0,2) x (-1/2,1/2), alpha = 0\n', Re, lam);
fprintf('%5s %9s | %11s %11s %11s %11s | %6s %6s %6s %6s | %6s %9s\n', 'N', 'h', 'max|du|', 'L2(u,v)', 'L2(p)', '|dJ|/J', 'ord', 'ord', 'ord', 'ord', 'Newton', 'residual');
for g = 1:ng
    if g == 1, o = nan(1,4); else, o = log(E(g-1,:)./E(g,:))/log(hs(g-1)/hs(g)); end
    fprintf('%5d %9.5f | %11.3e %11.3e %11.3e %11.3e | %6.2f %6.2f %6.2f %6.2f | %6d %9.1e\n', Nlist(g), hs(g), E(g,:), o, nwt(g), res(g));
end
