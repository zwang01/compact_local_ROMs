% TEST_MMS_ORDER  Code verification of the penalized solver by a manufactured
% solution (Table C.1 of the paper, rows 2-3), with the Brinkman term active.
% All-Dirichlet box (0,2) x (-1/2,1/2), nu = 0.1, the divergence-free velocity
% and pressure
%   u = 2 pi sin(pi x) cos(2 pi y),  v = -pi cos(pi x) sin(2 pi y),
%   p = cos(pi x) cos(2 pi y),
% the smooth phase field phi = (1 + tanh((|x - (1,0)| - 0.3)/0.1))/2, and the
% RAMP coefficient with alpha_max = 2.5e4*nu. The forcing
% f = -nu Lap u [+ (u.grad)u] + alpha u + grad p is added to the momentum rows
% and the discrete problem is solved as in the optimization (face-averaged
% alpha, pinned pressure). Prints the velocity, mean-free pressure and objective
% errors (J = int nu|grad u|^2 + alpha|u|^2, exact value by fine quadrature)
% with their observed orders, h = 1/16 ... 1/128, for Stokes and Navier-Stokes.
%   matlab -batch "setup_paths; test_mms_order"
nu = 0.1; w = 0.1; r0 = 0.3; amax = 2.5e4*nu;
uex = @(x,y) 2*pi*sin(pi*x).*cos(2*pi*y);
vex = @(x,y) -pi*cos(pi*x).*sin(2*pi*y);
pex = @(x,y) cos(pi*x).*cos(2*pi*y);
ux = @(x,y) 2*pi^2*cos(pi*x).*cos(2*pi*y);  uy = @(x,y) -4*pi^2*sin(pi*x).*sin(2*pi*y);
vx = @(x,y) pi^2*sin(pi*x).*sin(2*pi*y);    vy = @(x,y) -2*pi^2*cos(pi*x).*cos(2*pi*y);
px = @(x,y) -pi*sin(pi*x).*cos(2*pi*y);     py = @(x,y) -2*pi*cos(pi*x).*sin(2*pi*y);
phif = @(x,y) 0.5*(1 + tanh((sqrt((x-1).^2 + y.^2) - r0)/w));
q = 0.1; alphaf = @(x,y) amax*q*(1 - phif(x,y))./(q + phif(x,y));
Nlist = [32 64 128 256];

% exact objective by fine midpoint quadrature (well beyond the grids)
Mq = 2048; xq = (0.5:Mq-0.5)/Mq*2; yq = -0.5 + (0.5:Mq/2-0.5)/(Mq/2);
[XQ, YQ] = ndgrid(xq, yq); dAq = (2/Mq)*(1/(Mq/2));
gradsq = ux(XQ,YQ).^2 + uy(XQ,YQ).^2 + vx(XQ,YQ).^2 + vy(XQ,YQ).^2;
usq = uex(XQ,YQ).^2 + vex(XQ,YQ).^2;
Jex = sum(nu*gradsq(:) + alphaf(XQ(:),YQ(:)).*usq(:))*dAq;

for model = {'stokes', 'nse'}
    conv = strcmp(model{1}, 'nse');
    fprintf('\n=== %s, alpha_max = %g, nu = %g, band width %g (J_exact = %.8e) ===\n', model{1}, amax, nu, w, Jex);
    fprintf('%5s %10s %10s %10s %10s %10s | %6s %6s %6s %6s\n', 'N', 'h', 'max|du|', 'L2(u,v)', 'L2(p)', '|dJ|/J', 'ord', 'ord', 'ord', 'ord');
    E = zeros(numel(Nlist), 4); hs = zeros(numel(Nlist), 1);
    for g = 1:numel(Nlist)
        N = Nlist(g);
        par = struct('Nx', N, 'Ny', N/2, 'x0', 0, 'y0', -0.5, 'h', 2/N, 'nu', nu, 'alpha_max', amax, 'qramp', q);
        par.ubc = uex; par.vbc = vex;
        par.openRight = false(par.Ny,1); par.openBottom = false(par.Nx,1);
        par.objective = 'dissipation';
        S = brinkman_base(par); h = par.h;
        xu = par.x0 + (0:N)*h;        yu = par.y0 + ((1:N/2)-0.5)*h;   [XU, YU] = ndgrid(xu, yu);
        xv = par.x0 + ((1:N)-0.5)*h;  yv = par.y0 + (0:N/2)*h;         [XV, YV] = ndgrid(xv, yv);
        xc = par.x0 + ((1:N)-0.5)*h;  yc = par.y0 + ((1:N/2)-0.5)*h;   [XC, YC] = ndgrid(xc, yc);
        fu = 5*pi^2*nu*uex(XU,YU) + alphaf(XU,YU).*uex(XU,YU) + px(XU,YU);
        fv = 5*pi^2*nu*vex(XV,YV) + alphaf(XV,YV).*vex(XV,YV) + py(XV,YV);
        if conv
            fu = fu + uex(XU,YU).*ux(XU,YU) + vex(XU,YU).*uy(XU,YU);
            fv = fv + uex(XV,YV).*vx(XV,YV) + vex(XV,YV).*vy(XV,YV);
        end
        S.F(1:S.Nu) = S.F(1:S.Nu) + fu(:).*S.momu(:);
        S.F(S.Nu+1:S.Nu+S.Nv) = S.F(S.Nu+1:S.Nu+S.Nv) + fv(:).*S.momv(:);
        phi = phif(XC, YC);
        if conv, [J, ~, out] = brinkman_nse_solve_J(par, S, phi);
        else,    [J, ~, out] = brinkman_stokes_solve_J(par, S, phi); end
        eu = out.u - uex(XU,YU); ev = out.v - vex(XV,YV);
        ep = out.p - pex(XC,YC); ep = ep - mean(ep(:));            % pressure up to a constant
        E(g,1) = max(max(abs(eu(:))), max(abs(ev(:))));
        E(g,2) = sqrt(h^2*(sum(eu(:).^2) + sum(ev(:).^2)));
        E(g,3) = sqrt(h^2*sum(ep(:).^2));
        E(g,4) = abs(J - Jex)/abs(Jex);
        hs(g) = h;
        if g == 1, o = nan(1,4); else, o = log(E(g-1,:)./E(g,:))/log(hs(g-1)/hs(g)); end
        fprintf('%5d %10.5f %10.3e %10.3e %10.3e %10.3e | %6.2f %6.2f %6.2f %6.2f\n', N, h, E(g,:), o);
    end
end
