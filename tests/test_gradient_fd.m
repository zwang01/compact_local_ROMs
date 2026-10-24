% TEST_GRADIENT_FD  Central finite-difference check of the discrete-adjoint
% sensitivity dJ/dphi [eq:sens_flow] (Appendix C of the paper, code
% verification), on small grids with a smooth phase field in (0,1):
%   (1) 2D Stokes, dissipation          (2) 2D Navier-Stokes, dissipation
%   (3) 2D Navier-Stokes, drag          (4) 3D Stokes, drag (as in Section 5.3)
% Five cells per case (band and fluid cells in 3D), step 1e-6; prints the
% relative differences and fails if one exceeds 1e-4.
%   matlab -batch "setup_paths; test_gradient_fd"
delta = 1e-6; worst = 0;

% ---------- 2D: channel (-1/2, 3/2) x (-1/2, 1/2), open outlet ----------
par = struct('Nx', 24, 'Ny', 12, 'x0', -0.5, 'y0', -0.5, 'h', 2/24, 'nu', 5e-3, 'alpha_max', 500);
par.ubc = @(x,y) (abs(x+0.5) < 1e-12).*(0.05 - 0.2*y.^2);
par.vbc = @(x,y) 0*x;
par.openRight = true(par.Ny,1); par.openBottom = false(par.Nx,1);
S = brinkman_base(par);
[Xc, Yc] = ndgrid(par.x0+((1:par.Nx)-0.5)*par.h, par.y0+((1:par.Ny)-0.5)*par.h);
phi = 0.55 + 0.35*sin(3*Xc).*cos(2*Yc);
cases = {'Stokes, dissipation', 'dissipation', @brinkman_stokes_solve_J; ...
         'Navier-Stokes, dissipation', 'dissipation', @brinkman_nse_solve_J; ...
         'Navier-Stokes, drag', 'drag', @brinkman_nse_solve_J};
rng(1); cells = randperm(par.Nx*par.Ny, 5);
for c = 1:size(cases,1)
    par.objective = cases{c,2}; solve = cases{c,3};
    [J0, sens] = solve(par, S, phi);
    fprintf('\n2D %s (J = %.6e)\n', cases{c,1}, J0);
    for i = cells
        e = zeros(size(phi)); e(i) = delta;
        fd = (solve(par, S, phi+e) - solve(par, S, phi-e))/(2*delta);
        rel = abs(fd - sens(i))/max(abs(fd), 1e-14); worst = max(worst, rel);
        fprintf('  cell %4d: FD %+.8e  adjoint %+.8e  rel %.1e\n', i, fd, sens(i), rel);
    end
end

% ---------- 3D: box (-1/2, 3/2) x (-1/2, 1/2)^2, smooth sphere ----------
N = 16;
p3 = struct('Nx', N, 'Ny', N/2, 'Nz', N/2, 'h', 2/N, 'x0', -0.5, 'y0', -0.5, 'z0', -0.5, ...
            'nu', 1, 'alpha_max', 2.5e4, 'qramp', 0.1);
p3.uin = @(y,z) 0.05*(1-4*y.^2).*(1-4*z.^2);
S3 = brinkman_stokes3d_base(p3);
[X3, Y3, Z3] = ndgrid(p3.x0+((1:p3.Nx)-0.5)*p3.h, p3.y0+((1:p3.Ny)-0.5)*p3.h, p3.z0+((1:p3.Nz)-0.5)*p3.h);
hat_eps = 2*p3.h/(3*sqrt(2));
phi3 = 0.5*(1 + tanh((sqrt(X3.^2 + Y3.^2 + Z3.^2) - 0.2)/(2*sqrt(2)*hat_eps)));
band = find(phi3 > 0.2 & phi3 < 0.8); fl = find(phi3 > 0.95);
cells3 = [band(round(linspace(1, numel(band), 3)))', fl(round(linspace(1, numel(fl), 2)))'];
for obj = {'drag'}
    p3.objective = obj{1};
    [J0, sens] = brinkman3d_stokes_solve_J(p3, S3, phi3);
    fprintf('\n3D Stokes, %s (J = %.6e)\n', obj{1}, J0);
    for i = cells3
        e = zeros(size(phi3)); e(i) = delta;
        fd = (brinkman3d_stokes_solve_J(p3, S3, phi3+e) - brinkman3d_stokes_solve_J(p3, S3, phi3-e))/(2*delta);
        rel = abs(fd - sens(i))/max(abs(fd), 1e-14); worst = max(worst, rel);
        fprintf('  cell %5d (phi %.2f): FD %+.8e  adjoint %+.8e  rel %.1e\n', i, phi3(i), fd, sens(i), rel);
    end
end
fprintf('\nlargest relative difference %.1e\n', worst);
assert(worst < 1e-4, 'adjoint gradient check failed: %.2e', worst);
fprintf('GRADIENT CHECKS PASSED\n');
