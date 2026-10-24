% FIG_SEC53  Pironneau-body geometry and flow panels from run_sec53_body3d.
% Also writes an older two-method history panel; run fig_pc_main afterward
% to replace it with the published preconditioned-FOM/local/compact panel.
% Files go to
% results/figures/ under the names of the paper:
%   sub_body3d_box_init.png, sub_body3d_box_final.png : surfaces phi = 1/2 in D
%   body3d_Jhist40.pdf                                 : drag histories
%   sub_body3d_flow3d_final.png                        : final flow, speed on z = 0
% The quarter-duct fields are reflected to the full duct for display.
%   matlab -batch "setup_paths; fig_sec53"
set(0, 'DefaultAxesFontSize', 12);
F = load(results_dir('sec53_fom.mat')); F = F.res;
R = load(results_dir('sec53_rom.mat')); R = R.res;
table_rows('Fresh-LU comparison: Pironneau''s body', F, R);
fprintf('  full order: %.1f s per design iteration, %.0f%% of it in factorizations\n', F.time/numel(F.J), 100*F.t_fact/F.time);
if isfile(results_dir('sec53_rom_rule.mat'))
    Q = load(results_dir('sec53_rom_rule.mat')); Q = Q.res; n = numel(Q.J);
    fprintf('  Algorithm 2 under its own stopping rule: %d iterations, drag %.7e', n, Q.J(end));
    if n <= numel(F.J), fprintf(', gap to the full-order run at that iteration %.1e', abs(Q.J(end) - F.J(n))/abs(F.J(n))); end
    fprintf('\n');
end
fprintf('  drag %.6e -> %.6e (%.1f%% reduction)\n', F.J(1), F.J(end), 100*(1 - F.J(end)/F.J(1)));

% ---------------- drag histories ----------------
fh = figure('Position', [100 100 520 340], 'Visible', 'off');
plot(0:numel(R.J)-1, R.J, 'r--', 'LineWidth', 1.6); hold on
plot(0:numel(F.J)-1, F.J, 'k-', 'LineWidth', 1.0);
grid on; xlabel('iteration $k$', 'Interpreter', 'latex'); ylabel('drag $\mathcal{D}(\phi_k)$', 'Interpreter', 'latex');
legend({'ROM', 'full order'}, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 13);
exportgraphics(fh, results_dir('figures', 'body3d_Jhist40.pdf'), 'ContentType', 'vector'); close(fh);

% ---------------- shapes: full duct by reflection in y and z ----------------
par = F.par; nIt = numel(F.J);
refl = @(p, sy, sz) permute(cat(3, sz*flip(cat(2, sy*flip(p,2), p),3), cat(2, sy*flip(p,2), p)), [2 1 3]);
xc = F.xc(:).'; yf = [-fliplr(F.yc), F.yc]; zf = [-fliplr(F.zc), F.zc];
[X, Y, Z] = meshgrid(xc, yf, zf);
P0 = refl(F.phi0, 1, 1); Pf = refl(F.phi, 1, 1);
body = [0.85 0.33 0.10];
for s = 1:2
    if s == 1, P = P0; ttl = 'initial design'; fn = 'sub_body3d_box_init.png';
    else,      P = Pf; ttl = sprintf('final design, %d iterations', nIt); fn = 'sub_body3d_box_final.png'; end
    fv = isosurface(X, Y, Z, P, 0.5);
    L = max(fv.vertices(:,1)) - min(fv.vertices(:,1)); Dm = max(fv.vertices(:,2)) - min(fv.vertices(:,2));
    fprintf('  %-14s length %.3f, diameter %.3f, aspect ratio %.2f\n', ttl(1:min(end,14)), L, Dm, L/Dm);
    g = figure('Position', [100 100 700 400], 'Visible', 'off', 'Color', 'w'); axes('Position', [0 0 1 1]);
    q = patch(fv); isonormals(X, Y, Z, P, q); set(q, 'FaceColor', body, 'EdgeColor', 'none'); hold on
    draw_duct(par.x0, par.x0 + par.Nx*par.h);
    daspect([1 1 1]); axis([-0.72 1.32 -0.52 0.52 -0.78 0.72]); axis off; view(-32, 16);
    camlight('headlight'); lighting gouraud; material dull
    text(0.3, 0.5, -0.70, ttl, 'FontSize', 16, 'Interpreter', 'latex', 'HorizontalAlignment', 'center');
    exportgraphics(g, results_dir('figures', fn), 'Resolution', 300); close(g);
end
% half-angles of nose and tail: phi = 1/2 contour on the plane next to z = 0,
% straight-line fit within 0.06 of each tip (Pironneau: 60 degrees)
sec = [flip(F.phi(:,:,1), 2), F.phi(:,:,1)];
C = contourc(xc, yf, sec', [0.5 0.5]); q = 1; Pc = [];
while q < size(C,2), np = C(2,q); Pc = [Pc, C(:,q+1:q+np)]; q = q + np + 1; end %#ok<AGROW>
up = Pc(:, Pc(2,:) > 0.01); [~, i1] = min(up(1,:)); [~, i2] = max(up(1,:));
for t = [i1 i2]
    ii = abs(up(1,:) - up(1,t)) <= 0.06 & up(2,:) > 0.01;
    pp = polyfit(up(1,ii), up(2,ii), 1);
    fprintf('  half-angle at x = %+.3f: %.0f degrees\n', up(1,t), atand(abs(pp(1))));
end

% ---------------- final flow: speed on z = 0, streamlines, body ----------------
S = brinkman_stokes3d_base(par);
[~, ~, o] = brinkman3d_stokes_solve_J(par, S, F.phi);
Nx = par.Nx; Ny = par.Ny; Nz = par.Nz; h = par.h;
U = refl(0.5*(o.u(1:Nx,:,:) + o.u(2:Nx+1,:,:)), 1, 1);
V = refl(0.5*(o.v(:,1:Ny,:) + o.v(:,2:Ny+1,:)), -1, 1);
W = refl(0.5*(o.w(:,:,1:Nz) + o.w(:,:,2:Nz+1)), 1, -1);
SPD = sqrt(U.^2 + V.^2 + W.^2);
fh = figure('Position', [100 100 900 520], 'Visible', 'off', 'Color', 'w'); axes('Position', [0.02 0.02 0.9 0.96]);
hsl = slice(X, Y, Z, SPD, [], [], 0); set(hsl, 'EdgeColor', 'none', 'FaceColor', 'interp'); hold on
colormap(fh, parula); clim([0 0.05]);                 % inlet centreline speed 0.05
fv = isosurface(X, Y, Z, Pf, 0.5); q = patch(fv); isonormals(X, Y, Z, Pf, q);
set(q, 'FaceColor', body, 'EdgeColor', 'none');
[sy, sz] = meshgrid(linspace(-0.44, 0.44, 12), [0.03 0.11 0.21 0.33]);   % seeds just inside the inlet
sx = (par.x0 + 2*h)*ones(size(sy));
verts = stream3(X, Y, Z, U, V, W, sx(:), sy(:), sz(:), [0.25, 30000]);
for i = 1:numel(verts)                               % cut each streamline where it enters the body
    v = verts{i}; if isempty(v), continue, end
    idx = find(interp3(X, Y, Z, Pf, v(:,1), v(:,2), v(:,3)) < 0.5, 1, 'first');
    if ~isempty(idx), verts{i} = v(1:max(idx-1, 1), :); end
end
hs = streamline(verts); set(hs, 'Color', [0.15 0.15 0.15], 'LineWidth', 0.7);
draw_duct(par.x0, par.x0 + Nx*h, false);
daspect([1 1 1]); axis([par.x0 par.x0+Nx*h -0.52 0.52 -0.52 0.52]); axis off; view(-30, 22);
camlight('headlight'); lighting gouraud; material dull
cb = colorbar('Position', [0.9 0.25 0.02 0.5]); cb.Label.Interpreter = 'latex';
cb.Label.String = '$|\mathbf{u}|$ on $z = 0$'; cb.Label.FontSize = 14; cb.FontSize = 11;
text(0.3, 0, -0.68, sprintf('final design, %d iterations', nIt), 'FontSize', 15, 'Interpreter', 'latex', 'HorizontalAlignment', 'center');
exportgraphics(fh, results_dir('figures', 'sub_body3d_flow3d_final.png'), 'Resolution', 300); close(fh);

function draw_duct(x0, x1, inflow)
% edges of the duct (-1/2,1/2)^2, the inlet face shaded and an inflow arrow
if nargin < 3, inflow = true; end
c = [-0.5 0.5]; ed = [0.35 0.35 0.35];
for y = c, for z = c, plot3([x0 x1], [y y], [z z], '-', 'Color', ed, 'LineWidth', 0.8); end, end
for x = [x0 x1]
    plot3([x x x x x], [c(1) c(2) c(2) c(1) c(1)], [c(1) c(1) c(2) c(2) c(1)], '-', 'Color', ed, 'LineWidth', 0.8);
end
if inflow
    patch([x0 x0 x0 x0], [c(1) c(2) c(2) c(1)], [c(1) c(1) c(2) c(2)], [0.55 0.7 0.95], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    quiver3(x0-0.02, 0, 0.62, 0.3, 0, 0, 'k', 'LineWidth', 1.2, 'MaxHeadSize', 0.8, 'AutoScale', 'off');
    text(x0+0.1, 0, 0.74, 'inflow', 'FontSize', 13, 'Interpreter', 'latex', 'HorizontalAlignment', 'center');
end
end
