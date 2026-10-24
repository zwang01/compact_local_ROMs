function res = shapeopt_body2d(cfg)
% SHAPEOPT_BODY2D  Drag minimization for the isolated-body benchmark of the
% paper (Section 5.2; Gao and Ma, Test 2): steady Navier-Stokes-Brinkman flow
% in the channel D = (x0, xR) x (-1/2, 1/2), parabolic inflow
% u = umax*(1 - 4y^2), do-nothing outflow, no-slip walls; the initial design
% is a disk of radius r0 about the origin and the solid volume is held at
% Vsolid; Re = 2*r*umax/nu with r = 0.2, so nu = 0.4*umax/Re.
%
%   res = shapeopt_body2d(cfg)
%
% cfg.tol_eta = 0  : full-order optimization (Algorithm 1 of the paper);
% cfg.tol_eta > 0  : residual-screened local ROM (Algorithm 2), initial
%                    screening tolerance tau_eta^0 = cfg.tol_eta.
%
% Options (defaults = the settings of the paper):
%  problem   Re (40), umax (2.5), r0 (0.2), Vsolid (0.04*pi), x0 (-1), xR (3),
%            N (80: h = 1/N), eps_fac (2: eps = eps_fac*h), amax_fac (1:
%            alpha_max = amax_fac*2.5e4*nu), beta (1e3)
%  design    vol_project (true; false = clip after the update, as in the runs
%            of Sections 5.1-5.2, see ac_update_vol.m)
%  stopping  max_itr (12000: iterations k = 0..max_itr), stop_tol (1e-3: stop
%            when ||phi_{k+1} - phi_k||_inf < stop_tol and ...), stop_dJ_tol
%            (1e-4: ... |J_k - J_{k-50}|/|J_k| < stop_dJ_tol; 0 = test off);
%            stop_tol = 0 runs exactly max_itr+1 iterations
%  screening tol_eta (1e-5), tol_eta_max (3e-5), tol_eta_min (1e-9), gate ('cal'
%            calibrated indicators | 'eucl' Euclidean residuals), cal_warm (3),
%            cal_mem (10), n_warmup (2)
%  verify    tol_J (1e-4), verify ('dwr': estimate-triggered with trigger
%            tol_J_dwr = theta_v*tau_J = 5e-5 and backstop n_check_max = 400 |
%            'backstop': only after n_check consecutive accepted iterations),
%            n_check (25: also the relaxation period N_c), n_check_max (400)
%  solvers   fallback ('precond': GMRES preconditioned with the retained
%            factorization | 'refactor': factor at every full-order solve;
%            default 'precond' for the ROM, 'refactor' for tol_eta = 0),
%            pc_tol (1e-10), pc_maxit (20), pc_refresh_its (m_max = 20),
%            adj_reuse (reuse the last Newton factorization for the adjoint;
%            default true for the ROM, false for tol_eta = 0)
%  basis     w_max (40), r_max (30), pod_tol (1e-10), lsq ('qr' | 'normal':
%            least squares of the reduced iteration by QR or by the normal
%            equations, as in the runs of Sections 5.1-5.3; see rom_lsq.m)
%  output    errstudy (false; true: recompute every es_every-th accepted
%            reduced iteration at full order and record true errors and
%            indicators in res.es), es_every (1), record (snapshot recorder
%            for Figure 1, see snapshot_recorder.m), save_phi_at (iterations
%            k whose design phi_k is kept in res.phi_saved), verbose (true)
%
% res: histories over k (J, Jvisc, eta, eta_s, eta_a, etaJ, tol, mode: 1 =
% full-order solve, 0 = accepted reduced, r, rA), counts (iters, n_fom =
% warm-up + fallbacks, n_rom = accepted, n_verify, n_verify_dwr, n_pc,
% n_pc_fail, n_refresh, nfact = factorizations of the saddle system), wall
% times (time, t_pc, t_fact), designs (phi0, phi), grid (xc, yc), par, cfg.

def = struct('Re', 40, 'umax', 2.5, 'r0', 0.2, 'Vsolid', 0.04*pi, 'x0', -1, 'xR', 3, ...
    'N', 80, 'eps_fac', 2, 'amax_fac', 1, 'beta', 1e3, 'vol_project', true, ...
    'max_itr', 12000, 'stop_tol', 1e-3, 'stop_dJ_tol', 1e-4, ...
    'tol_eta', 1e-5, 'tol_eta_max', 3e-5, 'tol_eta_min', 1e-9, 'gate', 'cal', ...
    'cal_warm', 3, 'cal_mem', 10, 'n_warmup', 2, ...
    'tol_J', 1e-4, 'verify', 'dwr', 'tol_J_dwr', 5e-5, 'n_check', 25, 'n_check_max', 400, ...
    'fallback', [], 'pc_tol', 1e-10, 'pc_maxit', 20, 'pc_refresh_its', 20, 'adj_reuse', [], ...
    'w_max', 40, 'r_max', 30, 'pod_tol', 1e-10, 'lsq', 'qr', ...
    'errstudy', false, 'es_every', 1, 'record', [], 'save_phi_at', [], 'verbose', true);
if nargin < 1, cfg = struct(); end
cfg = set_defaults(cfg, def);
use_rom = cfg.tol_eta > 0;
if isempty(cfg.fallback),  if use_rom, cfg.fallback = 'precond'; else, cfg.fallback = 'refactor'; end, end
if isempty(cfg.adj_reuse), cfg.adj_reuse = use_rom; end

% ---------------- problem ----------------
par = struct();
par.h = 1/cfg.N;
par.Nx = round((cfg.xR - cfg.x0)/par.h); par.Ny = round(1/par.h);
par.x0 = cfg.x0; par.y0 = -0.5;
par.nu = 0.4*cfg.umax/cfg.Re;                     % Re = 2 r umax / nu, r = 0.2
par.alpha_max = cfg.amax_fac*(2.5e4*par.nu);
x0 = par.x0; umax = cfg.umax;
par.ubc = @(x,y) (abs(x-x0) < 1e-12).*(umax*(1 - 4*y.^2));
par.vbc = @(x,y) 0*x;
par.openRight = true(par.Ny,1); par.openBottom = false(par.Nx,1);
par.objective = 'drag';
par.adj_reuse_jac = cfg.adj_reuse;
par.lu_explicit = strcmp(cfg.fallback, 'precond');   % explicit LU factors when retained as K_ref
xc = par.x0 + ((1:par.Nx)-0.5)*par.h;
yc = par.y0 + ((1:par.Ny)-0.5)*par.h;
[Xc, Yc] = ndgrid(xc, yc);
phi0 = double(Xc.^2 + Yc.^2 >= cfg.r0^2);            % sharp disk
freeze = (Xc < par.x0+0.1) | (Xc > cfg.xR-0.1) | (Yc < par.y0+0.1) | (Yc > par.y0+1-0.1);
phiFrozen = ones(par.Nx, par.Ny);                    % phi = 1 in a band along the channel boundary
Vcells = ((cfg.xR - cfg.x0)*1 - cfg.Vsolid)/par.h^2; % fluid volume in cells

% ---------------- design update (Allen-Cahn) ----------------
hat_eps = cfg.eps_fac*par.h/(3*sqrt(2));
ac = struct('tau', 0.01, 'beta', cfg.beta, 'Kin', 5, 'vol_project', cfg.vol_project);
[Lac, bvec] = dirichlet_laplacian(par.Nx, par.Ny, par.h, ...
    ones(par.Ny,1), ones(par.Ny,1), ones(par.Nx,1), ones(par.Nx,1));
ac.Adec = decomposition(speye(par.Nx*par.Ny) - ac.tau*hat_eps^2*Lac);
ac.cvec = ac.tau*hat_eps^2*bvec;

% ---------------- flow operator, counters, histories ----------------
S = brinkman_base(par);
Ntot = S.Ntot;
brk_count_fact(0);
pc_opts = struct('tol', cfg.pc_tol, 'maxit', cfg.pc_maxit, 'refresh_its', cfg.pc_refresh_its);

% first state at Re > 100 by continuation in nu (Re/8, Re/4, Re/2)
Uwarm = [];
if cfg.Re > 100
    par_c = par;
    for nuc = par.nu*[8 4 2]
        par_c.nu = nuc; Sc = brinkman_base(par_c);
        [~, ~, outc] = brinkman_nse_solve_J(par_c, Sc, phi0, Uwarm);
        Uwarm = outc.U;
    end
end

tol_eta = cfg.tol_eta;
rom = struct('SnapU',zeros(Ntot,0),'SnapW',zeros(Ntot,0),'Phi',zeros(Ntot,0),'PhiA',zeros(Ntot,0), ...
    'K0Phi',zeros(Ntot,0),'K0tPhiA',zeros(Ntot,0),'r',0,'rA',0, ...
    'w_max',cfg.w_max,'r_max',cfg.r_max,'pod_tol',cfg.pod_tol,'lsq',cfg.lsq);
cal = struct('rs', [], 'ra', [], 'ratio', [], 'warm', cfg.cal_warm, 'mem', cfg.cal_mem);
Kref = []; Unorm = NaN;                               % retained factorization and ||U_ref||
if strcmp(cfg.verify, 'dwr'), n_check_cur = cfg.n_check_max; else, n_check_cur = cfg.n_check; end
n_fom = 0; n_rom = 0; n_verify = 0; n_verify_dwr = 0; n_relax = 0;
n_pc = 0; n_pc_fail = 0; n_refresh = 0; pc_its = []; t_pc = 0; t_fact = 0;
n_since_check = 0; n_since_relax = 0; etaJ_max = 0;
nan1 = nan(cfg.max_itr+1, 1);
H = struct('J',nan1,'Jvisc',nan1,'eta',nan1,'eta_s',nan1,'eta_a',nan1,'etaJ',nan1, ...
           'tol',nan1,'mode',nan1,'r',nan1,'rA',nan1);
volerr_max = 0; bvio_max = 0; es = []; rec = cfg.record;
phi_saved = struct('k', {}, 'phi', {});
phi = phi0; U = []; W = [];
tStart = tic;

for k = 0:cfg.max_itr
    if any(k == cfg.save_phi_at), phi_saved(end+1) = struct('k', k, 'phi', phi); end %#ok<AGROW>
    [dvec, au, av, dalpha] = brinkman_coeffs(par, S, phi);

    % ---- reduced state and adjoint, screening (Algorithm 2, lines 4-6)
    accepted = false; tried = false; eta_s = NaN; eta_a = NaN; etaJ = NaN;
    if use_rom && k >= cfg.n_warmup && rom.r > 0 && ~isempty(Uwarm)
        [U, W, J, Jvisc, sens, etaE_s, etaE_a, ~, rs, ra, dJdU_r, etaJ] = ...
            brinkman_nse_rom_step(par, S, rom, dvec, au, av, dalpha, Uwarm);
        tried = true;
        [eta_s, eta_a, cal] = indicator_screen(cal, cfg.gate, etaE_s, etaE_a, rs, ra, dJdU_r, Kref, Unorm);
        accepted = max(eta_s, eta_a) <= tol_eta;
        if ~accepted && strcmp(cfg.gate,'cal') && ~isempty(Kref) && numel(cal.rs) >= cal.warm
            cal = indicator_calibrate(cal, etaE_s, etaE_a, rs, ra, dJdU_r, Kref, Unorm, eta_s);
        end
    end

    if ~accepted
        % ---- full-order solve: warm-up, fallback (lines 20-22), or Algorithm 1
        pc_ok = false;
        if strcmp(cfg.fallback,'precond') && ~isempty(Kref)
            if tried, U0 = U; W0 = W; else, U0 = Uwarm; W0 = []; end
            t0 = tic;
            [Jp, sensp, outp, pinfo] = brinkman_nse_solve_pc(par, S, phi, U0, W0, Kref, pc_opts);
            t_pc = t_pc + toc(t0);
            pc_ok = pinfo.ok;
            if pc_ok
                J = Jp; sens = sensp; out = outp;
                n_pc = n_pc + 1; pc_its(end+1) = pinfo.its; %#ok<AGROW>
                Kref = outp.dJac;
                if pinfo.refreshed, n_refresh = n_refresh + 1; end
            else
                n_pc_fail = n_pc_fail + 1;
            end
        end
        if ~pc_ok                                     % factor (Newton with LU)
            t0 = tic;
            [J, sens, out] = brinkman_nse_solve_J(par, S, phi, Uwarm);
            t_fact = t_fact + toc(t0);
            Kref = out.dJac;
        end
        U = out.U; W = out.W; Jvisc = out.Jvisc;
        Unorm = norm(U);
        if use_rom, rom = brinkman_rom_enrich(rom, S, U, W); end
        Uwarm = U;
        n_fom = n_fom + 1; n_since_check = 0; n_since_relax = 0; etaJ_max = 0;
    else
        n_rom = n_rom + 1; n_since_check = n_since_check + 1; n_since_relax = n_since_relax + 1;
        etaJ_max = max(etaJ_max, abs(etaJ)/max(abs(J),1e-30));
        Uwarm = U;                                    % next solves start from the accepted reduced state
    end

    % ---- error study (measurement only): full-order solve at the same design
    if cfg.errstudy && accepted && mod(n_rom, cfg.es_every) == 0
        nf = brk_count_fact([]);
        parX = par; parX.adj_reuse_jac = false;       % reference adjoint factored at the converged state
        [JF, ~, outF] = brinkman_nse_solve_J(parX, S, phi, Uwarm);
        brk_count_fact(nf);                           % not part of the method's cost
        [ePs, ePa] = indicator_exact(rs, ra, dJdU_r, Kref, Unorm);
        es = es_append(es, k, etaE_s, etaE_a, ePs, ePa, eta_s, eta_a, etaJ, ...
                       norm(outF.U-U)/norm(outF.U), norm(outF.W-W)/norm(outF.W), J, JF);
    end

    % ---- relaxation of tau_eta (lines 8-10): N_c accepted iterations with
    %      |eta_J|/|J| < theta_r*tau_J throughout
    if accepted && n_since_relax >= cfg.n_check
        if etaJ_max < cfg.tol_J/10, tol_eta = min(2*tol_eta, cfg.tol_eta_max); n_relax = n_relax + 1; end
        n_since_relax = 0; etaJ_max = 0;
    end
    % ---- objective verification (lines 11-19)
    dwr_trig = strcmp(cfg.verify,'dwr') && isfinite(etaJ) && abs(etaJ)/max(abs(J),1e-30) > cfg.tol_J_dwr;
    if accepted && (n_since_check >= n_check_cur || dwr_trig)
        if strcmp(cfg.gate,'cal') && ~isempty(Kref)
            cal = indicator_calibrate(cal, etaE_s, etaE_a, rs, ra, dJdU_r, Kref, Unorm, eta_s);
        end
        pc_ok = false; refreshed = false;
        if strcmp(cfg.fallback,'precond') && ~isempty(Kref)
            t0 = tic;
            [JF, sensF, outF, pinfo] = brinkman_nse_solve_pc(par, S, phi, U, W, Kref, pc_opts);
            t_pc = t_pc + toc(t0);
            pc_ok = pinfo.ok;
            if pc_ok
                n_pc = n_pc + 1; pc_its(end+1) = pinfo.its; refreshed = pinfo.refreshed; %#ok<AGROW>
                if refreshed, n_refresh = n_refresh + 1; end
            else
                n_pc_fail = n_pc_fail + 1;
            end
        end
        if ~pc_ok
            t0 = tic;
            [JF, sensF, outF] = brinkman_nse_solve_J(par, S, phi, Uwarm);
            t_fact = t_fact + toc(t0);
        end
        relJ = abs(J - JF)/max(abs(JF), 1e-14);
        if relJ > cfg.tol_J                           % objective not tracked: tighten
            tol_eta = max(tol_eta/2, cfg.tol_eta_min);
            if cfg.verbose, fprintf('[verify] k = %d: objective mismatch %.2e, tol_eta -> %.2e\n', k, relJ, tol_eta); end
        elseif relJ < cfg.tol_J/10                    % room to spare: relax
            tol_eta = min(2*tol_eta, cfg.tol_eta_max);
        end
        if dwr_trig, n_verify_dwr = n_verify_dwr + 1; end
        Kref = outF.dJac; Unorm = norm(outF.U);
        rom = brinkman_rom_enrich(rom, S, outF.U, outF.W);
        % continue from the verified full-order pair
        J = JF; Jvisc = outF.Jvisc; U = outF.U; W = outF.W; sens = sensF; Uwarm = U;
        n_verify = n_verify + 1; n_since_check = 0; n_since_relax = 0; etaJ_max = 0;
    end

    % ---- histories
    if ~isempty(rec), rec = snapshot_recorder(rec, k, U); end
    H.J(k+1) = J; H.Jvisc(k+1) = Jvisc; H.eta(k+1) = max(eta_s, eta_a);
    H.eta_s(k+1) = eta_s; H.eta_a(k+1) = eta_a; H.etaJ(k+1) = etaJ; H.tol(k+1) = tol_eta;
    H.mode(k+1) = double(~accepted); H.r(k+1) = rom.r; H.rA(k+1) = rom.rA;

    % ---- design update (line 24)
    [phinew, acinfo] = ac_update_vol(phi, sens, freeze, phiFrozen, Vcells, ac);
    volerr_max = max(volerr_max, max(abs(acinfo.vol_err)));
    bvio_max = max(bvio_max, max(max(acinfo.phi_max) - 1, -min(acinfo.phi_min)));
    dphi = max(abs(phinew(:) - phi(:)));
    if cfg.verbose && mod(k,25) == 0
        if accepted, lbl = 'ROM'; else, lbl = 'FOM'; end
        fprintf('k = %5d  J = %.7e  [%s] eta = %8.2e  r = %2d/%2d  |dphi| = %.2e\n', ...
            k, J, lbl, max(eta_s,eta_a), rom.r, rom.rA, dphi);
    end
    phi = phinew;
    if stop_rule(cfg, k, dphi, H.J), break, end
end

res = struct('cfg', cfg, 'par', par, 'iters', k+1, 'k_last', k, ...
    'n_fom', n_fom, 'n_rom', n_rom, 'n_verify', n_verify, 'n_verify_dwr', n_verify_dwr, 'n_relax', n_relax, ...
    'n_pc', n_pc, 'n_pc_fail', n_pc_fail, 'n_refresh', n_refresh, 'pc_its', pc_its, ...
    'nfact', brk_count_fact([]), 'time', toc(tStart), 't_pc', t_pc, 't_fact', t_fact, ...
    'phi0', phi0, 'phi', phi, 'xc', xc, 'yc', yc, 'Vcells', Vcells, ...
    'volerr_max', volerr_max, 'bvio_max', bvio_max, 'cal', cal, 'es', es);
res.phi_saved = phi_saved;
fn = fieldnames(H);
for i = 1:numel(fn), res.(fn{i}) = H.(fn{i})(1:k+1); end
if ~isempty(rec), res.snap = snapshot_recorder(rec); end
if cfg.verbose, print_summary(sprintf('body2d Re = %g', cfg.Re), res); end
end
