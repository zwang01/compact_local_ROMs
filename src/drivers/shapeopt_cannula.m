function res = shapeopt_cannula(cfg)
% SHAPEOPT_CANNULA  Dissipation minimization for the inflow cannula of the
% paper (Section 5.1; Gao and Ma, Test 1): Stokes-Brinkman flow in the box
% D = [0, 2.35]^2, parabolic inflow on the left opening y in [2, 2.35],
% do-nothing outflow on the bottom opening x in [2, 2.35], no slip elsewhere.
% phi = 1 is prescribed on the two openings and phi = 0 on the rest of the box
% boundary; the initial design is the sharp L-shaped duct of width 0.35 along
% the top and right sides, and the fluid area is held at its initial value.
% nu = d*(d/2)^2/Re with d = 0.35, Re = 0.1.
%
%   res = shapeopt_cannula(cfg)
%
% cfg.tol_eta = 0 : full-order optimization (Algorithm 1 of the paper);
% cfg.tol_eta > 0 : residual-screened local ROM (Algorithm 2).
%
% Options (defaults = the settings of the paper): Nx (94: h = 2.35/Nx),
% eps_fac (2), amax_fac (1), beta (1e4), vol_project (true; the runs of
% Section 5.1 used false), max_itr (12000), stop_tol (1e-3), stop_dJ_tol (1e-4),
% tol_eta (1e-5), tol_eta_max (3e-5), tol_eta_min (1e-9), gate ('cal' | 'eucl'),
% cal_warm (3), cal_mem (10), n_warmup (2), tol_J (1e-4), verify ('dwr' |
% 'backstop'), tol_J_dwr (5e-5), n_check (25), n_check_max (400), fallback
% ('precond' for the ROM, 'refactor' for tol_eta = 0), pc_tol (1e-10), pc_maxit
% (20), pc_refresh_its (m_max = 8), w_max (40), r_max (30), pod_tol (1e-10),
% lsq ('qr'; 'normal' as in the runs of Section 5.1), errstudy (false),
% record ([]), save_phi_at ([]), verbose (true).
% See shapeopt_body2d.m for their meaning and for the fields of res.

def = struct('Nx', 94, 'eps_fac', 2, 'amax_fac', 1, 'beta', 1e4, 'vol_project', true, ...
    'max_itr', 12000, 'stop_tol', 1e-3, 'stop_dJ_tol', 1e-4, ...
    'tol_eta', 1e-5, 'tol_eta_max', 3e-5, 'tol_eta_min', 1e-9, 'gate', 'cal', ...
    'cal_warm', 3, 'cal_mem', 10, 'n_warmup', 2, ...
    'tol_J', 1e-4, 'verify', 'dwr', 'tol_J_dwr', 5e-5, 'n_check', 25, 'n_check_max', 400, ...
    'fallback', [], 'pc_tol', 1e-10, 'pc_maxit', 20, 'pc_refresh_its', 8, ...
    'w_max', 40, 'r_max', 30, 'pod_tol', 1e-10, 'lsq', 'qr', ...
    'errstudy', false, 'record', [], 'save_phi_at', [], 'verbose', true);
if nargin < 1, cfg = struct(); end
cfg = set_defaults(cfg, def);
use_rom = cfg.tol_eta > 0;
if isempty(cfg.fallback), if use_rom, cfg.fallback = 'precond'; else, cfg.fallback = 'refactor'; end, end
lu_explicit = strcmp(cfg.fallback, 'precond');       % explicit LU factors when retained as K_ref

% ---------------- problem ----------------
par = struct();
par.Nx = cfg.Nx; par.Ny = par.Nx; par.x0 = 0; par.y0 = 0;
par.h = 2.35/par.Nx; d = 0.35; Re = 0.1; ym = (d/2)^2;
par.nu = d*ym/Re; par.alpha_max = cfg.amax_fac*(2.5e4*par.nu);
par.ubc = @(x,y) (abs(x) < 1e-12).*(y >= 2 & y <= 2.35).*(y-2).*(2.35-y);
par.vbc = @(x,y) 0*x;
xc = par.x0 + ((1:par.Nx)-0.5)*par.h; yc = par.y0 + ((1:par.Ny)-0.5)*par.h;
par.openRight = false(par.Ny,1);
par.openBottom = (xc(:) >= 2 & xc(:) <= 2.35);
par.objective = 'dissipation';
[Xc, Yc] = ndgrid(xc, yc);
phi0 = double((Yc >= 2 & Yc <= 2.35) | (Xc >= 2 & Xc <= 2.35));   % L-shaped duct
freeze = false(par.Nx, par.Ny);                       % no frozen cells: Dirichlet phi on the box boundary
phiFrozen = phi0;
Vcells = sum(phi0(:));                                % fluid area of the initial duct, in cells

% ---------------- design update (Allen-Cahn) ----------------
bcL = double(yc(:) >= 2 & yc(:) <= 2.35); bcR = zeros(par.Ny,1);   % phi = 1 on the openings,
bcB = double(xc(:) >= 2 & xc(:) <= 2.35); bcT = zeros(par.Nx,1);   % phi = 0 elsewhere
hat_eps = cfg.eps_fac*par.h/(3*sqrt(2));
ac = struct('tau', 0.01, 'beta', cfg.beta, 'Kin', 5, 'vol_project', cfg.vol_project);
[Lac, bvec] = dirichlet_laplacian(par.Nx, par.Ny, par.h, bcL, bcR, bcB, bcT);
ac.Adec = decomposition(speye(par.Nx*par.Ny) - ac.tau*hat_eps^2*Lac);
ac.cvec = ac.tau*hat_eps^2*bvec;

% ---------------- flow operator, counters, histories ----------------
S = brinkman_base(par);
Ntot = S.Ntot;
brk_count_fact(0);
pc_opts = struct('tol', cfg.pc_tol, 'maxit', cfg.pc_maxit, 'refresh_its', cfg.pc_refresh_its);
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
    accepted = false; eta_s = NaN; eta_a = NaN; etaJ = NaN;
    if use_rom && k >= cfg.n_warmup && rom.r > 0
        [U, W, J, Jvisc, sens, etaE_s, etaE_a, rs, ra, dJdU_r, etaJ] = ...
            brinkman_rom_step(par, S, rom, dvec, au, av, dalpha);
        [eta_s, eta_a, cal] = indicator_screen(cal, cfg.gate, etaE_s, etaE_a, rs, ra, dJdU_r, Kref, Unorm);
        accepted = max(eta_s, eta_a) <= tol_eta;
        if ~accepted && strcmp(cfg.gate,'cal') && ~isempty(Kref) && numel(cal.rs) >= cal.warm
            cal = indicator_calibrate(cal, etaE_s, etaE_a, rs, ra, dJdU_r, Kref, Unorm, eta_s);
        end
    end

    if ~accepted
        % ---- full-order solve: warm-up, fallback (lines 20-22), or Algorithm 1
        K = S.K0 + spdiags(dvec, 0, Ntot, Ntot);
        pc_ok = false;
        if strcmp(cfg.fallback,'precond') && ~isempty(Kref)
            % GMRES through the retained factorization, warm-started from the
            % rejected reduced pair (or the previous pair); refactor only on failure
            t0 = tic;
            [Up, Wp, Jp, Jvp, ~, pinfo] = brinkman_stokes_solve_pc(par, S, K, au, av, Kref, U, W, pc_opts);
            t_pc = t_pc + toc(t0);
            pc_ok = pinfo.ok;
            if pc_ok
                U = Up; W = Wp; J = Jp; Jvisc = Jvp;
                n_pc = n_pc + 1; pc_its(end+1) = sum(pinfo.its); %#ok<AGROW>
                Kref = pinfo.Fref;
                if pinfo.refreshed, Unorm = norm(U); n_refresh = n_refresh + 1; end
            else
                n_pc_fail = n_pc_fail + 1;
            end
        end
        if ~pc_ok                                     % factor: state and adjoint by one LU
            t0 = tic;
            brk_count_fact;
            Kref = brk_lu(K, lu_explicit);
            U = Kref.solve(S.F); Unorm = norm(U);
            [J, Jvisc, dJdU] = brinkman_J_eval(par, S, au, av, U);
            W = Kref.solvet(dJdU);
            t_fact = t_fact + toc(t0);
        end
        sens = brinkman_sens_gather(par, S, dalpha, U, W);
        if use_rom, rom = brinkman_rom_enrich(rom, S, U, W); end
        n_fom = n_fom + 1; n_since_check = 0; n_since_relax = 0; etaJ_max = 0;
    else
        n_rom = n_rom + 1; n_since_check = n_since_check + 1; n_since_relax = n_since_relax + 1;
        etaJ_max = max(etaJ_max, abs(etaJ)/max(abs(J),1e-30));
    end

    % ---- error study (measurement only): full-order solve at the same design
    if cfg.errstudy && accepted
        K = S.K0 + spdiags(dvec, 0, Ntot, Ntot); dK = decomposition(K);   % not counted: not part of the method
        UF = dK \ S.F; [JF, ~, dJdUF] = brinkman_J_eval(par, S, au, av, UF); WF = dK' \ dJdUF;
        [ePs, ePa] = indicator_exact(rs, ra, dJdU_r, Kref, Unorm);
        es = es_append(es, k, etaE_s, etaE_a, ePs, ePa, eta_s, eta_a, etaJ, ...
                       norm(UF-U)/norm(UF), norm(WF-W)/norm(WF), J, JF);
    end

    % ---- relaxation of tau_eta (lines 8-10)
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
        K = S.K0 + spdiags(dvec, 0, Ntot, Ntot);
        pc_ok = false;
        if strcmp(cfg.fallback,'precond') && ~isempty(Kref)
            t0 = tic;
            [UF, WF, JF, JvF, ~, pinfo] = brinkman_stokes_solve_pc(par, S, K, au, av, Kref, U, W, pc_opts);
            t_pc = t_pc + toc(t0);
            pc_ok = pinfo.ok;
            if pc_ok
                n_pc = n_pc + 1; pc_its(end+1) = sum(pinfo.its); %#ok<AGROW>
                Kref = pinfo.Fref;
                if pinfo.refreshed, Unorm = norm(UF); n_refresh = n_refresh + 1; end
            else
                n_pc_fail = n_pc_fail + 1;
            end
        end
        if ~pc_ok
            t0 = tic;
            brk_count_fact;
            Kref = brk_lu(K, lu_explicit);
            UF = Kref.solve(S.F); Unorm = norm(UF);
            [JF, JvF, dJdUF] = brinkman_J_eval(par, S, au, av, UF);
            WF = Kref.solvet(dJdUF);
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
        rom = brinkman_rom_enrich(rom, S, UF, WF);
        % continue from the verified full-order pair
        sens = brinkman_sens_gather(par, S, dalpha, UF, WF);
        J = JF; Jvisc = JvF; U = UF; W = WF;
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
if cfg.verbose, print_summary('cannula', res); end
end
