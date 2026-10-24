function res = eff_compact_body3d(cfg)
% EFF_COMPACT_BODY3D  Isolated compact-ROM experiment for Pironneau's body.
% Based on SHAPEOPT_BODY3D: Pironneau's body of least drag in Stokes flow (paper,
% Section 5.3): the duct D = (x0, x0+Lx) x (-1/2,1/2)^2 with inflow
% u = 0.05(1-4y^2)(1-4z^2), nu = 1, no-slip walls and a do-nothing outlet,
% computed on the quarter duct y, z >= 0 (symmetry planes y = 0, z = 0). The
% initial design is a sphere of radius 0.2 with the equilibrium interfacial
% profile [eq:acprofile], the solid volume is held at its initial value, the
% drag is minimized, and phi = 1 is frozen within 0.1 of the inlet, the outlet
% and the walls.
%
%   res = eff_compact_body3d(cfg)
%
% cfg.tol_eta = 0 : full-order optimization (Algorithm 1 of the paper);
% cfg.tol_eta > 0 : residual-screened local ROM (Algorithm 2).
%
% Options (defaults = the settings of the paper): N (40: h = 1/N), Lx (2),
% x0 (-0.7), eps_fac (2), amax_fac (1), beta (1e3), vol_project (true; the
% sphere starts from the smooth profile, clipping only removes round-off
% excursions, and both settings give the same run), max_itr (6000), stop_tol
% (1e-3), stop_dJ_tol (1e-4), tol_eta (1e-5), tol_eta_max (3e-5), tol_eta_min
% (1e-9), gate ('cal' | 'eucl'), cal_warm (3), cal_mem (10), n_warmup (2), tol_J
% (1e-4), verify ('dwr' | 'backstop'), tol_J_dwr (5e-5), n_check (25),
% n_check_max (400), fallback ('precond' for the ROM, 'refactor' for
% tol_eta = 0), pc_tol (1e-10), pc_maxit (20), pc_refresh_its (m_max = 8),
% w_max (40), r_max (30), pod_tol (1e-10), lsq ('qr'; 'normal' as in the runs
% of Section 5.3), verbose (true). See shapeopt_body2d.m for their meaning and
% for res (here with grid vectors xc, yc, zc of the quarter duct).

def = struct('eff_rescue',2,'eff_goal_rescue',2, ...
    'N', 40, 'Lx', 2, 'x0', -0.7, 'eps_fac', 2, 'amax_fac', 1, 'beta', 1e3, 'vol_project', true, ...
    'max_itr', 6000, 'stop_tol', 1e-3, 'stop_dJ_tol', 1e-4, ...
    'tol_eta', 1e-5, 'tol_eta_max', 3e-5, 'tol_eta_min', 1e-9, 'gate', 'cal', ...
    'cal_warm', 3, 'cal_mem', 10, 'n_warmup', 2, ...
    'tol_J', 1e-4, 'verify', 'dwr', 'tol_J_dwr', 5e-5, 'n_check', 25, 'n_check_max', 400, ...
    'fallback', [], 'pc_tol', 1e-10, 'pc_maxit', 20, 'pc_refresh_its', 8, ...
    'w_max', 40, 'r_max', 30, 'pod_tol', 1e-10, 'lsq', 'qr', ...
    'errstudy',false,'es_every',1,'verbose', true);
if nargin < 1, cfg = struct(); end
cfg = set_defaults(cfg, def);
use_rom = cfg.tol_eta > 0;
if isempty(cfg.fallback), if use_rom, cfg.fallback = 'precond'; else, cfg.fallback = 'refactor'; end, end
lu_explicit = strcmp(cfg.fallback, 'precond');

% ---------------- problem (quarter duct) ----------------
par = struct();
par.h = 1/cfg.N; par.Nx = round(cfg.Lx*cfg.N); par.Ny = cfg.N/2; par.Nz = cfg.N/2;
par.x0 = cfg.x0; par.y0 = 0; par.z0 = 0; par.sym = true;
par.nu = 1; par.alpha_max = cfg.amax_fac*(2.5e4*par.nu); par.qramp = 0.1;
par.uin = @(y,z) 0.05*(1-4*y.^2).*(1-4*z.^2);
par.objective = 'drag';
xc = par.x0+((1:par.Nx)-0.5)*par.h; yc = ((1:par.Ny)-0.5)*par.h; zc = ((1:par.Nz)-0.5)*par.h;
[Xc, Yc, Zc] = ndgrid(xc, yc, zc);
r0 = 0.2;
hat_eps = cfg.eps_fac*par.h/(3*sqrt(2));
rr = sqrt(Xc.^2 + Yc.^2 + Zc.^2);
phi0 = 0.5*(1 + tanh((rr - r0)/(2*sqrt(2)*hat_eps)));    % sphere with the equilibrium profile
freeze = (Xc < par.x0+0.1) | (Xc > par.x0+cfg.Lx-0.1) | (Yc > 0.5-0.1) | (Zc > 0.5-0.1);
phiFrozen = ones(size(phi0));
Vcells = sum(phi0(:));

% ---------------- design update (Allen-Cahn, Neumann Laplacian) ----------------
ac = struct('tau', 0.01, 'beta', cfg.beta, 'Kin', 5, 'vol_project', cfg.vol_project);
ac.Adec = decomposition(speye(numel(phi0)) - ac.tau*hat_eps^2*neumann_laplacian3(par.Nx,par.Ny,par.Nz,par.h));

% ---------------- flow operator, counters, histories ----------------
S = brinkman_stokes3d_base(par);
Ntot = S.Ntot;
brk_count_fact(0);
pc_opts = struct('tol', cfg.pc_tol, 'maxit', cfg.pc_maxit, 'refresh_its', cfg.pc_refresh_its);
if cfg.verbose
    fprintf('[body3d] quarter-duct grid %dx%dx%d (h = %.4f), n = %d, tol_eta %g, fallback %s\n', ...
        par.Nx, par.Ny, par.Nz, par.h, Ntot, cfg.tol_eta, cfg.fallback);
end
tol_eta = cfg.tol_eta;
rom = struct('SnapU',zeros(Ntot,0),'SnapW',zeros(Ntot,0),'Phi',zeros(Ntot,0),'PhiA',zeros(Ntot,0), ...
    'K0Phi',zeros(Ntot,0),'K0tPhiA',zeros(Ntot,0),'r',0,'rA',0, ...
    'w_max',cfg.w_max,'r_max',cfg.r_max,'pod_tol',cfg.pod_tol,'lsq',cfg.lsq);
cal = struct('rs', [], 'ra', [], 'ratio', [], 'warm', cfg.cal_warm, 'mem', cfg.cal_mem);
Kref = []; Unorm = NaN;
if strcmp(cfg.verify, 'dwr'), n_check_cur = cfg.n_check_max; else, n_check_cur = cfg.n_check; end
n_fom = 0; n_rom = 0; n_verify = 0; n_verify_dwr = 0; n_relax = 0;
n_pc = 0; n_pc_fail = 0; n_refresh = 0; pc_its = []; t_pc = 0; t_fact = 0;
n_since_check = 0; n_since_relax = 0; etaJ_max = 0;
nan1 = nan(cfg.max_itr+1, 1);
H = struct('J',nan1,'Jvisc',nan1,'eta',nan1,'eta_s',nan1,'eta_a',nan1,'etaJ',nan1, ...
           'tol',nan1,'mode',nan1,'r',nan1,'rA',nan1);
volerr_max = 0; bvio_max = 0;
phi = phi0; U = []; W = [];
t_rom=0; t_enrich=0; t_design=0; t_rescue=0; t_goal=0;
n_rescue=0; n_rescued=0; n_goal_attempts=0; n_goal_corrected=0;
es=[];
tStart = tic;

for k = 0:cfg.max_itr
    [dvec, au, av, aw, dalpha] = brinkman3d_coeffs(par, S, phi);

    % ---- reduced state and adjoint, screening (Algorithm 2, lines 4-6)
    accepted = false; eta_s = NaN; eta_a = NaN; etaJ = NaN;
    if use_rom && k >= cfg.n_warmup && rom.r > 0
        tr=tic;
        [U, W, J, Jvisc, sens, etaE_s, etaE_a, rs, ra, dJdU_r, etaJ] = ...
            brinkman3d_rom_step(par, S, rom, dvec, au, av, aw, dalpha);
        t_rom=t_rom+toc(tr);
        [eta_s, eta_a, cal] = indicator_screen(cal, cfg.gate, etaE_s, etaE_a, rs, ra, dJdU_r, Kref, Unorm);
        accepted = max(eta_s, eta_a) <= tol_eta;
        for rescue_try=1:cfg.eff_rescue
            if accepted || isempty(Kref), break; end
            tc=tic;
            [U,W,J,Jvisc,sens,etaE_s,etaE_a,rs,ra,dJdU_r,etaJ] = ...
                eff_stokes_correct(par,S,rom,dvec,au,av,aw,dalpha,U,W,rs,Kref);
            [eta_s,eta_a,cal]=indicator_screen(cal,cfg.gate,etaE_s,etaE_a,rs,ra,dJdU_r,Kref,Unorm);
            accepted=max(eta_s,eta_a)<=tol_eta;
            n_rescue=n_rescue+1; n_rescued=n_rescued+double(accepted);
            if accepted
                te=tic; rom=brinkman_rom_enrich(rom,S,U,W); t_enrich=t_enrich+toc(te);
            end
            t_rescue=t_rescue+toc(tc);
        end
        if ~accepted && strcmp(cfg.gate,'cal') && ~isempty(Kref) && numel(cal.rs) >= cal.warm
            cal = indicator_calibrate(cal, etaE_s, etaE_a, rs, ra, dJdU_r, Kref, Unorm, eta_s);
        end
    end

    if ~accepted
        % ---- full-order solve: warm-up, fallback (lines 20-22), or Algorithm 1
        K = S.K0 + spdiags(dvec, 0, Ntot, Ntot);
        pc_ok = false;
        if strcmp(cfg.fallback,'precond') && ~isempty(Kref)
            t0 = tic;
            [Up, Wp, Jp, Jvp, ~, pinfo] = brinkman3d_stokes_solve_pc(par, S, K, au, av, aw, Kref, U, W, pc_opts);
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
        if ~pc_ok
            t0 = tic;
            brk_count_fact;
            Kref = brk_lu(K, lu_explicit);
            U = Kref.solve(S.F); Unorm = norm(U);
            [J, Jvisc, dJdU] = brinkman3d_J_eval(par, S, au, av, aw, U);
            W = Kref.solvet(dJdU);
            t_fact = t_fact + toc(t0);
        end
        sens = brinkman3d_sens_gather(par, S, dalpha, U, W);
        if use_rom, te=tic; rom = brinkman_rom_enrich(rom, S, U, W); t_enrich=t_enrich+toc(te); end
        n_fom = n_fom + 1; n_since_check = 0; n_since_relax = 0; etaJ_max = 0;
    else
        n_rom = n_rom + 1; n_since_check = n_since_check + 1; n_since_relax = n_since_relax + 1;
        etaJ_max = max(etaJ_max, abs(etaJ)/max(abs(J),1e-30));
    end

    % ---- relaxation of tau_eta (lines 8-10)
    if accepted && n_since_relax >= cfg.n_check
        if etaJ_max < cfg.tol_J/10, tol_eta = min(2*tol_eta, cfg.tol_eta_max); n_relax = n_relax + 1; end
        n_since_relax = 0; etaJ_max = 0;
    end
    % ---- objective verification (lines 11-19)
    dwr_trig = strcmp(cfg.verify,'dwr') && isfinite(etaJ) && abs(etaJ)/max(abs(J),1e-30) > cfg.tol_J_dwr;
    if accepted && dwr_trig && n_since_check<n_check_cur && ~isempty(Kref)
        tg=tic; Uc=U; Wc=W; rsc=rs;
        for goal_try=1:cfg.eff_goal_rescue
            [Uc,Wc,Jc,Jvc,sensc,eSc,eAc,rsc,rac,gc,etaJc] = ...
                eff_stokes_correct(par,S,rom,dvec,au,av,aw,dalpha,Uc,Wc,rsc,Kref);
            [esc,eac,cal]=indicator_screen(cal,cfg.gate,eSc,eAc,rsc,rac,gc,Kref,Unorm);
            n_goal_attempts=n_goal_attempts+1;
            if max(esc,eac)<=tol_eta && abs(etaJc)/max(abs(Jc),1e-30)<=cfg.tol_J_dwr
                U=Uc; W=Wc; J=Jc; Jvisc=Jvc; sens=sensc;
                rs=rsc; ra=rac; dJdU_r=gc; eta_s=esc; eta_a=eac;
                etaE_s=eSc; etaE_a=eAc; etaJ=etaJc; dwr_trig=false;
                te=tic; rom=brinkman_rom_enrich(rom,S,U,W); t_enrich=t_enrich+toc(te);
                n_goal_corrected=n_goal_corrected+1;
                break
            end
        end
        t_goal=t_goal+toc(tg);
    end
    % Independent same-design measurement; exclude audit runs from timings.
    if cfg.errstudy && accepted && mod(n_rom,cfg.es_every)==0
        K = S.K0 + spdiags(dvec,0,Ntot,Ntot); dK=decomposition(K);
        UF=dK\S.F;
        [JF,~,gF]=brinkman3d_J_eval(par,S,au,av,aw,UF);
        WF=dK'\gF;
        [ePs,ePa]=indicator_exact(rs,ra,dJdU_r,Kref,Unorm);
        es=es_append(es,k,etaE_s,etaE_a,ePs,ePa,eta_s,eta_a,etaJ, ...
            norm(UF-U)/norm(UF),norm(WF-W)/norm(WF),J,JF);
    end
    if accepted && (n_since_check >= n_check_cur || dwr_trig)
        if strcmp(cfg.gate,'cal') && ~isempty(Kref)
            cal = indicator_calibrate(cal, etaE_s, etaE_a, rs, ra, dJdU_r, Kref, Unorm, eta_s);
        end
        K = S.K0 + spdiags(dvec, 0, Ntot, Ntot);
        pc_ok = false;
        if strcmp(cfg.fallback,'precond') && ~isempty(Kref)
            t0 = tic;
            [UF, WF, JF, JvF, ~, pinfo] = brinkman3d_stokes_solve_pc(par, S, K, au, av, aw, Kref, U, W, pc_opts);
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
            [JF, JvF, dJdUF] = brinkman3d_J_eval(par, S, au, av, aw, UF);
            WF = Kref.solvet(dJdUF);
            t_fact = t_fact + toc(t0);
        end
        relJ = abs(J - JF)/max(abs(JF), 1e-14);
        if relJ > cfg.tol_J
            tol_eta = max(tol_eta/2, cfg.tol_eta_min);
            if cfg.verbose, fprintf('[verify] k = %d: objective mismatch %.2e, tol_eta -> %.2e\n', k, relJ, tol_eta); end
        elseif relJ < cfg.tol_J/10
            tol_eta = min(2*tol_eta, cfg.tol_eta_max);
        end
        if dwr_trig, n_verify_dwr = n_verify_dwr + 1; end
        te=tic; rom = brinkman_rom_enrich(rom, S, UF, WF); t_enrich=t_enrich+toc(te);
        sens = brinkman3d_sens_gather(par, S, dalpha, UF, WF);
        J = JF; Jvisc = JvF; U = UF; W = WF;
        n_verify = n_verify + 1; n_since_check = 0; n_since_relax = 0; etaJ_max = 0;
    end

    % ---- histories
    H.J(k+1) = J; H.Jvisc(k+1) = Jvisc; H.eta(k+1) = max(eta_s, eta_a);
    H.eta_s(k+1) = eta_s; H.eta_a(k+1) = eta_a; H.etaJ(k+1) = etaJ; H.tol(k+1) = tol_eta;
    H.mode(k+1) = double(~accepted); H.r(k+1) = rom.r; H.rA(k+1) = rom.rA;

    % ---- design update (line 24)
    td=tic; [phinew, acinfo] = ac_update_vol(phi, sens, freeze, phiFrozen, Vcells, ac); t_design=t_design+toc(td);
    volerr_max = max(volerr_max, max(abs(acinfo.vol_err)));
    bvio_max = max(bvio_max, max(max(acinfo.phi_max) - 1, -min(acinfo.phi_min)));
    dphi = max(abs(phinew(:) - phi(:)));
    if cfg.verbose && mod(k,25) == 0
        if accepted, lbl = 'ROM'; else, lbl = 'FOM'; end
        fprintf('k = %5d  J = %.7e  [%s] eta = %8.2e  r = %2d/%2d  |dphi| = %.2e  (%.0f s)\n', ...
            k, J, lbl, max(eta_s,eta_a), rom.r, rom.rA, dphi, toc(tStart));
    end
    phi = phinew;
    if stop_rule(cfg, k, dphi, H.J), break, end
end

res = struct('cfg', cfg, 'par', par, 'iters', k+1, 'k_last', k, ...
    'n_fom', n_fom, 'n_rom', n_rom, 'n_verify', n_verify, 'n_verify_dwr', n_verify_dwr, 'n_relax', n_relax, ...
    'n_pc', n_pc, 'n_pc_fail', n_pc_fail, 'n_refresh', n_refresh, 'pc_its', pc_its, ...
    'nfact', brk_count_fact([]), 'time', toc(tStart), 't_pc', t_pc, 't_fact', t_fact, ...
    'phi0', phi0, 'phi', phi, 'xc', xc, 'yc', yc, 'zc', zc, 'Vcells', Vcells, ...
    'volerr_max', volerr_max, 'bvio_max', bvio_max, 'cal', cal, 'es', es);
res.t_rom=t_rom; res.t_enrich=t_enrich; res.t_design=t_design;
res.t_rescue=t_rescue; res.t_goal=t_goal;
res.n_rescue=n_rescue; res.n_rescued=n_rescued;
res.n_goal_attempts=n_goal_attempts; res.n_goal_corrected=n_goal_corrected;
fn = fieldnames(H);
for i = 1:numel(fn), res.(fn{i}) = H.(fn{i})(1:k+1); end
if cfg.verbose, print_summary('body3d', res); end
end
