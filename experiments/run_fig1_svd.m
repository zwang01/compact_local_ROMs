% RUN_FIG1_SVD  Data of Figure 1 (Section 4 of the paper): normalized singular
% values of full-order state snapshots along a full-order design trajectory,
% for a GLOBAL snapshot set (every 4th iteration of the first 1200 for the
% cannula, every 3rd of the first 800 for the isolated body: 300 and 267
% snapshots) and for WINDOWS of w_max = 40 consecutive iterations starting at
% k = 100, 500, 1000 (cannula) and k = 100, 400, 700 (isolated body).
% Snapshots are normalized to unit length as in brinkman_rom_enrich.m.
% The trajectories are those of the full-order runs of Sections 5.1 and 5.2
% (same settings, the first 1200 and 800 iterations).
% Results: results/fig1_svd_<case>.mat; then fig1_svd draws the panels.
%
%   matlab -batch "setup_paths; run_fig1_svd"
%   matlab -batch "setup_paths; svdcases = {'nse40'}; run_fig1_svd"
if ~exist('svdcases', 'var'), svdcases = {'cannula', 'nse40', 'nse200'}; end
wlen = 40;
common = struct('tol_eta', 0, 'stop_tol', 0, 'vol_project', false);   % as the runs of Sections 5.1-5.2
for c = 1:numel(svdcases)
    sc = svdcases{c};
    switch sc
        case 'cannula'
            rec = struct('every', 4, 'nglob', 300, 'wstart', [100 500 1000], 'wlen', wlen);
            res = shapeopt_cannula(merge_struct(common, struct('max_itr', 1199, 'record', rec)));
        case {'nse40', 'nse200'}
            rec = struct('every', 3, 'nglob', 267, 'wstart', [100 400 700], 'wlen', wlen);
            res = shapeopt_body2d(merge_struct(common, struct('Re', sscanf(sc, 'nse%d'), 'max_itr', 799, 'record', rec)));
    end
    nrm = @(A) A ./ max(vecnorm(A), 1e-300);
    sg = svd(nrm(res.snap.Sg), 'econ'); sg = sg/sg(1);
    kw = rec.wstart; sw = cell(1, numel(kw));
    for q = 1:numel(kw), s = svd(nrm(res.snap.Sw{q}), 'econ'); sw{q} = s/s(1); end
    save(results_dir(sprintf('fig1_svd_%s.mat', sc)), 'sg', 'sw', 'kw', 'wlen', 'sc');
    need = @(s) find(s < 1e-10, 1, 'first');
    fprintf('%s: modes to reach 1e-10*sigma_1: global %d (%d snapshots); windows', sc, need(sg), size(res.snap.Sg, 2));
    for q = 1:numel(kw), fprintf(' k=%d: %d', kw(q), need(sw{q})); end
    fprintf('\n');
end
