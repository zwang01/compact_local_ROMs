function rec = snapshot_recorder(rec, k, U)
% SNAPSHOT_RECORDER  Collect full-order states along a design trajectory for
% the singular-value study of Figure 1 (run_fig1_svd.m). The drivers call it
% once per design iteration when cfg.record is set.
%
%   rec = snapshot_recorder(rec, k, U)   offer the state U of iteration k
%   rec = snapshot_recorder(rec)         finish: trim the global set
%
% rec.every  : keep every rec.every-th iteration in the GLOBAL set (rec.Sg, rec.kg)
% rec.nglob  : at most this many global snapshots
% rec.wstart : first iterations of the WINDOWS rec.Sw{q}
% rec.wlen   : window length (consecutive iterations)
if nargin == 1
    if isfield(rec, 'ig'), rec.Sg = rec.Sg(:, 1:rec.ig-1); rec.kg = rec.kg(1:rec.ig-1); end
    return
end
if ~isfield(rec, 'ig')                           % first call: allocate
    rec.Sg = zeros(numel(U), rec.nglob); rec.kg = zeros(1, rec.nglob); rec.ig = 1;
    rec.Sw = cell(1, numel(rec.wstart));
    for q = 1:numel(rec.wstart), rec.Sw{q} = zeros(numel(U), rec.wlen); end
end
if mod(k, rec.every) == 0 && rec.ig <= rec.nglob
    rec.Sg(:, rec.ig) = U; rec.kg(rec.ig) = k; rec.ig = rec.ig + 1;
end
for q = 1:numel(rec.wstart)
    j = k - rec.wstart(q) + 1;
    if j >= 1 && j <= rec.wlen, rec.Sw{q}(:, j) = U; end
end
end
