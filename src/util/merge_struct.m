function a = merge_struct(a, varargin)
% MERGE_STRUCT  Copy the fields of each further struct into a (later ones win):
%   cfg = merge_struct(common, struct('tol_eta', 0), extra)
for j = 1:numel(varargin)
    b = varargin{j}; fn = fieldnames(b);
    for i = 1:numel(fn), a.(fn{i}) = b.(fn{i}); end
end
end
