function d = results_dir(varargin)
% RESULTS_DIR  Path of the output folder, optionally joined with further path
% parts: results_dir('sec51_fom.mat'), results_dir('figures', 'svd_flow_nse40.pdf').
% The folder is results/ in the package root (created by setup_paths), or the
% folder named by the environment variable BRK_RESULTS if that is set.
base = getenv('BRK_RESULTS');
if isempty(base), base = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'results'); end
if ~exist(fullfile(base, 'figures'), 'dir'), mkdir(fullfile(base, 'figures')); end
d = fullfile(base, varargin{:});
end
