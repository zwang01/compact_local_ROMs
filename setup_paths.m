function root = setup_paths()
% SETUP_PATHS  Put the folders of code_brinkman_shapeopt on the MATLAB path.
% Results go to results/ in this folder (created on first use), or to the
% folder named by the environment variable BRK_RESULTS (see results_dir.m).
%
%   matlab -batch "setup_paths; run_tests"
root = fileparts(mfilename('fullpath'));
sub = {fullfile('src','model'), fullfile('src','rom'), fullfile('src','design'), ...
       fullfile('src','drivers'), fullfile('src','util'), 'experiments', 'figures', 'tests'};
for i = 1:numel(sub), addpath(fullfile(root, sub{i})); end
end
