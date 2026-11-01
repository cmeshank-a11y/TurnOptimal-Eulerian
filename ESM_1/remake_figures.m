function remake_figures()
%REMAKE_FIGURES  Redraw all twelve figures from results/results.mat without
%   recomputing (run run_turnaware_study first).
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'src'));
L = load(fullfile(here, 'results', 'results.mat'));
R = L.R;
if isfield(R, 'P'), P = R.P; else
    P = struct('tS', 0, 'tL', 6, 'tR', 18, 'tU', 40, 'degS', 30, 'degU', 150, ...
        'vS', 6/3.6, 'vD', 20/3.6, 'vA', 30/3.6, 'tUnload', 900);
end
make_figures(R, R.W, R.G, R.S, P, R.sample, R.splitBase, R.splitOpt, ...
    fullfile(here, 'figs'), isfield(R, 'bench'));
end
