
%% Add paths

repoRoot = fileparts(which('main.m'));
assert(~isempty(repoRoot), 'Cannot find main.m on MATLAB path. cd to repo root or add it to path.');

addpath(genpath(fullfile(repoRoot,'001_speckle_field_optimization')));
addpath(genpath(fullfile(repoRoot,'002_slm_optimization')));
% addpath(genpath(fullfile(repoRoot,'example')));  % if needed


%% Run 3d speckle field optimization

run_tailor3Dspeckle

% optimized field is stored as 'E0'


%% Run slm optimization

run_slm_optimization

% optimized SLM pattern is stored as 'SLMphase'


