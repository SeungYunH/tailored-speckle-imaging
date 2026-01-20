
%%% Add path
thisDir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(thisDir,'utils')));    % Utility functions

%% Parameters

%%% wavelength
wl = 0.488;

%%% Nikon 0.95/40
objNA = 0.95;
objMag = 40;

%%% additional magnification
f1 = 500e3;  % um
f2 = 62.9e3; % um
addMag = f1/f2;

%%% def. ramp
rampVec = [ 0, -1.5*0.95/0.488*objMag/40]; % um^-1

%%% SLM info
SLMpix = 20; % um;
SLMsize = [600, 800];

%%% calibration result
SLM2pival = 145;

%%% effective pixel size
effPix   = SLMpix/addMag;   


%% Pre-defined pattern optimization

%%% Simple patterns
% patTag = 'normal';
% patTag = 'point';
% patTag = 'quarter';
% patTag = 'edge_32';        % Size of edge and central rect.

%%% Speckles
% patTag = 'RayleighSpeckle_0.8';
% patTag = 'BesselSpeckle_0.95';
% patTag = 'FourierRing_0.9';

%%% SIM patterns (2DSIM_%.2f_%03d_%03d)
% patTag = '2DSIM_0.50_000_000';  % sin(theta)/NA, modulation angle, addPhase.
% patTag = '3DSIM_0.90_090_000';  % sin(theta)/NA, modulation angle, addPhase.

%%% Optimization
% [SLMphase, Eout] = EtoSLMpattern(patTag, SLMsize, SLM2pival, effPix,...
%     objNA, objMag, wl,'reg_param',0.2,'displayTF',0, 'rampVec', rampVec,'forceOpt',1);


%% Arbitrary field optimization

%%% Field input
patTag = 'Optimized';

if ~exist('E0','var') || isempty(E0)
    error('Tailor3D:MissingE0', ...
        'E0 is required but was not found in the workspace. Load/define E0 before running this section.');
end

% Crop field size
E0_cropped = centerCrop(E0, SLMsize);

%%% Optimization
[SLMphase, Eout] = EtoSLMpattern(patTag, SLMsize, SLM2pival, effPix,...
    objNA, objMag, wl,'E0', E0_cropped, 'reg_param',0.2,'displayTF',0, 'rampVec', rampVec,'forceOpt',1);








