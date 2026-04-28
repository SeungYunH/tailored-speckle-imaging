%% Setup paths and output directory
thisDir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(thisDir,'core')));     % Main algorithm functions
addpath(genpath(fullfile(thisDir,'utils')));    % Utility functions

%% =====================================================================
%  USER PARAMETERS  (edit this block for your setup)
%  =====================================================================

% ---- SLM / target geometry -------------------------------------------
SLMsize     = [600 800];   % Full SLM size [H, W] in pixels
target_size = [600 800];   % Window size for tailoring speckles [H, W] in pixels

% ---- Optical system --------------------------------------------------
% This code assumes SLM is place on the plane conjugate to the imaging
% plane via 4 lenses (two 2f-2f system), where the first two lenses 
% (f1, f2) have a relay magnification of f1/f2 (potentially providing the 
% 1st-order filtering for 002_slm_optimization) and the last two lenses are 
% tube-objective lense pair. 
wavelength = 488e-3;       % Wavelength in um
NA         = 0.95;         % Objective NA
objMag     = 40;           % Objective magnification
f1         = 500e3;        % Relay lens 1 focal length (um)
f2         = 62.9e3;       % Relay lens 2 focal length (um)  -> addMag = f1/f2
SLMpix     = 20;           % SLM pixel pitch (um)

% ---- Optimization settings -------------------------------------------
NAportion  = 0.8;          % Fraction of pupil radius used (0 < NAportion <= 1); prevents loss of high frequencies in experiment
NumIter    = 200;          % Gerchberg–Saxton iterations, up to 10000 for better results
batch_size = 1;            % Number of independent speckle realizations

% ---- Target PDFs per axial plane -------------------------------------
% Choose any combination of:
%   generateDeltaPDF(target_size)              -> uniform intensity   ('D')
%   generateSuperRayleighPDF(1, target_size)   -> Rayleigh speckle    ('R')
%   generateSuperRayleighPDF(n, target_size)   -> super-Rayleigh, n>1 ('Sup<n>')
%   generateBinaryPDF(target_size, ratio)      -> binary (ratio=1 is 50/50)
% The ORDER of this list must match z_targets below.
D   = generateDeltaPDF(target_size);
Bin = generateBinaryPDF(target_size, 1);
% R = generateSuperRayleighPDF(1, target_size);
% Sup = generateSuperRayleighPDF(3, target_size);

pdf_target_struct = [D, D, Bin, D, D];

% ---- Axial plane locations -------------------------------------------

% If using axial decorrelation length as axial distance unit:
effPix = SLMpix / (f1/f2) / objMag;                                % um
E_Ray  = generateRayleighSpeckle(SLMsize, NA, wavelength, effPix, NAportion, 1); % sample Rayleigh speckle
axialUnit = getAxialDecorrelationFWHM(E_Ray, effPix, wavelength, NA);
z_targets = [-2 -1 0 1 2] * 0.7;
z_target_physical = z_targets * axialUnit;

% % If using axial resolution as axial distance unit:
% axialUnit = wavelength / (1 - sqrt(1 - NA^2));       % um
% z_targets = [-2 -1 0 1 2] * 1.5;
% z_target_physical = z_targets * axialUnit;


% ---- Output ----------------------------------------------------------
SaveFile = true;                                   % Save result .mat
saveDir  = fullfile(fileparts(mfilename('fullpath')), 'data');  % output folder

if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

%% =====================================================================
%  Execution  (normally no edits needed below)
%  =====================================================================

%% Derived quantities
addMag = f1 / f2;                     % Additional magnification from relay
effPix = SLMpix / addMag / objMag;    % Effective pixel size on sample plane (um)

%% Assemble parameter struct
param = struct();
param.SLMsize     = SLMsize;
param.target_size = target_size;
param.NAportion   = NAportion;
param.pad_ratio   = 1.2; % Zero-padding factor for the working grid (>=1); accounts for diffraction leakage outside the FoV
param.NumIter     = NumIter;
param.batch_size  = batch_size;
param.SaveFile    = SaveFile;
param.saveDir     = saveDir;

param.wavelength = wavelength;
param.NA         = NA;
param.f1         = f1;
param.f2         = f2;
param.SLMpix     = SLMpix;
param.objMag     = objMag;
param.effPix     = effPix;

param.axialUnit         = axialUnit;          % axial distance unit (um)
param.z_targets         = z_targets;          % dimensionless z (in units of axialUnit), used for filename
param.z_target_physical = z_target_physical;  % physical axial positions (um), used in optimization

%% Sanity check: target PDF list length vs. axial plane count
if numel(param.z_target_physical) ~= numel(pdf_target_struct)
    warning('Length of z_target_physical (%d) must match the number of target PDFs (%d).', ...
            numel(param.z_target_physical), numel(pdf_target_struct));
end

%% Run speckle customization and save result
[E0, TargetImages] = CustomizeSpecklesAnalSave(pdf_target_struct, param);

%% Check the histogram
% figure(31), clf;
% for j = 1:length(pdf_target_struct)
%     hist_err(squeeze(TargetImages(:,:,j,:)),0,1);
% end
