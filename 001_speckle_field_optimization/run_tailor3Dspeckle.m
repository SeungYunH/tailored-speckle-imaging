%% Setup

thisDir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(thisDir,'core')));     % Main algorithm functions
addpath(genpath(fullfile(thisDir,'utils')));    % Utility functions

% Output directory (inside the same subfolder)
saveDir = fullfile(thisDir, 'data');
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

% SLM settings
SLMsize = [600 800];              % Non-square SLM (HxW)
target_size = [300 400];         


%% Generate target PDFs
% SuperR_order = 3;                                             % Super-Rayleigh exponent
% Sup = generateSuperRayleighPDF(SuperR_order, target_size);    % Super-Rayleigh speckles
R   = generateSuperRayleighPDF(1, target_size);                 % Rayleigh speckles
U  = generateDeltaPDF(target_size);                             % Uniform intensity
Bin = generateBinaryPDF(target_size, 1);                        % Binary


%% Set parameters
param = struct();
param.SLMsize = SLMsize;            % [H, W] in pixels
param.target_size = target_size;
param.NAportion = 0.8;              % Define NA portion for the parameters
param.pad_ratio = 1.2;              % pad ratio factor for working grid
param.NumIter = 500;                % Iterations for GS loop
param.batch_size = 1;               % Number of speckles in batch
param.SaveFile = true;              % Save result
param.saveDir = saveDir;            % Save folder

% additional magnification
f1 = 500e3;  % um
f2 = 62.9e3; % um
addMag = f1/f2;

% SLM info
SLMpix = 20; % um;

% Effective pixel size on the sample plane
objMag = 40;
effPix = SLMpix / addMag / objMag;

% Imaging parameters
param.wavelength = 488e-3;        % Wavelength in um
param.NA = 0.95;                  % Objective NA
param.f1 = f1;
param.f2 = f2;
param.SLMpix = SLMpix;
param.objMag = objMag;
param.effPix = effPix;            % SLM pixel size in sample plane, in um


%% Assemble target list and locations

pdf_target_struct = [U, U, Bin, U, U];

% Axial plane locations (in units of axial resolution)
param.z_targets = [-2 -1 0 1 2] * 1.5;

% Check d_list length vs. pdf_target_struct
if isfield(param, 'd_list') && ~isempty(param.z_targets)
    if numel(param.z_targets) ~= numel(pdf_target_struct)
        warning('Length of param.d_list (%d) must match the number of target PDFs (%d).', ...
                numel(param.z_targets), numel(pdf_target_struct));
    end
end


%% Run speckle customization and save result
[E0, TargetImages] = CustomizeSpecklesAnalSave(pdf_target_struct, param);


%% Check the histogram
% figure(31), clf;
% for j = 1:length(pdf_target_struct)
%     hist_err(squeeze(TargetImages(:,:,j,:)),0,1);
% end



