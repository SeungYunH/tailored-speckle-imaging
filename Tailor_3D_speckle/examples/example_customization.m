% example_customization.m
addpath('../core');
addpath('../utils');

% 1) Build PDFs
Dt = generateDeltaPDF([800 600]);
R  = generateSuperRayleighPDF(1, [800 600]);

pdf_list = [Dt, R, Dt];          % three planes
z_targets = [-1, 0, 1];          % dimensionless

% 2) Build param struct
param = struct();
param.SLMsize     = [800 600];
param.target_size = [80 60];
param.NA          = 0.95;
param.NAportion   = 0.9;
param.wavelength  = 0.488;   % um
param.pad_ratio   = 2;
param.batch_size  = 2;
param.SLMpix      = 20;      % um
param.objMag      = 40;
param.f1          = 500e3;
param.f2          = 62.9e3;
param.NumIter     = 50;
param.z_targets   = z_targets;

% 3) Run customization
[E_update, OptData, param] = CustomizeSpeckles(pdf_list, param);

% 4) Visualize output phase
figure; imagesc(OptData.SLMphase(:,:,1)); axis image off;
title('Optimized SLM Phase (first realization)');