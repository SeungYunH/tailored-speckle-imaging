% example_speckle_generation.m
addpath('../utils');

SLMsize = [800 600];
SLMpix   = 20;        % μm
f1       = 500e3;     % μm
f2       = 62.9e3;    % μm
objNA    = 0.95;
wl       = 0.488;     % μm

addMag = f1/f2;
effPix = SLMpix/addMag;  % effective pixel size in μm

E = GenerateRayleighSpeckle(SLMsize, SLMpix, addMag, objNA, wl, 1.0);

figure;
imagesc(abs(E).^2);
axis image off; colormap hot;
title('Example Rayleigh Speckle');