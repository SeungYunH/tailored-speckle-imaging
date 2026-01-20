% example_propagation.m
addpath('../core');
addpath('../utils');

% Load or generate a test field
Ny = 200; Nx = 300; K = 2;
E0 = randn(Ny,Nx,K) + 1i*randn(Ny,Nx,K);

% Set up propagation
effPix = 0.2;   % µm
wl     = 0.488; % µm
krho2  = computeKrho2(Nx,Ny,effPix);
k_parax= wl * krho2;
zVec   = [-1, 0, 1];  % µm

E_stack = angularSpectrumPropagation_paraxial(E0, zVec, k_parax);

% Display first batch slice
figure;
for iz = 1:numel(zVec)
    subplot(1, numel(zVec), iz);
    imagesc(abs(E_stack(:,:,iz,1)).^2);
    axis image off; title(sprintf('z=%.1fµm', zVec(iz)));
end
sgtitle('Paraxial Propagation Example');