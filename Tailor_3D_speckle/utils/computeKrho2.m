function krho2 = computeKrho2(Nx, Ny, effPix)
% computeKrho2
% Computes the squared transverse spatial frequency grid k_rho^2
%
% Inputs:
%   Nx      - number of pixels in x
%   Ny      - number of pixels in y
%   effPix  - effective pixel size [um]
%
% Output:
%   krho2   - (Ny x Nx), squared radial spatial frequency in (1/um)^2

% Define frequency grids (1/um)
[kx, ky] = meshgrid( ...
    ifftshift((-floor(Nx/2):ceil(Nx/2)-1) / (Nx * effPix)), ...
    ifftshift((-floor(Ny/2):ceil(Ny/2)-1) / (Ny * effPix)) ...
);

% Compute squared radial frequency
krho2 = kx.^2 + ky.^2;

end
