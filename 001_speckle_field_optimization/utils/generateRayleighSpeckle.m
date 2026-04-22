function E = generateRayleighSpeckle(gridSize, NA, wavelength, effPix, varargin)
% generateRayleighSpeckle  Pupil-limited circular-complex-Gaussian (Rayleigh) speckle.
%
% Usage:
%   E = generateRayleighSpeckle(gridSize, NA, wavelength, effPix)
%   E = generateRayleighSpeckle(gridSize, NA, wavelength, effPix, NAportion, batch_size)
%
% Inputs:
%   gridSize   - [Ny, Nx] field size in pixels
%   NA         - objective numerical aperture
%   wavelength - wavelength in um
%   effPix     - effective pixel pitch at the sample plane in um
%   NAportion  - (optional) fraction of the pupil radius used (default 1)
%   batch_size - (optional) number of independent realizations (default 1)
%
% Output:
%   E          - complex speckle field, size [Ny, Nx, batch_size], with
%                mean intensity normalized to 1.
%
% The field is built by passing white circular-complex Gaussian noise
% through a hard NA pupil mask (radius NA*NAportion / lambda in 1/um).
% By the central limit theorem, its intensity follows an exponential
% (Rayleigh) distribution.

    % ----- Optional arguments -----
    NAportion  = 1;
    batch_size = 1;
    if numel(varargin) >= 1 && ~isempty(varargin{1}), NAportion  = varargin{1}; end
    if numel(varargin) >= 2 && ~isempty(varargin{2}), batch_size = varargin{2}; end

    Ny = gridSize(1);
    Nx = gridSize(2);

    % ----- Pupil mask in spatial frequency (ifftshifted to match fft2) -----
    kgrid = ndgrid_matSizeIn([Ny Nx], 1, 'centerZero_ifftshift');   % pix^-1
    kx = kgrid{1} / effPix;     % 1/um
    ky = kgrid{2} / effPix;     % 1/um
    Kradius = NA * NAportion / wavelength;                          % 1/um
    pupil = (kx.^2 + ky.^2) <= Kradius^2;

    % ----- White complex-Gaussian noise, pupil-filtered -----
    E = randn(Ny, Nx, batch_size) + 1i*randn(Ny, Nx, batch_size);
    E = ifft2(fft2(E) .* repmat(pupil, 1, 1, batch_size));

    % ----- Normalize mean intensity to 1 -----
    E = E ./ sqrt(mean(abs(E).^2, 'all'));
end
