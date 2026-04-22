function [FWHM, z, corr_axial] = getAxialDecorrelationFWHM(E0, effPix, wavelength, NA)
% getAxialDecorrelationFWHM  Full-width-at-half-maximum of the axial
% intensity auto-correlation of a pupil-limited speckle field.
%
% Usage:
%   FWHM = getAxialDecorrelationFWHM(E0, effPix, wavelength, NA)
%
% Inputs:
%   E0         - complex field at z = 0, [Ny x Nx] or [Ny x Nx x B]
%   effPix     - effective pixel pitch at the sample plane (um)
%   wavelength - wavelength (um)
%   NA         - objective numerical aperture
%
% Outputs:
%   FWHM       - axial decorrelation FWHM in um
%   z          - internal z sampling vector used (um)
%   corr_axial - normalized axial correlation C(z)/C(0)
%
% The z test range is sampled internally from -3 to +3 axial resolutions
% in steps of 0.1 axial resolutions, where
%   axialRes = wavelength / (1 - sqrt(1 - NA^2)).

    % ----- Internal z sampling grid (um) -----
    axialRes = wavelength / (1 - sqrt(1 - NA^2));
    z_test   = (-3:0.1:3) * axialRes;

    % ----- Exact (non-paraxial) propagation kernel -----
    [H0, W0, B] = size(E0);
    krho2 = computeKrho2(W0, H0, effPix);
    k0    = 2*pi / wavelength;
    kz    = sqrt(max(k0^2 - 4*pi^2*krho2, 0));
    k_parax = (k0 - kz) / pi;

    Nz  = numel(z_test);
    iz0 = find(z_test == 0, 1);
    if isempty(iz0)
        [~, iz0] = min(abs(z_test));   % safest fallback if 0 not exactly hit
    end

    % Propagate and compute intensity stack
    E = angularSpectrumPropagation(E0, z_test, k_parax);
    I = abs(E).^2;
    I = centerCrop(I, [300 400]);

    % Reference intensity at z = 0
    I0 = I(:,:,iz0,:);
    I0 = I0 - mean(mean(I0,1),2);
    denom0 = sum(sum(I0.^2,1),2);   % 1 x 1 x 1 x B

    corr_all = zeros(B, Nz);
    for iz = 1:Nz
        Iz = I(:,:,iz,:);
        Iz = Iz - mean(mean(Iz,1),2);
        denomz = sum(sum(Iz.^2,1),2);
        numer  = sum(sum(I0 .* Iz,1),2);
        corr_all(:,iz) = reshape(numer ./ sqrt(denom0 .* denomz), [B 1]);
    end

    corr_axial = mean(corr_all, 1, 'omitnan');
    corr_axial = corr_axial / max(corr_axial);

    z = z_test(:).';

    % Linear interpolation for the 0.5 crossing on the z >= 0 branch
    pos   = z >= 0;
    z_pos = z(pos);
    c_pos = corr_axial(pos);
    idx   = find(c_pos <= 0.5, 1, 'first');

    if isempty(idx) || idx == 1
        FWHM = NaN;
        warning('Half-maximum crossing not found within sampled z range.');
    else
        z1 = z_pos(idx-1); z2 = z_pos(idx);
        c1 = c_pos(idx-1); c2 = c_pos(idx);
        z_half = z1 + (0.5 - c1) * (z2 - z1) / (c2 - c1);
        FWHM = 2 * z_half;
    end
end
