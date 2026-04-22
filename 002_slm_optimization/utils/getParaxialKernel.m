function k_parax = getParaxialKernel(H, W)
% k_parax = getParaxialKernel(H, W)
%
% Computes the EXACT angular-spectrum propagation kernel (non-paraxial)
% packaged so that the existing propagator
%     H_prop = exp(-1i * pi * z * k_parax)
% equals exp(1i * z * (kz - k0)), i.e. exact free-space propagation with
% the constant piston k0*z dropped.
%
% Relation to the paraxial approximation:
%     paraxial  :  k_parax ≈ wavelength * krho2           (small angles)
%     exact     :  k_parax = (k0 - sqrt(k0^2 - krho^2)) / pi
% where k0 = 2*pi/wavelength and krho^2 = (2*pi)^2 * (fx^2 + fy^2).
% Evanescent frequencies (krho > k0) are clamped to 0 via max(·,0); the
% NA pupil zeros them out anyway.
%
% Inputs:
%   H, W     : image height and width (in pixels)
% Output:
%   k_parax  : exact propagation kernel (unit: µm⁻¹)

% System parameters
wavelength = 488e-3;   % µm
SLMpix = 20;           % µm
f1 = 500e3;            % µm
f2 = 62.9e3;           % µm
objMag = 40;

% Compute effective pixel size at sample plane
addMag = f1 / f2;
effPix = SLMpix / addMag / objMag;  % in µm

% Create spatial frequency grid
fx = ifftshift((-floor(W/2):ceil(W/2)-1) / (W * effPix));  % µm⁻¹
fy = ifftshift((-floor(H/2):ceil(H/2)-1) / (H * effPix));  % µm⁻¹
[FX, FY] = meshgrid(fx, fy);

% Compute squared transverse spatial frequency
krho2 = FX.^2 + FY.^2;  % in µm⁻²

% Exact (non-paraxial) propagation kernel.
k0 = 2*pi / wavelength;                                 % µm⁻¹
kz = sqrt(max(k0^2 - 4*pi^2*krho2, 0));                 % evanescent -> 0
k_parax = (k0 - kz) / pi;                               % µm⁻¹
end
