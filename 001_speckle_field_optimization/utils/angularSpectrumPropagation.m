function E_stack = angularSpectrumPropagation(E0, zVec, k_parax)
% angularSpectrumPropagation
% Propagates a batch of complex fields using the EXACT (non-paraxial)
% angular spectrum method. The transfer function applied per frequency is
%     H = exp(-1i * pi * z * k_parax)
% where k_parax is precomputed so that this equals exp(1i*z*(kz - k0)),
% i.e. exact free-space propagation with the constant piston k0*z removed.
% Paraxial limit: k_parax -> wavelength * krho2 (small angles).
%
% Inputs:
%   E0       - (Ny x Nx x Nbatch), complex input fields at z = 0
%   zVec     - (1 x Nz), vector of axial distances in um
%   k_parax  - (Ny x Nx), exact propagation kernel (1/um):
%                k_parax = (k0 - sqrt(k0^2 - krho^2)) / pi
%
% Output:
%   E_stack  - (Ny x Nx x Nz x Nbatch), propagated fields

[Ny, Nx, Nbatch] = size(E0);
Nz = numel(zVec);

% Fourier transform of input fields
E0_fft = fft2(E0);  % (Ny x Nx x Nbatch)

% Allocate output
E_stack = zeros(Ny, Nx, Nz, Nbatch, 'like', E0);

% Loop over z positions
for iz = 1:Nz
    z = zVec(iz);
    H = exp(-1i * pi * z * k_parax);            % (Ny x Nx)
    H = repmat(H, 1, 1, Nbatch);                % (Ny x Nx x Nbatch)
    E_stack(:,:,iz,:) = ifft2(E0_fft .* H);     % (Ny x Nx x Nbatch)
end

end
