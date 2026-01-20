function [I, speckle] = GenRaySpeck(phase)
% GenRaySpeck Generates Rayleigh speckles using FFT propagation
%
% Inputs:
%   phase - [Nx x Ny x K] or [Nx*Ny x K] matrix of random phases
%
% Outputs:
%   I        - intensity patterns [Nx x Ny x K]
%   speckle  - complex speckle fields [Nx x Ny x K]

    if ndims(phase) == 3
        [Nx, Ny, K] = size(phase);
        phaser = exp(1i * 2 * pi * phase);

    elseif ismatrix(phase)
        [N_flat, K] = size(phase);
        % Try to infer rectangular dimensions from a known SLM shape
        error('Matrix input is ambiguous. Please provide phase in [Nx x Ny x K] format for rectangular SLMs.');
        % Alternatively, you could allow passing a shape argument
    else
        error('Invalid phase input dimensions.');
    end

    speckle = zeros(Nx, Ny, K);
    for k = 1:K
        speckle(:,:,k) = fft2(phaser(:,:,k));
    end

    % Normalize power per realization
    speckle = speckle ./ sqrt(mean(abs(speckle).^2, [1 2]));
    I = abs(speckle).^2;
end
