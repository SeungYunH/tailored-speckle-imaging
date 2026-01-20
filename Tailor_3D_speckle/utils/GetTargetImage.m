function TargetImages = GetTargetImage(E0, param)
% GetTargetImage - Generates axial intensity stack from 3D speckle field
% using paraxial angular spectrum propagation and visualizes the first realization.
%
% Inputs:
%   E0    - complex speckle field at z = 0 [Ny x Nx x K]
%   param - structure with fields:
%              .z_targets   - list of dimensionless axial positions
%              .k_parax     - wavelength * krho2 matrix [Ny x Nx]
%              .effPix      - effective pixel size in um
%              .wavelength  - wavelength in um
%              .NA          - numerical aperture
%              .range_x     - cropping indices in x
%              .range_y     - cropping indices in y
%              .dim_target  - crop size for output
%              .batch_size  - number of speckle realizations (K)
%
% Output:
%   TargetImages - 4D intensity stack [dim_target x dim_target x Z x K]

    % axialRes = 2 * param.wavelength / param.NA^2;  % axial resolution in um
    axialRes = param.wavelength / (1-sqrt(1-(param.NA)^2));  % in um
    zVec = param.z_targets * axialRes;             % convert to physical units (um)

    % Propagate fields
    E_stack = angularSpectrumPropagation(E0, zVec, param.k_parax);

    % Compute intensity and crop target region
    TargetImages = abs(E_stack).^2;
    TargetImages = TargetImages(param.range_x, param.range_y, :, :);  % [HxW x Z x K]

    % Visualization of the first realization
    Nz = numel(param.z_targets);
    figure;
    for iz = 1:Nz
        subplot(1, Nz, iz);
        imagesc(TargetImages(:, :, iz, 1)); axis image off;
        title(sprintf('z = %.2f µm', zVec(iz)));
        colormap hot;
    end
    sgtitle('Axial Intensity Profiles (First Realization)');
end
