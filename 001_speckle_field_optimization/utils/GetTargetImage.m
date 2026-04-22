function TargetImages = GetTargetImage(E0, param)
% GetTargetImage - Generates axial intensity stack from 3D speckle field
% using paraxial angular spectrum propagation and visualizes the first realization.
%
% Inputs:
%   E0    - complex speckle field at z = 0 [Ny x Nx x K]
%   param - structure with fields:
%              .z_target_physical - physical axial positions (um)
%              .k_parax           - propagation kernel matrix [Ny x Nx]
%              .effPix            - effective pixel size in um
%              .range_x           - cropping indices in x
%              .range_y           - cropping indices in y
%              .batch_size        - number of speckle realizations (K)
%
% Output:
%   TargetImages - 4D intensity stack [dim_target x dim_target x Z x K]

    zVec = param.z_target_physical;

    % Propagate fields
    E_stack = angularSpectrumPropagation(E0, zVec, param.k_parax);

    % Compute intensity and crop target region
    TargetImages = abs(E_stack).^2;
    TargetImages = TargetImages(param.range_x, param.range_y, :, :);  % [HxW x Z x K]

    % Visualization of the first realization with tiled layout that fits
    % the figure window (choose rows/cols based on image aspect ratio).
    Nz = numel(zVec);
    [Hroi, Wroi, ~, ~] = size(TargetImages);
    imgAspect = Wroi / Hroi;

    % Screen-aware figure window
    scr = get(groot, 'ScreenSize');        % [x y w h]
    figW = round(0.9 * scr(3));
    figH = round(0.7 * scr(4));
    figX = round((scr(3) - figW) / 2);
    figY = round((scr(4) - figH) / 2);

    % Pick (nRows, nCols) so the per-tile aspect matches the figure aspect.
    figAspect = figW / figH;
    nCols = max(1, round(sqrt(Nz * figAspect / imgAspect)));
    nCols = min(nCols, Nz);
    nRows = ceil(Nz / nCols);

    figure('Position', [figX figY figW figH]);
    t = tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
    for iz = 1:Nz
        nexttile;
        imagesc(TargetImages(:, :, iz, 1)); axis image off;
        title(sprintf('z = %.2f %sm', zVec(iz), char(181)));
        colormap turbo;
    end
    title(t, 'Axial Intensity Profiles (First Realization)');
end
