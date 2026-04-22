function [contrast, z_fine] = ContrastFigure(E0, param)
% ContrastFigure - Computes and plots the axial speckle-contrast profile on
% a fine axial grid (0.05 * axialUnit spacing) covering the target planes
% with a 0.2 * axialUnit margin on each side.
%
% Inputs:
%   E0    - complex speckle field at z = 0 [Ny x Nx x K]
%   param - struct with fields:
%             .axialUnit         - axial distance unit (um)
%             .z_target_physical - physical axial positions (um)
%             .k_parax           - propagation kernel (1/um)
%             .range_x, .range_y - cropping indices for the target FoV
%
% Outputs:
%   contrast - contrast values at each z in z_fine
%   z_fine   - dense axial sampling grid (um)

    % ----- Dense axial sampling grid -----
    axialUnit = param.axialUnit;
    z_min = min(param.z_target_physical) - 0.2 * axialUnit;
    z_max = max(param.z_target_physical) + 0.2 * axialUnit;
    z_fine = (z_min : 0.05*axialUnit : z_max);
    if z_fine(end) < z_max - 1e-12
        z_fine(end+1) = z_max;
    end

    % ----- Propagate E0 to the fine grid, compute intensity in target FoV -----
    E = angularSpectrumPropagation(E0, z_fine, param.k_parax);
    I = abs(E).^2;
    I = I(param.range_x, param.range_y, :, :);

    % ----- Speckle contrast = std / mean over (x, y, realization) per z -----
    contrast = squeeze(std(I, 1, [1 2 4]) ./ mean(I, [1 2 4]));

    % ----- Plot -----
    figure; clf; hold on;
    plot(z_fine, contrast, 'LineWidth', 3, 'Color', default_color('r'));
    % Mark the tailored target planes
    zt = param.z_target_physical(:).';
    if ~isempty(zt)
        ct = interp1(z_fine, contrast, zt, 'linear', 'extrap');
        plot(zt, ct, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    end
    yline(1, 'k--', 'LineWidth', 1.5);
    FS(18);
    xlabel('z (\mum)');
    ylabel('Contrast');
    xlim([z_fine(1), z_fine(end)]);
    box on;
    title('Axial Contrast Profile');
end
