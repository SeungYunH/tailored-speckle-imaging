function contrast = ContrastFigure(TargetImages, param)
% ContrastFigure - Computes and plots the axial contrast profile
%
% Inputs:
%   TargetImages - 4D intensity stack [dim_target x dim_target x Z x K]
%   param        - structure with fields:
%                    .z_targets  - dimensionless axial positions
%                    .wavelength - in um
%                    .NA         - numerical aperture
%
% Output:
%   contrast     - contrast values at each axial position

    % Compute physical z values
    axialRes = param.wavelength / (1-sqrt(1-(param.NA)^2));  % axial resolution in um
    zVec = param.z_targets * axialRes;

    % Compute contrast across (x, y, realization) at each z
    contrast = squeeze(std(TargetImages, 1, [1 2 4]) ./ mean(TargetImages, [1 2 4]));

    % Plot contrast profile
    figure; clf; hold on;
    plot(zVec, contrast, 'LineWidth', 3, 'Color', default_color('g'));
    yline(1, 'k--', 'LineWidth', 1.5);
    FS(18); % Set font size if FS() is defined in your codebase
    xlabel('z (µm)');
    ylabel('Contrast');
    
    % After plotting your curve…
    if numel(zVec) > 1
        xlim([min(zVec), max(zVec)]);
    else
        % zVec is a single value – force the [0,2] range
        xlim([0, 2]);
    end

    box on;
    title('Axial Contrast Profile');
end
