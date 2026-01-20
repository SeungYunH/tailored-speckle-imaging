function data = generateSuperRayleighPDF(order, target_size, DispHistogram)
% generateSuperRayleighPDF Generates target PDF for super-Rayleigh speckles
%
% Inputs:
%   order         - exponent applied to the Rayleigh speckle field
%   SLMsize       - [Nx, Ny] size of SLM in pixels (default = [199, 199])
%   DispHistogram - (optional) flag to display histogram (default = 0)
%
% Output:
%   data.pdf      - intensity PDF vector (sorted, thinned)
%   data.name     - name string for output ('R', 'Sup2', etc.)

    if nargin < 2 || isempty(target_size)
        target_size = [199, 199];
    end
    if nargin < 3
        DispHistogram = 0;
    end

    Nx = target_size(1);
    Ny = target_size(2);
    K = 100;  % Number of realizations

    % Generate Rayleigh speckles using FFT
    [I, speckle] = GenRaySpeck(rand(Nx, Ny, K));

    % Super-Rayleigh manipulation
    supR_speckle = speckle.^order ./ sqrt(mean(abs(speckle.^order).^2, 'all'));
    supR_speckle = FFT2(abs(iFFT2(speckle)) .* exp(1i * angle(iFFT2(supR_speckle))));
    supR_I = abs(supR_speckle).^2;

    % Contrast output
    fprintf('Contrasts: SuperR=%.3f, Rayleigh=%.3f\n', std(supR_I,1,"all"), std(I,1,"all"))

    % Build PDF by sorting and thinning
    supR_pdf = reshape(supR_I, [Nx * Ny * K, 1]);
    supR_pdf = sort(supR_pdf);
    supR_pdf = supR_pdf(1:K:end);  % Reduce to 1D PDF vector

    % Optional histogram display
    if DispHistogram
        figure(31); clf; hold on;
        hist_err(I, 0, 2, 'XLim', 10, 'LineWidth', 2.5);
        hist_err(supR_I, 0, 2, 'XLim', 10, 'LineWidth', 2.5);
        set(gca, 'YScale', 'log'); FS; legend('Rayleigh', 'SuperR');
    end

    % Return struct
    data.pdf = supR_pdf;
    if order == 1
        data.name = 'R';
    else
        data.name = ['Sup' num2str(order)];
    end
end
