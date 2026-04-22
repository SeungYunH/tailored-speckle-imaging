function [E0, TargetImages] = CustomizeSpecklesAnalSave(pdf_target_struct, param)
% CustomizeSpecklesAnalSave - Customizes speckles and saves results
%
% Usage:
%   E0 = CustomizeSpecklesAnalSave(pdf_target_struct, param)
%
% Inputs:
%   pdf_target_struct : array of target PDF structs (with .pdf and .name)
%   param              : struct with customization parameters
%
% Output:
%   E0                 : optimized complex speckle field [dim x dim x K]

    % Default axial positions (um) if caller did not supply them
    if ~isfield(param,'z_target_physical') || isempty(param.z_target_physical)
        M = numel(pdf_target_struct);
        param.z_target_physical = zeros(1, M);
    end

    % Run speckle customization
    [E0, errpdf, param] = CustomizeSpeckles(pdf_target_struct, param);

    %% Target image reconstruction
    TargetImages = GetTargetImage(E0, param);

    %% Contrast figure (dense axial sweep: 0.1 * axialUnit spacing)
    [contrast, z_contrast] = ContrastFigure(E0, param);

    %% Save outputs
    z_target_physical = param.z_target_physical;  % physical axial positions (um)
    NumIter   = param.NumIter;
    NAportion = param.NAportion;

    % Get file name (uses dimensionless z_targets for a compact label)
    if isfield(param,'z_targets') && ~isempty(param.z_targets)
        z_label = param.z_targets;
    else
        z_label = round(z_target_physical, 2);  % fallback
    end
    param.FileName = generateFileName(pdf_target_struct, z_label, NAportion);
    % e.g. param.FileName -> "Dt_R_Dt_M=3_d=-1-0-1"

    % Timestamp
    DateTimeStr = string(datetime('now', 'Format', 'MMddyy_HHmmss'));
    fileNameBase = param.FileName + "_batch=" + num2str(param.batch_size) + "_" + DateTimeStr;

    % Save full results
    save(fullfile(param.saveDir, fileNameBase + ".mat"), ...
         "E0", "contrast", "z_contrast", "z_target_physical", "NumIter", "param", "errpdf");

    % % Save just the first frame & phase
    % E1 = E0(:, :, 1);
    % save(fullfile(param.saveDir, "OneE_" + fileNameBase + ".mat"), ...
    %      "E1", "param");

end
