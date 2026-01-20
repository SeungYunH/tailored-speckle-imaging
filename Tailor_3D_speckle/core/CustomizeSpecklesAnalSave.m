function [E0, TargetImages] = CustomizeSpecklesAnalSave(pdf_target_struct, param)
% CustomizeSpecklesAnalSave - Customizes speckles and saves results
%
% Usage:
%   E0 = CustomizeSpecklesAnalSave(pdf_target_struct, param)
%   E0 = CustomizeSpecklesAnalSave(pdf_target_struct, param, d_list)
%
% Inputs:
%   pdf_target_struct : array of target PDF structs (with .pdf and .name)
%   param              : struct with customization parameters
%   d_list             : (optional) list of axial distances (in units of axial resolution)
%
% Output:
%   E0                 : optimized complex speckle field [dim x dim x K]

    % Set d_targets if not already set
    if isempty(param.z_targets)
        M = numel(pdf_target_struct);
        param.z_targets = linspace(-2, 2, M-1);
    end

    % Run speckle customization
    [E0, errpdf, param] = CustomizeSpeckles(pdf_target_struct, param);

    %% Target image reconstruction
    TargetImages = GetTargetImage(E0, param);

    %% Contrast figure
    contrast = ContrastFigure(TargetImages, param);

    %% Save outputs
    z_targets = param.z_targets;        % match param field
    NumIter   = param.NumIter;
    NAportion = param.NAportion;

    % Get file name
    param.FileName = generateFileName(pdf_target_struct, z_targets, NAportion);
    % e.g. param.FileName -> "Dt_R_Dt_M=3_d=-1-0-1"

    % Timestamp
    DateTimeStr = string(datetime('now', 'Format', 'MMddyy_HHmmss'));
    fileNameBase = param.FileName + "_batch=" + num2str(param.batch_size) + "_" + DateTimeStr;

    % Save full results
    save(fullfile(param.saveDir, fileNameBase + ".mat"), ...
         "E0", "contrast", "z_targets", "NumIter", "param", "errpdf");

    % Save just the first frame & phase
    E1 = E0(:, :, 1);
    save(fullfile(param.saveDir, "OneE_" + fileNameBase + ".mat"), ...
         "E1", "param");

end
