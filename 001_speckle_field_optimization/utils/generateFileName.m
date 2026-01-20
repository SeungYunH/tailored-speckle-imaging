function fileName = generateFileName(pdf_target_struct, z_targets, NAportion)
% generateFileName  Builds a descriptive filename from target PDFs and axial positions
%
% Usage:
%   fileName = generateFileName(pdf_target_struct, z_targets)
%
% Inputs:
%   pdf_target_struct  - 1×M struct array with fields:
%                          .name (string or char)  
%                          .pdf  (vector, not used here but implies length)
%   z_targets          - 1×M numeric vector of dimensionless axial positions
%
% Output:
%   fileName           - character vector, e.g. 'Dt_R_Dt_M=3_d=0-1-2'

    M = numel(pdf_target_struct);
    if numel(z_targets) ~= M
        warning('Length of z_targets (%d) does not match number of PDFs (%d).', ...
                numel(z_targets), M);
    end

    % Combine PDF names
    names = {pdf_target_struct.name};
    combinedName = strjoin(names, '_');

    % before your if: build an "0p9"‐style string
    NAstr = strrep(sprintf('%.2f',NAportion),'.','p');

    % Build filename
    if M > 1
        distStr = strjoin(string(z_targets), '-');  % e.g. "0-1-2"
        fileName = sprintf('%s_M=%d_z=%s_%s', ...
            combinedName, M, distStr, NAstr);
    else
        fileName = sprintf('%s_M=%d_%s', ...
            combinedName, M, NAstr);
    end

end
