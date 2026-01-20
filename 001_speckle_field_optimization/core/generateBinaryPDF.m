function data = generateBinaryPDF(SLMsize, ratio)
% generateBinaryPDF Generates binary intensity PDF
%
% Inputs:
%   SLMsize - [Nx, Ny] size of SLM in pixels (default = [199, 199])
%   ratio   - ratio between 0 and 1 for the binary thresholding (default = 1)
%
% Output:
%   data.pdf   - binary PDF vector
%   data.name  - label string

    if nargin < 1
        SLMsize = [199, 199];
    end
    if nargin < 2
        ratio = 1;
    end

    N = prod(SLMsize);
    data.pdf = ones(N, 1);
    data.pdf(1:round(ratio/(ratio+1) * N)) = 0;

    if ratio == 1
        data.name = 'Bin';
    else
        data.name = ['Bin' num2str(ratio)];
    end
end
