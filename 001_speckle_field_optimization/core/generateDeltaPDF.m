function data = generateDeltaPDF(SLMsize)
% generateDeltaPDF Generates uniform intensity PDF (delta-type speckles)
%
% Input:
%   SLMsize - [Nx, Ny] size of SLM in pixels (default = [199, 199])
%
% Output:
%   data.pdf   - uniform PDF vector
%   data.name  - label string

    if nargin < 1
        SLMsize = [199, 199];
    end

    N = prod(SLMsize);
    data.pdf = ones(N, 1);
    data.name = 'D';
end
