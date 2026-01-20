function I = LIT_2d(I, target)
% LIT_2d  Local intensity transformation (vectorized)
% Remaps each trial's pixel intensities to match the target PDF.
%
% Inputs:
%   I      - [x-by-y-by-z] array of intensities
%   target - [NumPixel-by-1] vector with desired PDF (length x*y)
%
% Output:
%   I      - remapped intensity array, same size as input

[x, y, z] = size(I);
NumPixel = x*y;
NumTrial = z;

% Flatten spatial dimensions
I2 = reshape(I, NumPixel, NumTrial);

% Sort intensities along each trial, get indices
[~, idx] = sort(I2, 1, 'ascend');

% Sort target PDF once
t_sorted = sort(target);

% Build linear indices for vectorized assignment
offset = (0:NumTrial-1) * NumPixel;          % 1 x z
linearIdx = idx + offset;                     % NumPixel x z

% Prepare output array
I2_new = zeros(size(I2), 'like', I2);

% Tile sorted target for all trials
Tmat = repmat(t_sorted, 1, NumTrial);

% Assign new intensities in one vectorized step
I2_new(linearIdx(:)) = Tmat(:);

% Reshape back to [x y z]
I = reshape(I2_new, x, y, z);
end
