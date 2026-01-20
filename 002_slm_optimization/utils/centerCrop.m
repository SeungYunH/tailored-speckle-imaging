function E_cropped = centerCrop(E, cropSize)
% centerCropE0 Center‐crop array E0 to the given size
%   Ein0 = centerCropE0(E0, [h, w]) crops the first two dimensions of E0
%   to height h and width w, preserving any additional dimensions.

[h0, w0, ~] = size(E);
r0 = floor((h0 - cropSize(1))/2) + 1;
c0 = floor((w0 - cropSize(2))/2) + 1;
E_cropped = E(r0:r0+cropSize(1)-1, c0:c0+cropSize(2)-1, :);
end
