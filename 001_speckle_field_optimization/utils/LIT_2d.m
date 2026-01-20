% Local Intensity Transformation
function I = LIT_2d (I, target)

[x,y,z] = size(I);
I = reshape(I, [x*y, z]);

% I_sort_position stores intensity and its corresponding coordinate
[NumPixel, NumTrial] = size(I);
I_sort_position = zeros(NumPixel, NumTrial, 2);
I_sort_position(:,:,1) = I;
[~,I_sort_position(:,:,2)] = meshgrid(1:NumTrial,1:NumPixel);

for trial = 1:NumTrial
    I_sort_position2 = sortrows(squeeze(I_sort_position(:,trial,:)));
    target = sort(target);
    for i = 1:NumPixel
        I(I_sort_position2(i,2),trial) = target(i);
    end
end

I = reshape(I, [x, y, z]);

% % I_sort_position stores intensity and its corresponding coordinate
% [NumPixel, NumPixel2, NumTrial] = size(I);
% I_sort_position = zeros(NumPixel, NumPixel2, NumTrial, 2);
% I_sort_position(:,:,:,1) = I;
% [~,I_sort_position(:,:,:,2),~] = meshgrid(1:NumTrial,1:NumPixel,1:NumPixel2);
% 
% for trial = 1:NumTrial
%     I_sort_position2 = sortrows(squeeze(I_sort_position(:,:,trial,:)));
%     target = sort(target);
%     for i = 1:NumPixel
%         for j = 1:NumPixel2
%             I(I_sort_position2(i,j,2),trial) = target(i);
%         end
%     end
% end
end
