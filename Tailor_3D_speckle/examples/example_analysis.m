% example_analysis.m
addpath('../core');
addpath('../utils');

% Assume E0, OptData, param are available from customization

% Generate target images
TargetImages = GetTargetImage(E0, param);

% Plot axial slices of first realization
figure;
for iz = 1:numel(param.z_targets)
    subplot(1, numel(param.z_targets), iz);
    imagesc(TargetImages(:,:,iz,1));
    axis image off; title(sprintf('z = %.1f', param.z_targets(iz)));
end
sgtitle('Customized Speckle Profiles');

% Compute and plot contrast
contrast = ContrastFigure(TargetImages, param);