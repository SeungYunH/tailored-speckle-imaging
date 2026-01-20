function [E_update, errpdf, param] = CustomizeSpeckles(pdf_target_struct, param)
% CustomizeSpeckles (GPU-enabled)
% Generates speckle patterns with tailored intensity PDF using multi-plane GS optimization
% with optional GPU acceleration.

% Check for GPU availability
try
    useGPU = (gpuDeviceCount > 0);
catch
    useGPU = 0;  % gpuDeviceCount not available or errored
end

if useGPU
    g = gpuDevice(1);
    fprintf('Running on GPU: %s\n', g.Name);
else
    fprintf('GPU not available, running on CPU.\n');
end

% ------------ Parse and validate parameters ------------
if ~isfield(param, 'SLMsize'), error('param.SLMsize must be provided'); end
if ~isfield(param, 'NumIter'), param.NumIter = 1000; end
if ~isfield(param, 'NA'), param.NA = 0.95; end
if ~isfield(param, 'NAportion'), param.NAportion = 0.9; end
if ~isfield(param, 'wavelength'), param.wavelength = 0.488; end % in um
if ~isfield(param, 'pad_ratio'), param.pad_ratio = 2; end
if ~isfield(param, 'batch_size'), param.batch_size = 1; end
if ~isfield(param, 'SLMpix'), param.SLMpix = 20; end
if ~isfield(param, 'objMag'), param.objMag = 40; end
if ~isfield(param, 'f1'), param.f1 = 500e3; end
if ~isfield(param, 'f2'), param.f2 = 62.9e3; end
if ~isfield(param, 'target_mask_on'), param.target_mask_on = 0; end

% Derived parameters
param.effPix = param.SLMpix / (param.f1 / param.f2) / param.objMag;  % um
[SLM_H, SLM_W] = deal(param.SLMsize(1), param.SLMsize(2));
[dim_y, dim_x] = deal(ceil(param.pad_ratio * SLM_H), ceil(param.pad_ratio * SLM_W));

% Collect target PDFs
M = numel(pdf_target_struct);
pdf_target = zeros(numel(pdf_target_struct(1).pdf), M);
for i = 1:M
    pdf_target(:,i) = pdf_target_struct(i).pdf;
end

% Regions of interest
NyT = param.target_size(1);
NxT = param.target_size(2);
cx = floor(dim_y/2) + 1;
cy = floor(dim_x/2) + 1;
range_x = cx + (-floor(NyT/2) : ceil(NyT/2)-1);
range_y = cy + (-floor(NxT/2) : ceil(NxT/2)-1);
param.range_x = range_x;
param.range_y = range_y;

target_mask = false(dim_y, dim_x);
target_mask(range_x, range_y) = true;

% Propagation parameters
krho2 = computeKrho2(dim_x, dim_y, param.effPix); % um^-1
k_parax = krho2 * param.wavelength;
param.k_parax = k_parax;

if ~isfield(param, 'z_targets') || isempty(param.z_targets)
    param.z_targets = zeros(1, M);
end
% axialRes = 2 * param.wavelength / (param.NA)^2;  % in um
axialRes = param.wavelength / (1-sqrt(1-(param.NA)^2));  % in um
z_target_physical = param.z_targets * axialRes;  % um
param.z_target_physical = z_target_physical;

% Pupil filter
k_NA = param.NA / param.wavelength;
Kradius = k_NA * param.NAportion;
kgrid = ndgrid_matSizeIn([dim_y dim_x],1,'centerZero_ifftshift');
kgrid = cat(3, kgrid{1}, kgrid{2}) ./ param.effPix;
pupil = (kgrid(:,:,1)./Kradius).^2 + (kgrid(:,:,2)./Kradius).^2 <= 1;

% Initialize field
rng('shuffle');
E_stack = randn(dim_y, dim_x, param.batch_size) + 1i*randn(dim_y, dim_x, param.batch_size);
E_stack = ifft2(fft2(E_stack) .* repmat(pupil,1,1,param.batch_size));
if param.target_mask_on == 1
    E_stack(~target_mask) = 0;
end
MeanI = 1 * mean(abs(E_stack).^2, 'all');
E_stack = angularSpectrumPropagation(E_stack, z_target_physical, k_parax);

% Move data to GPU if available
if useGPU
    pdf_target = gpuArray(pdf_target);
    pupil      = gpuArray(pupil);
    E_stack    = gpuArray(E_stack);
end

% Prepare error tracking
errpdf = zeros(param.NumIter,1);
temp   = zeros(M,1);

tStart = tic;
for iter = 1:param.NumIter
    % Enforce intensity PDF
    for j = 1:M
        It1 = abs(E_stack(:,:,j,:)).^2;
        It = It1(range_x, range_y, :);
        if useGPU
            mG   = gather(mean(It1,'all'));
            mOut = (gather(sum(It1,'all')) - gather(sum(It,'all'))) / max(numel(It1)-numel(It),1);

            It_cpu = gather(It);
            It_cpu = LIT_2d(It_cpu, pdf_target(:,j));
            % It_cpu = It_cpu / mean(It_cpu, 'all') * MeanI;
            It_cpu = It_cpu / mean(It_cpu, 'all') * mG;
            It = gpuArray(It_cpu);
        else
             mG   = mean(It1,'all');
             mOut = (sum(It1,'all') - sum(It,'all')) / max(numel(It1)-numel(It),1);

            It = LIT_2d(It, pdf_target(:,j));
            % It = It / mean(It, 'all') * MeanI;
            It = It / mean(It, 'all') * mG;
        end
        % I_full = abs(E_stack(:,:,j,:)).^2;
        I_full = abs(E_stack(:,:,j,:)).^2 * (mG / max(mOut, eps));
        I_full(range_x, range_y, :) = It;
        E_stack(:,:,j,:) = sqrt(I_full) .* exp(1i * angle(E_stack(:,:,j,:)));
    end

    % Backward propagation and average
    E_update = zeros(dim_y, dim_x, M, param.batch_size, 'like', E_stack);
    for j = 1:M
        E_b = angularSpectrumPropagation(squeeze(E_stack(:,:,j,:)), -z_target_physical(j), k_parax);
        E_update(:,:,j,:) = E_b;
    end
    E_update = squeeze(mean(E_update,3));
    if param.target_mask_on == 1
        E_update(~target_mask) = 0;
    end

    % Forward apply pupil and propagate
    E_update = ifft2(fft2(E_update) .* pupil);
    E_stack  = angularSpectrumPropagation(E_update, z_target_physical, k_parax);

    % Compute PDF error (on CPU)
    for j = 1:M
        Ierr = gather(abs(E_stack(range_x, range_y, j, :)).^2);
        pdf_j = gather(pdf_target(:,j));
        temp(j) = hist_err(squeeze(Ierr), pdf_j);
    end
    errpdf(iter) = sum(temp);

    % Progress display
    elapsed = toc(tStart);
    eta = elapsed / iter * (param.NumIter - iter);
    fprintf('Iter %3d/%d: err = %.4f (%.1fs elapsed, ETA %.2fmin)\n', iter, param.NumIter, errpdf(iter), elapsed, eta/60);
    
    % fprintf('1 Mean intensity %.2f, current mean:%.2f, MeanI %.2f)\n', mean(abs(E_stack).^2,'all'), sqrt(mean(abs(E_stack).^2,'all')), sqrt(MeanI));
    E_stack = E_stack / sqrt(mean(abs(E_stack).^2,'all')) * sqrt(MeanI);
    % fprintf('2 Mean intensity %.2f, current mean:%.2f, MeanI %.2f)\n', mean(abs(E_update).^2,'all'), sqrt(mean(abs(E_update).^2,'all')), sqrt(MeanI));

end

% Gather results back to CPU
if useGPU
    E_update = gather(E_update);
end
end
