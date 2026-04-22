function [SLMphase, E_out, E_input] = EtoSLMpattern(patTag, SLMsize, SLM2pival, effPix, objNA, objMag, wl, varargin)
%%% patTag: char
%%% effpix: um
%%% k_NA:   um^-1

%% parsing
p = inputParser();

addParameter(p, 'gpuInd',  1, @isscalar);
addParameter(p, 'monitorInd',  1, @isscalar);
addParameter(p, 'relDataPath', 'SLMPatternStock'); %

addParameter(p, 'padRatio', 2, @isscalar);
addParameter(p, 'NAsafetyFactor', 1.4, @isscalar);
addParameter(p, 'invalidPix', 0, @isscalar); % SLM invalid pix (manual = 8 pix)
addParameter(p, 'rampVec', [0,0], @isvector); % unit: um^-1
addParameter(p, 'E0', 0, @ismatrix);

addParameter(p, 'displayTF', 0, @isscalar); % if 0, return outImg only.
addParameter(p, 'forceOpt', 0, @isscalar); % if 0, return outImg only.
addParameter(p, 'reg_param', 0.01, @isscalar);
addParameter(p, 'SIM_ramp_shift', 0, @isscalar);
addParameter(p, 'AxialPropZ', 0, @isscalar);

%%%
parse(p, varargin{:});
gpuInd   = p.Results.gpuInd;
monitorInd = p.Results.monitorInd;
relDataPath = p.Results.relDataPath;

padRatio = p.Results.padRatio;
NAsafetyFactor = p.Results.NAsafetyFactor;
invalidPix = p.Results.invalidPix;
rampVec = p.Results.rampVec;
E0 = p.Results.E0;
displayTF = p.Results.displayTF;
forceOpt = p.Results.forceOpt;
reg_param_input = p.Results.reg_param;
SIM_ramp_shift = p.Results.SIM_ramp_shift;
AxialPropZ = p.Results.AxialPropZ;

%% def. topPath
topPath = fileparts(mfilename("fullpath"));
topPath = fullfile(topPath, relDataPath);
if ~exist(topPath,"dir")
    mkdir(topPath);
end

%% SAVED file search
saveStr = fullfile( topPath, sprintf('SLMpat_%.2f_%d_%d_%s.mat',objNA,objMag,round(wl*1e3),patTag) );
if exist(saveStr,"file") && (~forceOpt) %%% if the file exists, just load and display
    
    if nargout <= 1
        load(saveStr,'SLMphase');
        
    elseif nargout == 2
        load(saveStr,'SLMphase','E_out');

    else
        error(' Invalid nargout error ')
    end

else %%% if the file doesn't exists, make one and save

    %% gpuTest
    try
        gpuDevice(gpuInd);
    catch
        gpuInd = false;
    end


    %% contants
    k_NA = objNA/wl;
    SLMsize(2) = SLMsize(2) - invalidPix;
    padSize  = padRatio * SLMsize;
    effPix = effPix/objMag;

    %%% SLMmask
    SLMmask = mpad( true(SLMsize), padSize );

    %%% xgrid
    xgrid = ndgrid_matSizeIn(padSize,0,'centerZero'); % pix
    xgrid = cat(3, xgrid{1}, xgrid{2});
    xgrid = xgrid .* effPix; % um

    %%% kgrid
    kgrid = ndgrid_matSizeIn(padSize,1,'centerZero_ifftshift'); % pix^-1
    kgrid = cat(3, kgrid{1}, kgrid{2});
    kgrid = kgrid./effPix; % um^-1
    if gpuInd
        xgrid = gpuArray(single(xgrid));
        kgrid = gpuArray(single(kgrid));
        SLMmask = gpuArray(single(SLMmask));
    end

    %% load phi0
    phase0File = dir(fullfile( topPath, sprintf('Phase0_%d*.mat',round(wl*1e3)) ));
    
    if isempty(phase0File)
        warning('No SLM phase correction pattern detected.')
        phi0 = zeros(padSize);
    else
        assert(isscalar(phase0File), 'The phase correction file must be unique.')

        d = load( fullfile( phase0File(1).folder, phase0File(1).name) );
        phi0 = d.SLMcorrection;
        phi0 = phi0(:, 1:(end-invalidPix));
        assert( all(size(phi0) == SLMsize), 'phi0 size does not match' )
        phi0 = mpad(phi0, padSize);
        fprintf('SLM phase correction pattern found at %s (loaded %dx%d)\n', ...
            fullfile(phase0File(1).folder, phase0File(1).name), size(d.SLMcorrection,1), size(d.SLMcorrection,2));
    end
    if gpuInd
        phi0 = gpuArray(single(phi0));
    end

    %% def. NAmask
    krhosq = (kgrid(:,:,1) - rampVec(1)).^2 + (kgrid(:,:,2) - rampVec(2)).^2;
    NAmask = krhosq <= ( k_NA )^2;
    NAmask_safe = krhosq <= ( k_NA * NAsafetyFactor )^2;
    % figure, imagesc(fftshift(krhosq)), axis image
    % figure, imagesc(fftshift(NAmask)), axis image
    % figure, imagesc(fftshift(NAmask_safe)), axis image

    %% Ein0
    Ein0  = tag2Ein0(patTag, E0, SLMsize, effPix, k_NA, SIM_ramp_shift);
    E_input = Ein0;

    %%% Fine tuning of
    k_parax = getParaxialKernel(size(Ein0,1),size(Ein0,2));
    Ein0 = angularSpectrumPropagation(Ein0, AxialPropZ, k_parax);

    %% Etarget
    %%% padding  masking
    Ein0 = mpad(Ein0, padSize);

    %%% phase ramp
    phaseRamp = exp( 1i*2*pi*(xgrid(:,:,1)*rampVec(1)+xgrid(:,:,2)*rampVec(2)) );
    Etarget = Ein0 .* phaseRamp;
    Etarget = ifft2(fft2(Etarget) .* NAmask);

    %%% phase0
    % Etarget = Etarget .* exp(-1i*phi0);

    %%% normalization
    % nonZeroFraction = sum(Etarget ~=0,'all') / prod(SLMsize);
    % totalInt  = sum(abs(Etarget).^2,'all');
    targetInt = prod(SLMsize);
    % Etarget = Etarget * sqrt( targetInt/totalInt  );

    % maxTargetAmp = max(abs(Etarget),[],'all');
    % Etarget = Etarget./maxTargetAmp;
    % figure, imagesc(abs(Etarget).^2), axis image
    % figure, imagesc(fftshift(abs(fft2(Etarget)))); axis image

    %%% for corr loss
    Etarget = Etarget/norm(Etarget(:));

    %% iteration
    tic;
    % figure,
    % phi_iter = zeros( padSize );
    phi_iter = rand( padSize )*2*pi;

    OptMethod = 'NAG1';
    eta = 1;
    reg_param = reg_param_input; %0.01;
    cont_param = 0.00000;
    iterMax = 100;
    phi2SLMconst = SLM2pival/2/pi;


    %%% GPU variable
    lossLib = zeros(iterMax,2);
    if gpuInd
        phi_iter = gpuArray(single(phi_iter));
        lossLib  = gpuArray(single(lossLib));
    end

    %%% Determine Optimization Method
    switch OptMethod
        case 'NAG1'
            mMap  = zeros(size(phi_iter),'like', phi_iter);
            NAGcoeff = zeros(iterMax,1);
            t0 = 1;
            for ii = 1:iterMax
                t1 = (1+sqrt(1+4*t0^2))/2;
                NAGcoeff(ii) = (t0-1)/t1;
                t0=t1;
            end

        case 'NAG2'
            mMap = zeros(size(phi_iter),'like', phi_iter);
            beta1 = 0.9;

        case 'Adam'
            mMap = zeros(size(phi_iter),'like',phi_iter);
            vMap = zeros(size(phi_iter),'like',phi_iter);
            % vMap_hat = 0;
            beta1 = 0.9;
            beta2 = 0.999;
            epsilon = 10^-8;

    end

    for ii = 1:iterMax

        %%% rounding
        switch OptMethod
            case 'NAG2'
                phiSLM_iter = round(real(phi_iter - beta1*mMap) * phi2SLMconst) /phi2SLMconst;
            otherwise
                phiSLM_iter = round(real(phi_iter) * phi2SLMconst) /phi2SLMconst;
        end
        % phiSLM_iter = phi_iter;

        Eiter = exp(1i*(phiSLM_iter + phi0)) .* SLMmask;
        Fiter = ifft2(fft2(Eiter) .* NAmask_safe);
        Fiter_NA = ifft2(fft2(Eiter) .* NAmask);

        % Fiter = Eiter;

        %%% backward
        %%% correlation loss
        normFactor = norm(Fiter(:));
        % Fiter_norm = Fiter / normFactor;
        % Ecorr = conj(Etarget).*Fiter_norm;
        % Ecorrsum = sum(Ecorr,'all');
        % g = real(1i * ifft2(fft2(Ecorrsum.*Etarget*normFactor.* (1 - 1/2/normFactor^2.*abs(Fiter).^2)).*NAmask_safe) .* conj(Eiter));
        % lossLib(ii) = -abs(Ecorrsum).^2;

        %%% L2 loss
        Esub = Fiter - Etarget.*normFactor;
        g = real(-1i * ifft2(fft2(Esub - reg_param*Fiter_NA).*NAmask_safe) .* conj(Eiter));
        lossLib(ii,1) = abs(Etarget(:)' * Fiter(:) / normFactor).^2;
        % lossLib(ii,1) = abs(Etarget(:)' * Fiter_NA(:) / norm(Fiter_NA(:))).^2;

        %%% minimzie discontinuity
        if cont_param > 0
            phi_SLM_iter_mod = mod(phiSLM_iter*phi2SLMconst, SLM2pival);
            g_disc_x = phi_SLM_iter_mod - circshift(phi_SLM_iter_mod,[0,1]);
            g_disc_y = phi_SLM_iter_mod - circshift(phi_SLM_iter_mod,[1,0]);
            g = g + cont_param*( g_disc_x + g_disc_y );
        end
        % Esub = Eiter - Etarget;
        % g = real(-1i * ifft2(fft2(Esub)) .* conj(Eiter));

        lossLib(ii,2) = sum(abs(Fiter_NA).^2,'all')/targetInt;
        % lossLib(ii,2) = normFactor.^2/targetInt;


        %%% opt

        switch OptMethod
            case 'NAG1'
                mMap_prev = mMap;
                mMap = phi_iter - eta * g;
                phi_iter = mMap + NAGcoeff(ii).*(mMap - mMap_prev);

            case 'NAG2'
                mMap = beta1*mMap + eta*g;
                phi_iter = phi_iter - mMap;

            case 'Adam'
                mMap   = beta1*mMap + (1-beta1)*g;
                vMap   = beta2*vMap  + (1-beta2)*abs(g).^2;
                mMap_hat = mMap / (1-beta1^ii);
                vMap_hat = vMap / (1-beta2^ii);
                phi_iter  = phi_iter - eta * mMap_hat ./ (sqrt(vMap_hat) + epsilon);
        end
        % phi_iter = phi_iter - stepSize*g;

        %%% Figure plot during optimization
        % subplot(221), imagesc(mod(phiSLM_iter*phi2SLMconst,SLM2pival)),axis image; colorbar
        % subplot(222), contrastVis(abs(fftshift(fft2(Eiter))),99.99); axis image; colorbar
        % subplot(223), imagesc(abs(Fiter.^2)),axis image; colorbar
        % subplot(224),
        % yyaxis left,  plot(lossLib(1:ii,1)); colorbar;
        % title(sprintf('corr: %.3f\n Preserved power: %.3f %%',lossLib(ii,1),lossLib(ii,2)*100))
        % yyaxis right, plot(lossLib(1:ii,2)); colorbar;
        % drawnow;

        fprintf('%03d/%03d\n',ii,iterMax);

    end
    SLMphase = uint8(mod(phiSLM_iter*phi2SLMconst, SLM2pival));
    SLMphase = gather(mcrop(SLMphase, SLMsize));
    invalidArea = zeros(SLMsize(1), invalidPix, 'uint8');
    % invalidArea = mod(round(SLM2pival*rand(SLMsize(1), 8+aa)),SLM2pival); %

    SLMphase = [SLMphase , invalidArea];

    toc;
    figure(1232),
    subplot(221), imagesc(SLMphase),axis image; colorbar;
    title(sprintf('SLM pattern\n invalid pix: %d',invalidPix))

    subplot(222), contrastVis(log10(abs(fftshift(fft2(Eiter)))),99.9999); axis image; colorbar;
    title(sprintf('log(Fourier plane)\n NAsafetyFactor: %.2f',NAsafetyFactor))
    colormap(subplot(222),'turbo')

    subplot(223), imagesc(mcrop(abs(Fiter_NA).^2, SLMsize)),axis image; colorbar
    title(sprintf('image plane (itensity)'))

    subplot(224),
    yyaxis left,  plot(lossLib(1:ii,1));  ylabel('corr')
    yyaxis right, plot(lossLib(1:ii,2));  ylabel('power fraction')
    title(sprintf('corr: %.6f\n Preserved power: %.6f %%',lossLib(ii,1),lossLib(ii,2)*100))
    xlabel('iteration #')
    drawnow;

    %% fieldOut
    E_out = gather( Fiter_NA.*conj(phaseRamp) );
    
    %% save    
    save(saveStr,'SLMphase','E_out','effPix','-v7.3',"-nocompression");
    

end

% if exist(Eiter, "var")
%     E_FFT = fftshift(fft2(Eiter));
% else
%     E_FFT = 0;
% end


%% display
if displayTF
    hAx = fullScreenFigure(999, monitorInd);
    imagesc(SLMphase, 'Parent',hAx);
    clim([0, 255]); axis image off, colormap('gray');
    drawnow;
end


%% display (Psychtoolox)
% [windowPtr,rect] = Screen('OpenWindow',1);
% Screen('PutImage', windowPtr, outImg);
% Screen('Flip', windowPtr);
% Screen('CloseAll');

end




%% subfunctions

function Ein0 = tag2Ein0(patTag, Ein0, SLMsize, effPix, k_NA, SIM_ramp_shift)

%%% split the tag
patTagCell = split(patTag,'_');
patName   = patTagCell{1};
patParams = double(string(patTagCell(2:end)));

switch patName
    case 'normal'%%% planewave
        Ein0 = ones(SLMsize);

    case 'point' %%% point
        Ein0 = zeros(SLMsize);
        Ein0(1,1) = 1;
        Ein0 = ifftshift(Ein0);

    case 'quarter' %%% 1st quadrant
        Ein0 = padarray( ones(round(SLMsize/2)), SLMsize-round(SLMsize/2), 0, "post" );

    case 'edge' %%% SLM edge highlight
        linewidth = patParams(1);
        Ein0 = ~mpad(~mpad(true(linewidth),SLMsize-2*linewidth),SLMsize);

    case 'RayleighSpeckle' %%% Rayleigh speckle
        NAportion = patParams(1);
        Kradius = k_NA*NAportion; % um^-1, minval = 2*pix;
        Ein0 = randn(SLMsize) +  1i*randn(SLMsize);

        %%% NAcrop
        kgrid0 = ndgrid_matSizeIn(SLMsize,1,'centerZero_ifftshift'); % pix^-1
        kgrid0 = cat(3, kgrid0{1}, kgrid0{2});
        kgrid0 = kgrid0./effPix; % um^-1

        NAmask = (kgrid0(:,:,1)./Kradius).^2 + (kgrid0(:,:,2)./Kradius).^2 <= 1;
        Ein0 = ifft2(fft2(Ein0).*NAmask);

    case 'BesselSpeckle' %%% Bessel speckle
        NAportion = patParams(1);
        Kradius = k_NA*NAportion; % um^-1, minval = 2*pix;
        Kwidth  = 1./(SLMsize*effPix);
        kgrid0 = ndgrid_matSizeIn(SLMsize,1,'centerZero'); % pix^-1
        kgrid0 = cat(3, kgrid0{1}, kgrid0{2});
        kgrid0 = kgrid0./effPix; % um^-1
        FEin0 = (kgrid0(:,:,1)./(Kradius - Kwidth(1))).^2 + (kgrid0(:,:,2)./(Kradius - Kwidth(2))).^2 >= 1 ...
            & (kgrid0(:,:,1)./(Kradius + Kwidth(1))).^2 + (kgrid0(:,:,2)./(Kradius + Kwidth(2))).^2 <= 1;
        Ein0 = fftshift(ifft2(ifftshift(FEin0.*exp(1i*2*pi*rand(SLMsize)))));

    case 'FourierRing' %%% Bessel speckle
        NAportion = patParams(1);
        Kradius = k_NA*NAportion; % um^-1, minval = 2*pix;
        Kwidth  = 1./(SLMsize*effPix);
        kgrid0 = ndgrid_matSizeIn(SLMsize,1,'centerZero'); % pix^-1
        kgrid0 = cat(3, kgrid0{1}, kgrid0{2});
        kgrid0 = kgrid0./effPix; % um^-1
        FEin0 = (kgrid0(:,:,1)./(Kradius - Kwidth(1))).^2 + (kgrid0(:,:,2)./(Kradius - Kwidth(2))).^2 >= 1 ...
            & (kgrid0(:,:,1)./(Kradius + Kwidth(1))).^2 + (kgrid0(:,:,2)./(Kradius + Kwidth(2))).^2 <= 1;
        Ein0 = fftshift(ifft2(ifftshift(FEin0.*exp(1i*2*pi*zeros(SLMsize)))));

    case {'2DSIM','3DSIM'}
        % %% sine pattern
        NAportion = patParams(1);
        modAngle  = patParams(2); % deg
        addPhase = patParams(3); % deg

        kmodVec = k_NA * NAportion .* [sind(modAngle), cosd(modAngle)]; % um^-1
        xgrid0 = ndgrid_matSizeIn(SLMsize,0,'centerZero'); % pix
        xgrid0 = cat(3, xgrid0{1}, xgrid0{2});
        xgrid0 = xgrid0.*effPix; % um
        firstOrder = exp( 1i*2*pi*(xgrid0(:,:,1)*kmodVec(1)+xgrid0(:,:,2)*kmodVec(2)) + 1i*deg2rad(addPhase) );
        Ein0 = real((firstOrder + conj(firstOrder))/2);
        
        if strcmpi(patName, '3DSIM')
            Ein0 = Ein0+1; % nomal angle addition
        end

    case {'2DSIMincoh','2DSIMincohSimul'}
        % %% sine pattern
        NAportion = patParams(1);
        modAngle  = patParams(2); % deg
        addPhase = patParams(3); % deg

        % ... your code above ...
        kmodVec = k_NA * NAportion .* [sind(modAngle), cosd(modAngle)]; % [um^-1]
        xgrid0 = ndgrid_matSizeIn(SLMsize,0,'centerZero'); % [pix]
        xgrid0 = cat(3, xgrid0{1}, xgrid0{2});
        xgrid0 = xgrid0 .* effPix;                           % [um]

        phi = 2*pi*( xgrid0(:,:,1)*kmodVec(1) + xgrid0(:,:,2)*kmodVec(2) ) + deg2rad(addPhase);
        firstOrder = exp(1i*phi);
        Ein0 = real((firstOrder + conj(firstOrder))/2);                  

        kRamp = kmodVec * SIM_ramp_shift;                           
        phiRamp = 2*pi*( xgrid0(:,:,1)*kRamp(1) + xgrid0(:,:,2)*kRamp(2) );

        % Pattern
        Ein0 = Ein0 .* exp(1i*phiRamp);     

        k_parax = getParaxialKernel(size(Ein0,1),size(Ein0,2));
        if strcmpi(patName, '2DSIMincoh')
            Ein0 = angularSpectrumPropagation(Ein0, 0.3, k_parax);
        end

    otherwise
        assert(all(size(Ein0)==SLMsize), 'Ein0 size mismatch');
end

%%% show
% figure,
% subplot(121),imagesc(abs(Ein0)); axis image
% subplot(122), imagesc(fftshift(abs(fft2(Ein0)))); axis image
end