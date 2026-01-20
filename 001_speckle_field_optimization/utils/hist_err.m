% histogram Error
function [e, histoValue, x, targetValue] = hist_err(I, target, graph, varargin)
% When only the histogram is needed, input I only.
if nargin < 2
    target = zeros(length(I));
end

if length(size(I))==3
    I = reshape(I, [size(I,1)*size(I,2), size(I,3)]);
end
if length(size(target))==3
    target = reshape(target, [size(target,1)*size(target,2), size(target,3)]);
end

% Set bin size
bin_size = 0.05;

Color = 0;
LineWidth = 1;
XLim = 3;
ColorDiff = 0;
for setting = 1:2:nargin-3
    switch varargin{setting}
        case 'Color'
            Color = varargin{setting + 1};
        case 'LineWidth'
            LineWidth = varargin{setting + 1};
        case 'XLim'
            XLim = varargin{setting + 1};
        case 'ColorDiff'
            ColorDiff = varargin{setting + 1};
        case 'bin_size'
            bin_size = varargin{setting + 1};
    end
end

% Normalize intensities
I = I/mean(mean(I));
target = target/mean(mean(target));
N1 = size(I,1)*size(I,2);
N2 = size(target,1)*size(target,2);


% Number of bins
Maxbin = XLim; %15;
Num_of_bins = round(Maxbin / bin_size);

% Set histogram variables
histoValue = zeros(1,Num_of_bins);
targetValue = zeros(1,Num_of_bins);

% Count the histogram except the last point
I_count = fix(I/bin_size);
target_count = fix(target/bin_size);
for i = 1:Num_of_bins  %-1
    histoValue(i) = sum(sum(I_count == i-1));
    targetValue(i) = sum(sum(target_count == i-1));
end

% % The last point includes the remaining
% histoValue(Num_of_bins) = N1-sum(histoValue(1:end-1));
% targetValue(Num_of_bins) = N2-sum(targetValue(1:end-1));

% Normalize the histogram
histoValue = histoValue/N1/bin_size;
targetValue = targetValue/N2/bin_size;

histoValue = histoValue/sum(sum(histoValue));
targetValue = targetValue/sum(sum(targetValue));

% draw histogram if graph == 1
if nargin > 2
    if graph == 1   % including target
        % draw histogram
        figure(31), hold on;
        if Color == 0
            a = plot(0:bin_size:Maxbin-bin_size, histoValue/bin_size, '-','LineWidth',LineWidth);
            if ColorDiff == 0
                plot(0:bin_size:Maxbin-bin_size, targetValue/bin_size, '--', 'Color', a.Color,'LineWidth',LineWidth)
            else
                plot(0:bin_size:Maxbin-bin_size, targetValue/bin_size, '--', 'LineWidth',LineWidth)
            end
        else
            plot(0:bin_size:Maxbin-bin_size, histoValue/bin_size, '--','Color',Color,'LineWidth',LineWidth)
            plot(0:bin_size:Maxbin-bin_size, targetValue/bin_size, '-','Color',Color,'LineWidth',LineWidth)
        end
        legend('I','target'), FS, title('histogram')
    end
    if graph > 1    % No including target
        % draw histogram
        figure(31), hold on;
        if Color == 0
            plot(0:bin_size:Maxbin-bin_size, histoValue/bin_size,'LineWidth',LineWidth)
        else
            plot(0:bin_size:Maxbin-bin_size, histoValue/bin_size, '-','Color',Color,'LineWidth',LineWidth)
        end
        legend('I'), FS, title('histogram')
    end
end

e = sum(abs(histoValue-targetValue))/2;
x = 0:bin_size:Maxbin-bin_size;

histoValue = histoValue/bin_size;
targetValue = targetValue/bin_size;


end



