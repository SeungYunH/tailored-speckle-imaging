function [hAx,sizeOut] = fullScreenFigure(figInd, Nmonitor)

if nargin < 1
    figInd = [];
end

if nargin < 2
    Nmonitor = 1;
end


MP = get(0, 'MonitorPositions');
sizeOut = MP(Nmonitor, [4,3]); % [V, H]
% N = size(MP, 1);
% newPosition = MP(1,:);
% if size(MP, 1) == 1  % 'Single monitor -- do nothing'
% else
%     % Multiple monitors - shift to the Nth monitor.
%     newPosition(1) = newPosition(1) + MP(N,1);
% end

if isempty(figInd)
    fh = figure();
else
    fh = figure(figInd);
end
fh.set('Position', MP(Nmonitor, :), 'units', 'pixels');
fh.WindowState = 'fullscreen';
% fh.OuterPosition=[0 0 1 1];
set(fh, 'Toolbar', 'none', 'Menu', 'none');

hAx  = gca;
set(hAx,'Unit','normalized','Position',[0 0 1 1]);
end