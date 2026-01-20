function FS(n)

if nargin == 0
    n = 14;
end

ax = gca;
ax.FontSize = n;
ax.FontName = 'Arial';
% ax.FontName = 'Times New Roman';

box on;

end